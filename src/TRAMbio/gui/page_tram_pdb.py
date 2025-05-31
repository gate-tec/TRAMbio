import os
import re
import shutil
import time
from typing import Optional

from nicegui import ui, run

from TRAMbio.gui.io.local_file_picker import local_file_picker
from TRAMbio.gui.util.shutdown import shutdown
from TRAMbio.gui.viewer.config import init_header
from TRAMbio.gui.parameters.interactions import create_interaction_tabs
from TRAMbio.gui.viewer.viewer import Viewer, ViewerStyle
from TRAMbio.services import ParameterRegistry, WorkflowServiceRegistry, IOServiceRegistry
from TRAMbio.services.parameter import GeneralWorkflowParameter


step1_next = None
re_load_interactions = True

atom_viewer = None
interaction_viewer = None
component_viewer = None

pdb_path = None
out_dir = None
out_prefix = None


def run_tram_pdb(
        parameter_id: str,
        pdb_path: str,
        out_dir: str,
        out_prefix: str,
        origin_registry
) -> str:
    """Workflow integration from CLI"""
    if origin_registry is not None:
        # load registry to enable parallel CPU-bound execution
        import TRAMbio.services.parameter.registry as registry_space
        registry_space.ParameterRegistry.load_registry(origin_registry)

    ParameterRegistry.get_parameter_set(parameter_id).set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

    pdb_workflow_service = WorkflowServiceRegistry.PDB.single_service()
    pymol_workflow_service = WorkflowServiceRegistry.PYMOL.single_service()

    out_file = os.path.join(out_dir, f"{out_prefix}_components.xml")
    temp_file = os.path.join(out_dir, f".temp_{out_prefix}_components.xml")

    out_pdb_file = f"{out_prefix}_components.pdb"
    out_pdb_path = os.path.join(out_dir, out_pdb_file)
    out_pml_path = os.path.join(out_dir, f"{out_prefix}_components.pml")

    structure_generator = pdb_workflow_service.pdb_to_components(
        pdb_path=pdb_path,
        edge_data_file=None,  # TODO: Change later to store bonds
        parameter_id=parameter_id
    )

    pdb_workflow_service.run_pipeline_on_generator(
        generator=structure_generator,
        out_file=out_file,
        temp_file=temp_file,
        parameter_id=parameter_id
    )

    pymol_workflow_service.convert_to_pymol_files(
        in_pdb_path=pdb_path,
        in_xtc_path=None,
        in_xml_path=out_file,
        in_bond_path=None,
        out_pdb_path=out_pdb_path,
        out_pml_path=out_pml_path,
        out_prefix=out_prefix,
        rel_out_pdb_path=out_pdb_file,
        parameter_id=parameter_id
    )

    return out_pdb_path


def create_page_tram_pdb(cli_pdb_path: Optional[str]):

    # Download script
    def fetch_pdb_code(pdb_code: str):
        pdb_code = pdb_code.lower()
        if re.match(r"^[0-9a-z]{4}$", pdb_code) is not None:
            pdb_io_service = IOServiceRegistry.PDB.single_service()
            try:
                pdb_data = pdb_io_service.read(pdb_code)
                save_path = os.path.join(os.curdir, f"{pdb_code}.pdb")
                with open(save_path, "w") as out_file:
                    pdb_stream = pdb_data.to_pdb_stream()
                    pdb_stream.seek(0)
                    shutil.copyfileobj(pdb_stream, out_file, -1)  # noqa

                handle_input_file(save_path)
            except ValueError as e:
                ui.notify(f"[Error] {e.args}")
        else:
            ui.notify("Please provide valid PDB-code")

    # PDB load script
    def handle_input_file(file_path: str):
        if not os.path.exists(file_path) or not file_path.endswith(".pdb"):
            ui.notify("[Error] Please select a valid PDB file.")
            return
        global pdb_path, out_dir, out_prefix, re_load_interactions
        pdb_path = file_path
        out_dir = os.path.dirname(pdb_path)
        out_prefix = os.path.basename(pdb_path)[:-4]

        step1_next.props["disable"] = False
        step1_next.update()
        with open(pdb_path) as file:
            atom_viewer.add_models([file])
        re_load_interactions = True

    @ui.page('/tram-pdb')
    def page_tram_pdb():
        init_header()  # required for importing 3Dmol.js

        parameter_id = f"TRAMbio_GUI_PDB_{time.perf_counter_ns()}"

        #####################
        # Page Settings #####
        #####################

        ui.page_title('TRAMbio - PDB Analysis')

        ####################
        # Data Actions #####
        ####################

        async def pick_file() -> None:
            result = await local_file_picker(os.path.abspath(os.curdir), multiple=False, upper_limit=None)
            if result is not None and len(result) > 0:
                handle_input_file(result[0])

        #######################
        # Stepper Actions #####
        #######################

        async def action1_next() -> None:
            global re_load_interactions
            stepper.next()
            tab_panels.set_value("tabBonds")
            step2_next.props["disable"] = False
            step2_next.update()

            if re_load_interactions:
                with open(pdb_path) as file:
                    interaction_viewer.add_models([file])
                re_load_interactions = False

        async def action2_prev() -> None:
            stepper.previous()
            tab_panels.set_value("tabConfig")

        async def action2_next() -> None:
            # re-init viewer
            skeleton.set_visibility(True)
            skeleton.update()

            stepper.next()
            tab_panels.set_value("tabComp")
            tab_panels.update()

            # load parameters
            callback()

            # start calculation
            global pdb_path, out_dir, out_prefix, component_viewer
            out_pdb_path = await run.cpu_bound(
                run_tram_pdb,
                parameter_id, pdb_path, out_dir, out_prefix, ParameterRegistry.get_parameter_set(parameter_id)
            )

            skeleton.set_visibility(False)
            skeleton.update()

            with open(out_pdb_path) as file:
                component_viewer.add_models([file])

        ######################
        # Control Layout #####
        ######################

        global step1_next

        with ui.header(elevated=True).classes('w-full hidden'):
            with ui.tabs() as tabs:
                tab1 = ui.tab("tabConfig", label="Configuration")
                tab2 = ui.tab("tabBonds", label="Interactions").props("disable")
                tab3 = ui.tab("tabComp", label="Dilution Analysis").props("disable")

        with ui.left_drawer(top_corner=True, bottom_corner=True).style("border-right: 1px solid var(--q-primary)"):
            with ui.row():
                ui.link("TRAMbio", "/").classes('text-h4 text-black no-underline hover:underline')
                button_settings = ui.button(text="Settings", icon="settings")
                ui.button("Shutdown", on_click=shutdown).props('outline')

            with ui.stepper().props("vertical").classes("h-full").style("color: black;") as stepper:
                with ui.step("Configuration"):
                    ui.label("Load PDB and select atoms")
                    with ui.stepper_navigation():
                        step1_next = ui.button('Next', on_click=action1_next).props("disable")
                with ui.step("Interactions"):
                    ui.label("Specify interaction model")
                    with ui.stepper_navigation():
                        ui.button("Back", on_click=action2_prev).props('flat')
                        step2_next = ui.button("Next", on_click=action2_next).props("disable")
                with ui.step("Analysis"):
                    ui.label("Run Dilution Analysis")
                    with ui.stepper_navigation():
                        with ui.button("Start New", on_click=ui.navigate.reload).props('flat'):
                            ui.tooltip("Restart the analysis process at step 1")

        ###################
        # Base Layout #####
        ###################

        global atom_viewer
        global interaction_viewer
        global component_viewer

        with ui.tab_panels(tabs, value=tab1).classes('w-full') as tab_panels:

            ################################
            # Loading & Atom Selection #####
            ################################

            with ui.tab_panel(tab1):
                with ui.row().classes("w-full items-center"):
                    ui.button('Select file...', on_click=pick_file).props('icon=folder')

                    ui.space()

                    async def toggle_fetch(e) -> None:
                        fetch_button.props["disable"] = len(e.value) != 4
                        fetch_button.update()

                    pdb_code_input = ui.input(
                        label="PDB-code", placeholder="1L2Y",
                        validation={"Input too long": lambda value: len(value) <= 4},
                        on_change=toggle_fetch
                    )
                    fetch_button = ui.button("Download", on_click=lambda: fetch_pdb_code(pdb_code_input.value)).props("disable")


                with ui.splitter(value=30, reverse=True).classes('w-full') as splitter:
                    with splitter.before:
                        atom_viewer = Viewer(multi_model=False)

                        atom_viewer.create_viewer()

                    with splitter.after:
                        with ui.card() as viewer_container_2:
                            viewer_container_2.add_slot("viewer")
                            ui.label('Atom selection options').classes('text-h4')

            ##############################
            # Interaction definition #####
            ##############################

            with ui.tab_panel(tab2):
                with ui.splitter(value=50, reverse=True).classes("w-full") as splitter:
                    with splitter.before:
                        interaction_viewer = Viewer(multi_model=False)

                        interaction_viewer.create_viewer()

                        with ui.row().classes('items-center justify-between'):
                            ui.button("Recalculate Interactions", on_click=lambda x: ui.notify("This feature is not yet implemented."))

                    with splitter.after:
                        callback = create_interaction_tabs(parameter_id)

            #####################
            # Result Viewer #####
            #####################

            with ui.tab_panel(tab3):
                with ui.grid(columns='1fr auto', rows='1fr auto').classes('w-full h-full gap-5'):
                    with ui.element("div").classes("w-full relative").style("min-height: 650px;"):  # top left
                        skeleton = ui.skeleton().classes("absolute w-full h-full z-50")
                        skeleton.set_visibility(False)

                        component_viewer = Viewer(multi_model=True, style=ViewerStyle.COLORED_COMPONENTS)

                        component_viewer.create_viewer()

                    with ui.column().style("min-width: 150px;"):  # top right
                        pass

                    component_viewer.create_slider()  # bottom left

                    component_viewer.create_pagination()  # bottom right


        if cli_pdb_path is not None:
            handle_input_file(cli_pdb_path)
