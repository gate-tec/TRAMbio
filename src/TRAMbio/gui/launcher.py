import os
import argparse

from nicegui import ui

from TRAMbio import set_log_level
from TRAMbio.gui.page_tram_pdb import create_page_tram_pdb
from TRAMbio.gui.util.shutdown import shutdown
from TRAMbio.gui.parameters.registry import create_registry_with_bindable_properties


@ui.page('/')
def home_page():
    set_log_level("INFO")

    # make parameters into bindable registry attributes
    create_registry_with_bindable_properties()

    #####################
    # Page Settings #####
    #####################

    ui.page_title('TRAMbio - Topological Rigidity Analysis in Molecular Biology')

    ####################
    # Page Content #####
    ####################

    with ui.left_drawer(top_corner=True, bottom_corner=True).style("border-right: 1px solid var(--q-primary)"):
        with ui.row():
            ui.label("TRAMbio").classes('text-h4')
            button_settings = ui.button(text="Settings", icon="settings")
            ui.button("Shutdown", on_click=shutdown).props('outline')

    with ui.row().classes("w-full h-full"):
        with ui.card().classes("items-center hover:cursor-pointer").props('flat bordered').on("click", lambda: ui.navigate.to("/tram-pdb")):
            ui.icon("biotech").classes('text-3xl')
            ui.label("PDB Analysis").classes('text-3xl')
        with ui.card().classes("items-center hover:cursor-not-allowed").props('flat bordered disable'):
            ui.icon("timeline").classes('text-3xl')
            ui.label("Trajectory Analysis").classes('text-3xl')


#####################
# Run Functions #####
#####################

def launch(reload=False, pdb_path=None):

    # create sub-pages
    create_page_tram_pdb(pdb_path)

    # TODO: insert tram-xtc page

    ui.run(reload=reload, show_welcome_message=False, show=True)  # reload False is required for production environment


def main():
    parser = argparse.ArgumentParser(prog="TRAMbio",
                                     description="GUI for TRAMbio workflows",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-p', '--pdb',
                        type=str, metavar='PDB_FILE', dest="pdb", required=False,
                        help="Protein input file in PDB v3.3 format.")
    parser.add_argument('-x', '--xtc',
                        type=str, metavar='XTC_FILE', dest="xtc", required=False,
                        help="Trajectory file in XTC format.")

    args = parser.parse_args()

    launch(reload=False, pdb_path=args.pdb)


if __name__ in {"__main__", "__mp_main__"}:
    launch(reload=True)  # development entry-point
