from nicegui import ui, events, binding

from TRAMbio.services.parameter import HydrogenBondParameter, HydrophobicInteractionParameter, \
    DisulphideBridgeParameter, AromaticInteractionParameter, CationPiInteractionParameter, \
    PdbEntryInteractionParameter
from TRAMbio.services import ParameterRegistry


class AndGate:
    active_positive = binding.BindableProperty()
    active_negative = binding.BindableProperty()

    def __init__(self, input_1: bool, input_2: bool):
        self.input_1: bool = input_1
        self.input_2: bool = input_2
        self.active_positive: bool = self.input_1 and self.input_2
        self.active_negative: bool = self.input_1 and not self.input_2

    def update_input_1(self, e: events.ValueChangeEventArguments):
        self.input_1 = e.value
        self.update()

    def update_input_2(self, e: events.ValueChangeEventArguments):
        self.input_2 = e.value
        self.update()

    def update(self):
        self.active_positive = self.input_1 and self.input_2
        self.active_negative = self.input_1 and not self.input_2


def create_interaction_tabs(parameter_id: str):
    parameter_registry = ParameterRegistry.get_parameter_set(parameter_id)

    with ui.tabs().classes("w-full").props("dense shrink outside-arrows") as interaction_tabs:
        tab_hbond = ui.tab("tab_hbond", label="Hydrogen bonds").style('white-space: pre-wrap')
        ui.separator().props('vertical')
        tab_hydrophobic = ui.tab("tab_hydro", label="Hydrophobic interactions").style('white-space: pre-wrap')
        ui.separator().props('vertical')
        tab_disulphide = ui.tab("tab_disulphide", label="Disulphide Bridges").style('white-space: pre-wrap')
        ui.separator().props('vertical')
        tab_aromatic = ui.tab("tab_aromatic", label="Aromatic Interactions").style('white-space: pre-wrap')
        ui.separator().props('vertical')
        tab_cation = ui.tab("tab_cation", label="Cation-\u03C0 Interactions").style('white-space: pre-wrap')
        ui.separator().props('vertical')
        tab_pdb = ui.tab("tab_pdb", label="PDB Entries")
        # TODO: Add option for custom interactions (e.g., CSV upload)

    with ui.tab_panels(interaction_tabs, value=tab_hbond).classes('w-full'):

        ######################
        # HYDROGEN BONDS #####
        ######################

        with ui.tab_panel(tab_hbond):
            hbond_extra = AndGate(
                input_1=parameter_registry(HydrogenBondParameter.INCLUDE.value),
                input_2=parameter_registry(HydrogenBondParameter.STRONG_ENERGY_THRESHOLD.value) < parameter_registry(HydrogenBondParameter.ENERGY_THRESHOLD.value)
            )

            ui.separator().classes("w-full")

            # minimal length and bar-count
            with ui.row().classes("w-full"):
                ui.checkbox("Include", value=parameter_registry(HydrogenBondParameter.INCLUDE.value), on_change=hbond_extra.update_input_1) \
                    .bind_value(parameter_registry, HydrogenBondParameter.INCLUDE.value)

                ui.space()

                with ui.column():
                    with ui.number(label="Minimum Distance", min=0, value=parameter_registry(HydrogenBondParameter.MINIMUM_LENGTH.value), precision=3, step=0.1, suffix="\u212B") \
                            .classes("w-full") \
                            .bind_enabled_from(parameter_registry, HydrogenBondParameter.INCLUDE.value) \
                            .bind_value(parameter_registry, HydrogenBondParameter.MINIMUM_LENGTH.value):
                        ui.tooltip("Minimum length cut-off for hydrogen bonds")

                    with ui.row().classes("items-center"):
                        with ui.label("Bar count:"):
                            ui.tooltip("The number of bars to model (strong hydrogen bonds)")

                        ui.toggle({x: x for x in range(1, 7)}, value=parameter_registry(HydrogenBondParameter.BAR_COUNT.value)) \
                            .bind_enabled_from(parameter_registry, HydrogenBondParameter.INCLUDE.value) \
                            .bind_value(parameter_registry, HydrogenBondParameter.BAR_COUNT.value)

            ui.separator().classes("w-full")

            # energy threshold

            ui.label("Energy Threshold").classes("text-lg")

            async def synchronize_sliders(e: events.ValueChangeEventArguments):
                strong_slider.props["inner-max"] = e.value
                strong_slider.update()

            with ui.row().classes("w-full items-center"):
                with ui.column().classes("grow"):
                    ui.space()
                    ui.slider(min=-10, max=10, step=0.01, value=parameter_registry(HydrogenBondParameter.ENERGY_THRESHOLD.value), on_change=synchronize_sliders) \
                        .props('label label-always') \
                        .bind_enabled_from(parameter_registry, HydrogenBondParameter.INCLUDE.value) \
                        .bind_value(parameter_registry, HydrogenBondParameter.ENERGY_THRESHOLD.value)
                ui.number(label='Threshold', value=parameter_registry(HydrogenBondParameter.ENERGY_THRESHOLD.value)) \
                    .bind_enabled_from(parameter_registry, HydrogenBondParameter.INCLUDE.value) \
                    .bind_value(parameter_registry, HydrogenBondParameter.ENERGY_THRESHOLD.value)

            ui.separator().classes("w-full")

            # strong energy threshold

            with ui.switch("Strength Threshold", value=hbond_extra.active_positive, on_change=hbond_extra.update_input_2) \
                    .bind_enabled_from(parameter_registry, HydrogenBondParameter.INCLUDE.value) \
                    .bind_value(hbond_extra, "strong_energy"):
                ui.tooltip("A threshold separating strong and weak hydrogen bonds")

            with ui.row().classes("w-full items-center"):
                with ui.column().style("flex-grow: 1;"):
                    ui.space()
                    strong_slider = ui.slider(min=-10, max=10, step=0.01, value=parameter_registry(HydrogenBondParameter.STRONG_ENERGY_THRESHOLD.value)) \
                        .props(f'label label-always inner-max={parameter_registry(HydrogenBondParameter.ENERGY_THRESHOLD.value)}') \
                        .bind_enabled_from(hbond_extra, "active_positive") \
                        .bind_value(parameter_registry, HydrogenBondParameter.STRONG_ENERGY_THRESHOLD.value)
                ui.number(label='Strength Threshold', value=parameter_registry(HydrogenBondParameter.STRONG_ENERGY_THRESHOLD.value)) \
                    .bind_enabled_from(hbond_extra, "active_positive") \
                    .bind_value(parameter_registry, HydrogenBondParameter.STRONG_ENERGY_THRESHOLD.value)

        ################################
        # HYDROPHOBIC INTERACTIONS #####
        ################################

        with ui.tab_panel(tab_hydrophobic):
            hydrophobic_extra = AndGate(
                input_1=parameter_registry(HydrophobicInteractionParameter.INCLUDE.value),
                input_2=parameter_registry(HydrophobicInteractionParameter.POTENTIAL.value)
            )

            ui.separator().classes("w-full")

            # minimal length and bar-count

            with ui.row().classes("w-full"):
                ui.checkbox("Include", value=parameter_registry(HydrophobicInteractionParameter.INCLUDE.value), on_change=hydrophobic_extra.update_input_1) \
                    .bind_value(parameter_registry, HydrophobicInteractionParameter.INCLUDE.value)

                ui.space()

                with ui.column():
                    with ui.switch("Minimal length") \
                            .bind_enabled_from(parameter_registry, HydrophobicInteractionParameter.INCLUDE.value) \
                            .bind_value(parameter_registry, HydrophobicInteractionParameter.MINIMAL_LENGTH.value):
                        ui.tooltip("Only use the hydrophobic interactions with minimal length")

                    with ui.row().classes("items-center"):
                        with ui.label("Bar count:"):
                            ui.tooltip("The number of bars to model hydrophobic interactions")

                        ui.toggle({x: x for x in range(1, 6)},
                                  value=parameter_registry(HydrophobicInteractionParameter.BAR_COUNT.value)) \
                            .bind_enabled_from(parameter_registry, HydrophobicInteractionParameter.INCLUDE.value) \
                            .bind_value(parameter_registry, HydrophobicInteractionParameter.BAR_COUNT.value)

            ui.separator().classes("w-full")

            with ui.row().classes("w-full items-center justify-around"):
                with ui.radio(
                        options={False: "Surface Interaction", True: "Potential Interaction"},
                        value=hydrophobic_extra.active_positive,
                        on_change=hydrophobic_extra.update_input_2
                    ) \
                        .props("inline") \
                        .bind_enabled_from(parameter_registry, HydrophobicInteractionParameter.INCLUDE.value) \
                                .bind_value(parameter_registry, HydrophobicInteractionParameter.POTENTIAL.value):
                    ui.tooltip("Select the type of interaction calculation")

            ui.separator().classes("w-full")

            ui.label("Surface Interaction").classes("text-lg")

            with ui.number(label="Cut-off Distance", min=0,
                           value=parameter_registry(HydrophobicInteractionParameter.SURFACE_CUTOFF_DISTANCE.value),
                           precision=3, step=0.1, suffix="\u212B") \
                    .bind_enabled_from(hydrophobic_extra, "active_negative") \
                    .bind_value(parameter_registry, HydrophobicInteractionParameter.SURFACE_CUTOFF_DISTANCE.value):
                ui.tooltip("Maximum distance cut-off between the van-der-Waals surfaces for a hydrophobic interaction")

            ui.separator().classes("w-full")

            ui.label("Potential Interaction").classes("text-lg")

            with ui.number(label="Cut-off Distance", min=0,
                           value=parameter_registry(HydrophobicInteractionParameter.POTENTIAL_CUTOFF_DISTANCE.value),
                           precision=3, step=0.1, suffix="\u212B") \
                    .bind_enabled_from(hydrophobic_extra, "active_positive") \
                    .bind_value(parameter_registry, HydrophobicInteractionParameter.POTENTIAL_CUTOFF_DISTANCE.value):
                ui.tooltip("Maximum length cut-off for hydrophobic interactions")

            with ui.row().classes("w-full items-center justify-between"):
                with ui.number(label="Scale 1-4", min=0,
                               value=parameter_registry(
                                   HydrophobicInteractionParameter.SCALE_14.value),
                               precision=3, step=0.1) \
                        .bind_enabled_from(hydrophobic_extra, "active_positive") \
                        .bind_value(parameter_registry,
                                    HydrophobicInteractionParameter.SCALE_14.value):
                    ui.tooltip("Potential multiplier for 1-4 covalent neighbors")

                with ui.number(label="Scale 1-5", min=0,
                               value=parameter_registry(
                                   HydrophobicInteractionParameter.SCALE_15.value),
                               precision=3, step=0.1) \
                        .bind_enabled_from(hydrophobic_extra, "active_positive") \
                        .bind_value(parameter_registry,
                                    HydrophobicInteractionParameter.SCALE_15.value):
                    ui.tooltip("Potential multiplier for 1-5 covalent neighbors")

                with ui.number(label="Scale Unbounded", min=0,
                               value=parameter_registry(
                                   HydrophobicInteractionParameter.SCALE_UNBOUNDED.value),
                               precision=3, step=0.1) \
                        .bind_enabled_from(hydrophobic_extra, "active_positive") \
                        .bind_value(parameter_registry,
                                    HydrophobicInteractionParameter.SCALE_UNBOUNDED.value):
                    ui.tooltip("Potential multiplier for 1-6 covalent neighbors and further")

            with ui.row().classes("w-full items-center"):
                with ui.column().style("flex-grow: 1;"):
                    ui.space()
                    ui.slider(min=-10, max=10, step=0.01, value=parameter_registry(HydrophobicInteractionParameter.ENERGY_THRESHOLD.value)) \
                        .props(f'label label-always') \
                        .bind_enabled_from(hydrophobic_extra, "active_positive") \
                        .bind_value(parameter_registry, HydrophobicInteractionParameter.ENERGY_THRESHOLD.value)
                ui.number(label='Energy Threshold', value=parameter_registry(HydrophobicInteractionParameter.ENERGY_THRESHOLD.value)) \
                    .bind_enabled_from(hydrophobic_extra, "active_positive") \
                    .bind_value(parameter_registry, HydrophobicInteractionParameter.ENERGY_THRESHOLD.value)

        ##########################
        # DISULPHIDE BRIDGES #####
        ##########################

        with ui.tab_panel(tab_disulphide):

            ui.separator().classes("w-full")

            with ui.row().classes("w-full"):
                ui.checkbox("Include", value=parameter_registry(DisulphideBridgeParameter.INCLUDE.value), on_change=hydrophobic_extra.update_input_1) \
                    .bind_value(parameter_registry, DisulphideBridgeParameter.INCLUDE.value)

                ui.space()

                with ui.column():
                    with ui.number(label="Cut-off Distance", min=0,
                                   value=parameter_registry(HydrogenBondParameter.MINIMUM_LENGTH.value), precision=3,
                                   step=0.1, suffix="\u212B") \
                            .classes("w-full") \
                            .bind_enabled_from(parameter_registry, HydrogenBondParameter.INCLUDE.value) \
                            .bind_value(parameter_registry, HydrogenBondParameter.MINIMUM_LENGTH.value):
                        ui.tooltip("Maximum cut-off distance for disulphide bridges")

        #############################
        # AROMATIC INTERACTIONS #####
        #############################

        with ui.tab_panel(tab_aromatic):
            ui.separator().classes("w-full")

            with ui.row().classes("w-full"):
                ui.checkbox("Include", value=parameter_registry(AromaticInteractionParameter.INCLUDE.value),
                            on_change=hydrophobic_extra.update_input_1) \
                    .bind_value(parameter_registry, AromaticInteractionParameter.INCLUDE.value)


                ui.space()

                with ui.column():
                    with ui.number(label="Angle variance", min=0,
                                   value=parameter_registry(AromaticInteractionParameter.ANGLE_VARIANCE.value), precision=3,
                                   step=0.1, suffix="°") \
                            .classes("w-full") \
                            .bind_enabled_from(parameter_registry, AromaticInteractionParameter.INCLUDE.value) \
                            .bind_value(parameter_registry, AromaticInteractionParameter.ANGLE_VARIANCE.value):
                        ui.tooltip("Allowed angle deviations for aromatic interactions")

                    with ui.row().classes("items-center"):
                        with ui.label("Bar count:"):
                            ui.tooltip("The number of bars to model aromatic interactions")

                        ui.toggle({x: x for x in range(1, 6)},
                                  value=parameter_registry(AromaticInteractionParameter.BAR_COUNT.value)) \
                            .bind_enabled_from(parameter_registry, AromaticInteractionParameter.INCLUDE.value) \
                            .bind_value(parameter_registry, AromaticInteractionParameter.BAR_COUNT.value)

            ui.separator()

            with ui.row().classes("w-full"):
                with ui.number(label="\u03C0-Stacking cut-off", min=0,
                               value=parameter_registry(AromaticInteractionParameter.CUTOFF_DISTANCE_PI.value), precision=3,
                               step=0.1, suffix="\u212B") \
                        .bind_enabled_from(parameter_registry, AromaticInteractionParameter.INCLUDE.value) \
                        .bind_value(parameter_registry, AromaticInteractionParameter.CUTOFF_DISTANCE_PI.value):
                    ui.tooltip("Distance cut-off for \u03C0-stacking")

                with ui.number(label="T-Stacking cut-off", min=0,
                               value=parameter_registry(AromaticInteractionParameter.CUTOFF_DISTANCE_T.value), precision=3,
                               step=0.1, suffix="\u212B") \
                        .bind_enabled_from(parameter_registry, AromaticInteractionParameter.INCLUDE.value) \
                        .bind_value(parameter_registry, AromaticInteractionParameter.CUTOFF_DISTANCE_T.value):
                    ui.tooltip("Distance cut-off for T-stacking")

        ##############################
        # CATION-PI INTERACTIONS #####
        ##############################

        with ui.tab_panel(tab_cation):
            ui.separator().classes("w-full")

            with ui.row().classes("w-full"):
                ui.checkbox("Include", value=parameter_registry(CationPiInteractionParameter.INCLUDE.value),
                            on_change=hydrophobic_extra.update_input_1) \
                    .bind_value(parameter_registry, CationPiInteractionParameter.INCLUDE.value)

                ui.space()

                with ui.column():
                    with ui.number(label="Cut-off Distance", min=0,
                                   value=parameter_registry(CationPiInteractionParameter.CUTOFF_DISTANCE.value), precision=3,
                                   step=0.1, suffix="\u212B") \
                            .classes("w-full") \
                            .bind_enabled_from(parameter_registry, CationPiInteractionParameter.INCLUDE.value) \
                            .bind_value(parameter_registry, CationPiInteractionParameter.CUTOFF_DISTANCE.value):
                        ui.tooltip("Maximum cut-off distance for cation-\u03C0 interactions")

                    with ui.row().classes("items-center"):
                        with ui.label("Bar count:"):
                            ui.tooltip("The number of bars to model cation-\u03C0 interactions")

                        ui.toggle({x: x for x in range(1, 6)},
                                  value=parameter_registry(CationPiInteractionParameter.BAR_COUNT.value)) \
                            .bind_enabled_from(parameter_registry, CationPiInteractionParameter.INCLUDE.value) \
                            .bind_value(parameter_registry, CationPiInteractionParameter.BAR_COUNT.value)

        ###################
        # PDB ENTRIES #####
        ###################

        with ui.tab_panel(tab_pdb):
            ui.separator().classes("w-full")

            ui.checkbox("Include SS-BOND records", value=parameter_registry(PdbEntryInteractionParameter.SSBOND_INCLUDE.value),
                        on_change=hydrophobic_extra.update_input_1) \
                .bind_value(parameter_registry, PdbEntryInteractionParameter.SSBOND_INCLUDE.value)

            ui.checkbox("Include LINK records", value=parameter_registry(PdbEntryInteractionParameter.LINK_INCLUDE.value),
                        on_change=hydrophobic_extra.update_input_1) \
                .bind_value(parameter_registry, PdbEntryInteractionParameter.LINK_INCLUDE.value)

            ui.checkbox("Include CONECT records", value=parameter_registry(PdbEntryInteractionParameter.CONECT_INCLUDE.value),
                        on_change=hydrophobic_extra.update_input_1) \
                .bind_value(parameter_registry, PdbEntryInteractionParameter.CONECT_INCLUDE.value)

    def callback():
        # format hydrogen bonds
        if not hbond_extra.active_positive:
            parameter_registry.set_parameter(
                HydrogenBondParameter.STRONG_ENERGY_THRESHOLD.value,
                parameter_registry(HydrogenBondParameter.ENERGY_THRESHOLD.value) + 1.0
            )
        pass

    return callback