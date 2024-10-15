import os

from TRAMbio.services import IOServiceRegistry as IOServiceRegistry, ParameterRegistry as ParameterRegistry, \
    lock_registry
from TRAMbio.services.parameter.registry import verbosity_from_parameter
from TRAMbio.services.workflow import WorkflowServiceRegistry, IPebbleGameWorkflowService
from TRAMbio.pebble_game.pebble_game_mp import run_component_pebble_game_mp
from TRAMbio.pebble_game.pebble_game import get_evaluation_from_results
from TRAMbio.services.parameter import PebbleGameParameter


__all__ = []


class PebbleGameWorkflowService(IPebbleGameWorkflowService):

    @property
    def name(self):
        return "PebbleGameWorkflowService"

    @lock_registry(kwargs_name="parameter_id")
    @verbosity_from_parameter(parameter_name="parameter_id", verbose_name="verbose")
    def analyze_graph(
            self,
            graph_ml_path: str,
            out_path: str,
            verbose: bool = False,
            parameter_id: str = ''
    ) -> None:
        xml_io_service = IOServiceRegistry.XML.single_service()
        parameter_registry = ParameterRegistry.get_parameter_set(parameter_id=parameter_id)

        graph = xml_io_service.read_graphml(graphml_path=graph_ml_path)

        param_k = parameter_registry.get_parameter(PebbleGameParameter.K.value)
        param_l = parameter_registry.get_parameter(PebbleGameParameter.L.value)
        num_threads = min(max(parameter_registry.get_parameter(PebbleGameParameter.THREADS.value), 1), os.cpu_count())

        # parameter check
        if param_k < 1:
            raise ValueError  # TODO
        if not (0 <= param_l < 2 * param_k):
            raise ValueError  # TODO

        # evaluate graph
        rho, pebble_excess, components = run_component_pebble_game_mp(
            graph=graph, k=param_k, l=param_l, threads=num_threads, verbose=verbose
        )
        category = get_evaluation_from_results(rho=rho, pebble_excess=pebble_excess)

        xml_io_service.write_pebble_game_results(
            xml_out_path=out_path,
            category=category,
            components=components,
            parameter_id=parameter_id
        )


WorkflowServiceRegistry.PEBBLE.register_service(PebbleGameWorkflowService())
