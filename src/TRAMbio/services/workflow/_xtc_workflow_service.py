import time
from typing import Tuple, Optional, List, Union
import os

from concurrent.futures import ProcessPoolExecutor, Future
import numpy as np

import pandas as pd
from TRAMbio.services import StructureServiceRegistry, IOServiceRegistry, lock_registry, ParameterRegistry
from TRAMbio.services.parameter.registry import verbosity_from_parameter
from TRAMbio.services.workflow import WorkflowServiceRegistry, BaseWorkflowService, IXtcWorkflowService
from TRAMbio.services.workflow.util import init_mp_pool, create_graph_from_frame, apply_pebble_game_mp
from TRAMbio.services.parameter import XtcParameter, GeneralWorkflowParameter
from TRAMbio.util.structure_library.components import PebbleGameResult
from TRAMbio.util.structure_library.graph_struct import export_bond_frame
from TRAMbio.util.structure_library.generator import CustomGenerator, as_custom_generator

from tqdm import tqdm
from loguru import logger
from TRAMbio import set_log_level

_tqdm_logger = logger.bind(task="tqdm")

__all__ = []


def sigmoid(x):
    return 1 / (1 + np.exp(-x))


class MultiProcessingException(Exception):
    pass


class MockFuture:

    def __init__(self, fn, /, *args, **kwargs):
        self.__fn = fn
        self.__args = args
        self.__kwargs = kwargs
        self.__result = None

    def timed_run(self):
        t_start = time.perf_counter()
        self.__result = self.__fn(*self.__args, **self.__kwargs)
        return time.perf_counter() - t_start

    def done(self) -> bool:
        return self.__result is not None

    def result(self):
        if self.__result is None:
            raise MultiProcessingException
        return self.__result

    def cancel(self):
        pass


class XtcWorkflowService(BaseWorkflowService, IXtcWorkflowService):

    @property
    def name(self):
        return "XtcWorkflowService"

    @as_custom_generator(Tuple[str, List[PebbleGameResult]], pd.DataFrame)
    @lock_registry(kwargs_name="parameter_id")
    @verbosity_from_parameter("parameter_id", "verbose")
    def trajectory_to_components(
            self,
            xtc_path: str,
            pdb_path: str,
            edge_data_file: Optional[str] = None,
            cores: int = os.cpu_count() - 1,
            verbose: bool = False,
            parameter_id: str = ''
    ) -> CustomGenerator[Tuple[str, List[PebbleGameResult]], pd.DataFrame]:
        pdb_service = StructureServiceRegistry.PDB.single_service()
        parameter_registry = ParameterRegistry.get_parameter_set(parameter_id=parameter_id)
        xtc_io_service = IOServiceRegistry.XTC.query_service(parameter_registry(XtcParameter.MODULE.value))
        bond_io_service = IOServiceRegistry.BND.single_service()
        # deactivate verbosity for sub-tasks
        sub_parameter_id = parameter_id + "_sub"
        parameter_registry.clone(sub_parameter_id).set_parameter(GeneralWorkflowParameter.VERBOSE.value, False)

        logger.info('Building graph from initial frame...')

        num_frames, frame_generator = xtc_io_service.read(
            xtc_path=xtc_path, pdb_path=pdb_path, stride=parameter_registry(XtcParameter.STRIDE.value)
        )
        dynamic_scaling = parameter_registry(XtcParameter.DYNAMIC_SCALING.value)

        ###########################
        # Process first frame #####
        ###########################

        t_start = time.perf_counter()
        first_frame, raw_df = next(frame_generator)
        delta_frame = time.perf_counter() - t_start

        atom_df = pdb_service.export_atom_df(raw_df=raw_df, check_ids=True, parameter_id=sub_parameter_id)
        others_df = pdb_service.export_others_df(raw_df=raw_df, ter_only=False, parameter_id=sub_parameter_id)

        base_protein_graph = pdb_service.create_graph_struct(atom_df=atom_df, others_df=others_df,
                                                             parameter_id=sub_parameter_id)

        protein_graph = pdb_service.copy_graph_for_frame(atom_df=atom_df, others_df=others_df,
                                                         protein_graph=base_protein_graph, parameter_id=sub_parameter_id)

        # insert interactions into first frame
        pdb_service.apply_non_covalent_interactions(protein_graph=protein_graph, parameter_id=sub_parameter_id)

        # Tag frame number
        protein_graph.graphs['atom'].graph['frame'] = first_frame

        columns = []
        if edge_data_file is not None:
            # Store complete edge data for first frame (including constant base edges)
            old_edge_data = export_bond_frame(protein_graph.graphs['full'], include_base_edges=True)
            old_edge_data['_merge'] = 'right_only'
            old_edge_data['frame'] = first_frame
            columns = old_edge_data.columns
            # save edge data in git style tsv
            bond_io_service.store_bonds(edge_data_file, old_edge_data, mode='w')

        # Export edge data for first frame as reference
        old_edge_data = export_bond_frame(protein_graph.graphs['full'], include_base_edges=False) \
            .loc[:, ['node1', 'node2', 'type']]

        # run first frame
        logger.info('Calculating stability for initial frame...')
        mock_future = MockFuture(apply_pebble_game_mp, (first_frame, protein_graph.graphs,))
        delta_pebble = mock_future.timed_run()

        ######################
        # Process frames #####
        ######################

        available_threads = max(min(os.cpu_count() - 1, cores), 3)
        if dynamic_scaling:
            delta_scale = max(min(delta_pebble / delta_frame, 5), 0.2)  # clamp to interval [0.2, 5]
            # transform to [-4,4]
            if delta_scale < 1:
                delta_scale = (-1 / delta_scale) + 1
            else:
                delta_scale -= 1
            scaling_factor = sigmoid(delta_scale)

            base_workers = available_threads // 3
            max_pebble_workers = base_workers + int(scaling_factor * (available_threads - 2 * base_workers))
            max_frame_workers = available_threads - max_pebble_workers
        else:
            max_frame_workers = max(min(num_frames, (available_threads // 3) * 2), 2)  # apr. 2/3 of cores
            max_pebble_workers = max(available_threads - max_frame_workers, 1)  # apr. 1/3 of cores

        pre_fetch = 2  # factor of tasks to preload into working context in order to reduce idle time
        max_running_1 = max_frame_workers + int(pre_fetch * max_frame_workers)
        max_running_2 = max_pebble_workers + int(pre_fetch * max_pebble_workers)
        max_movements = min(max(max_frame_workers // 2, 1), 5)  # maximum number of intermediate movements between queue updates

        logger.info(f"Multi-Processing remaining {num_frames - 1} frames.")
        logger.debug(f"Using {max_frame_workers} cores for queue 1 and {max_pebble_workers} cores for queue 2.")

        # Two-pool multiprocessing with:
        # - Pool 1: graph construction from PDB data (every frame)
        # - Intermediate (local) check for each graph (discard if no change to previous frame)
        # - Pool 2: pebble game run on graph (potentially every frame)
        sub_registry = ParameterRegistry.get_parameter_set(parameter_id=sub_parameter_id)
        try:
            with ProcessPoolExecutor(
                    max_workers=max_frame_workers,
                    initializer=init_mp_pool,
                    initargs=(base_protein_graph, sub_registry, sub_parameter_id, 'CRITICAL')
            ) as executor_1, \
            ProcessPoolExecutor(
                max_workers=max_pebble_workers,
            ) as executor_2, \
            tqdm(total=num_frames, desc='Reading frames', disable=not verbose, position=0) as progress_1, \
            tqdm(total=num_frames, desc='Stability tests', disable=not verbose, position=1) as progress_2:

                set_log_level('CRITICAL')

                # initial work load
                work_load_1: List[Future] = [
                    executor_1.submit(create_graph_from_frame, (next(frame_generator), sub_parameter_id,))
                    for _ in range(min(max_running_1, num_frames - 1))
                ]
                work_load_2: List[Union[MockFuture, Future]] = [
                    mock_future  # inject finished run for first frame
                ]
                progress_1.update()  # first frame

                running_1 = len(work_load_1)
                remaining_1 = num_frames - running_1 - 1  # first frame already processed
                running_2 = 1
                remaining_2 = num_frames  # only decrement on FINISHED tasks for pool 2

                new_frames = 1
                old_frames = 0

                try:
                    while remaining_2 > 0:

                        # free-up queue 2, if possible
                        yield_future = None
                        if running_2 > 0 and work_load_2[0].done():  # non-blocking
                            yield_future = work_load_2.pop(0)
                            progress_2.update()
                            running_2 -= 1
                            remaining_2 -= 1

                        # update queue 2
                        movements = 0
                        while (running_1 > 0 and movements < max_movements and
                                running_2 < max_running_2 and  # evaluate space in queue 2
                                work_load_1[0].done()  # non-blocking call
                        ):
                            future = work_load_1.pop(0)
                            progress_1.update()
                            try:
                                new_frame, new_protein_graphs = future.result()  # blocking call
                                running_1 -= 1

                                # check difference
                                new_frame_graphs = new_protein_graphs.graphs
                                # Export and compare edge data to previous frame
                                new_edge_data = export_bond_frame(new_frame_graphs['full'], include_base_edges=False)

                                difference = old_edge_data \
                                    .merge(new_edge_data, indicator=True, how='outer', on=['node1', 'node2', 'type'],
                                           suffixes=(False, False)) \
                                    .loc[lambda x: x["_merge"] != 'both'].reset_index(drop=True)

                                if len(difference) > 0:
                                    if edge_data_file is not None:
                                        # save edges in data file (only differences to previous frame)
                                        difference["frame"] = new_frame
                                        bond_io_service.store_bonds(edge_data_file, difference[columns], mode='a')
                                    old_edge_data = new_edge_data.loc[:, ['node1', 'node2', 'type']]

                                    work_load_2.append(executor_2.submit(
                                        apply_pebble_game_mp, (new_frame, new_frame_graphs,)
                                    ))
                                    running_2 += 1
                                    # mark new frame
                                    new_frames += 1
                                else:
                                    # discarded frame
                                    old_frames += 1
                                    # advance progress for queue 2
                                    progress_2.update()
                                    remaining_2 -= 1

                                progress_1.set_postfix({"New": new_frames, "Old": old_frames})

                            finally:
                                future.cancel()

                        # update queue 1
                        while remaining_1 > 0 and running_1 < max_running_1:
                            # append next frame to queue
                            try:
                                work_load_1.append(executor_1.submit(
                                    create_graph_from_frame, (next(frame_generator), sub_parameter_id,)
                                ))
                                remaining_1 -= 1
                                running_1 += 1
                            except StopIteration as e:
                                raise MultiProcessingException() from e

                        # yield next result, if possible
                        if yield_future is not None:
                            try:
                                frame, components = yield_future.result()
                                yield str(frame), components
                            finally:
                                yield_future.cancel()
                finally:
                    for future in work_load_1:
                        future.cancel()
                    for future in work_load_2:
                        future.cancel()

        except Exception as e:
            raise MultiProcessingException() from e

        return base_protein_graph.hydrogen_mapping


WorkflowServiceRegistry.XTC.register_service(XtcWorkflowService())
