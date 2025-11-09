import sys
from enum import Enum
from typing import List, Union


class MissingDependencyError(Exception):

    def __init__(self, module, dependency: Union[str, List[str]]):
        dependency_string = dependency
        ending = 'y'
        if isinstance(dependency, list):
            dependency_string = ", ".join(dependency)
            if len(dependency) > 1:
                ending = 'ies'
        message = f"{sys.modules['TRAMbio'].__name__} submodule \"{module}\" requires missing dependenc{ending}: {dependency_string}"
        super().__init__(message)
        self.message = message


class DependencyWithSafeTest(Enum):
    MDAnalysis = "MDAnalysis"
    mdtraj = "mdtraj"

class SafeDependencyTest:

    _tested_dependencies: dict[DependencyWithSafeTest, bool] = {}

    @classmethod
    def is_present(cls, dependency: DependencyWithSafeTest) -> bool:
        if dependency in cls._tested_dependencies.keys():
            return cls._tested_dependencies[dependency]

        result = False
        if dependency == DependencyWithSafeTest.MDAnalysis:
            result = True
            try:
                import MDAnalysis as mda
            except ModuleNotFoundError:
                result = False
        if dependency == DependencyWithSafeTest.mdtraj:
            result = True
            try:
                import mdtraj as md
            except ModuleNotFoundError:
                result = False

        cls._tested_dependencies[dependency] = result
        return result
