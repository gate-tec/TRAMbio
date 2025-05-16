import sys
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
