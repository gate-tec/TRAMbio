import warnings
from typing import Type, List, Union

try:
    from pandas.errors import SettingWithCopyWarning
    HAS_SETTING_WITH_COPY_WARNING = True
except ImportError:
    class SettingWithCopyWarning(Warning):
        pass
    HAS_SETTING_WITH_COPY_WARNING = False


class WarningWrapper:

    def __init__(self, warning_type: Union[List[Type[Warning]], Type[Warning]] = None):
        self._warning_list = [Warning] if warning_type is None else (
            warning_type if isinstance(warning_type, list) else [warning_type]
        )
        self.catch_warning_object = warnings.catch_warnings()

    def __enter__(self):
        self.catch_warning_object.__enter__()
        if HAS_SETTING_WITH_COPY_WARNING:
            warnings.simplefilter(action='ignore', category=SettingWithCopyWarning)
        for warning_type in self._warning_list:
            warnings.simplefilter(action='ignore', category=warning_type)

    def __exit__(self, *exc_info):
        self.catch_warning_object.__exit__(*exc_info)
