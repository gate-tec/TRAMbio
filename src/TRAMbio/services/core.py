from typing import TypeVar, Generic, Dict, Union, List
import abc

from TRAMbio.util.wrapper.base.list import first


class IBaseService(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        if hasattr(subclass, 'name'):
            return True
        return NotImplemented

    @property
    @abc.abstractmethod
    def name(self) -> str:
        raise NotImplementedError


_T = TypeVar('_T', bound=IBaseService)
_E = TypeVar('_E', bound=Exception)


class BaseServiceRegistry(Generic[_T]):

    def __init__(self):
        self.__services: Dict[str, _T] = {}

    def register_service(self, service: _T):
        if service is not None:
            self.__services[service.name] = service

    def list_service_names(self) -> List[str]:
        return list(self.__services.keys())

    def list_services(self) -> List[_T]:
        return list(self.__services.values())

    def query_service(self, service_name: str, error_on_none: bool = True) -> Union[_T, None]:
        if error_on_none and service_name not in self.list_service_names():
            raise KeyError(f"Service {service_name} is not registered")
        return self.__services.get(service_name, None)

    def single_service(self, error_on_none: bool = True):
        if len(self.list_service_names()) == 0:
            raise ValueError("No service registered")
        return self.query_service(first(self.list_service_names()), error_on_none=error_on_none)


class ExtendedServiceRegistry(Generic[_T, _E], BaseServiceRegistry[_T]):

    def __init__(self):
        super().__init__()
        self.__exceptions: Dict[str, _E] = {}

    def list_service_names(self) -> List[str]:
        return super().list_service_names() + list(self.__exceptions.keys())

    def register_exception(self, name: str, exception: _E):
        self.__exceptions[name] = exception

    def query_service(self, service_name: str, error_on_none: bool = True) -> Union[_T, None]:
        service = super().query_service(service_name=service_name, error_on_none=error_on_none)
        if service is not None:
            return service
        if service_name in self.__exceptions.keys():
            raise self.__exceptions.get(service_name)
        return None
