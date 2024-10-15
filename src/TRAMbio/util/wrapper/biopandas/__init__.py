from TRAMbio.util.errors.dependencies import MissingDependencyError

try:
    import biopandas as bp
except ImportError:
    exc = MissingDependencyError(
        module="tram.util.wrapper.biopandas",
        dependency="biopandas"
    )
    try:
        raise exc
    finally:
        exc.__context__ = None
