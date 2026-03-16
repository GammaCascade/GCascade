from .core import *  # noqa: F401,F403
from . import core as _core


__version__ = "0.1.0"

try:
    del EBLindex
except NameError:
    pass

try:
    del NUMBA_ENABLED
except NameError:
    pass

try:
    del PROGRESS_ENABLED
except NameError:
    pass


def __getattr__(name: str):
    if name == "EBLindex":
        return _core.EBLindex
    if name == "NUMBA_ENABLED":
        return _core.NUMBA_ENABLED
    if name == "PROGRESS_ENABLED":
        return _core.PROGRESS_ENABLED
    raise AttributeError(name)
