import importlib

from .core import *  # noqa: F401,F403
from . import config as _config
from . import legacy  # noqa: F401

_state_module = importlib.import_module(".state", __name__)


__version__ = "5.0"

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
        return _state_module.EBLindex
    if name == "NUMBA_ENABLED":
        return _config.NUMBA_ENABLED
    if name == "PROGRESS_ENABLED":
        return _config.PROGRESS_ENABLED
    raise AttributeError(name)
