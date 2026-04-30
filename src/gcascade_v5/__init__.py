"""GCascadeV5: electron-tracking gamma-ray cascade transport."""

import importlib

from .core import *  # noqa: F401,F403
from . import config as _config

_state_module = importlib.import_module(".state", __name__)


__version__ = "5.1"

try:
    del EBLindex
except NameError:
    pass

try:
    del PROGRESS_ENABLED
except NameError:
    pass


def __getattr__(name: str):
    """Expose live runtime state values directly from the package namespace."""
    if name == "EBLindex":
        return _state_module.EBLindex
    if name == "PROGRESS_ENABLED":
        return _config.PROGRESS_ENABLED
    raise AttributeError(name)

