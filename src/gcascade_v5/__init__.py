from .core import *  # noqa: F401,F403
from . import core as _core


__version__ = "0.1.0"

try:
    del EBLindex
except NameError:
    pass


def __getattr__(name: str):
    if name == "EBLindex":
        return _core.EBLindex
    raise AttributeError(name)
