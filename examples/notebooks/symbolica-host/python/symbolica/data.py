"""Locate installed resource wheels for synchronous native operations."""

from __future__ import annotations

import json
import os
from pathlib import Path

from . import core


def manifest() -> dict:
    """Describe this host's expected resources, including loaded status."""
    return json.loads(core._data_manifest())


def set_directory(directory: str | os.PathLike[str]) -> None:
    """Use one extracted resource bundle without eagerly reading its files."""
    set_directories(directory)


def set_directories(*directories: str | os.PathLike[str]) -> None:
    """Search extracted resource bundles in the supplied order."""
    core._set_data_directories([os.fspath(directory) for directory in directories])


def configure() -> None:
    """Prefer explicit environment paths, otherwise find installed data wheels."""
    if directory := os.environ.get("SYMBOLICA_DATA_DIRECTORY"):
        set_directories(*directory.split(os.pathsep))
    else:
        package_data = Path(__file__).parent / "_data"
        set_directories(*(package_data / bundle for bundle in ("fonts", "vakint-data")))
