"""Read parameter labels before the UFO loader parses expressions."""

import importlib
import json
import sys
from pathlib import Path

from ufo_model_loader.common import DATA_PATH


def parameter_labels(input_model_path):
    # Match the loader's local-path and bundled-model lookup order.
    path = Path(input_model_path)
    if not path.exists() and not path.is_absolute():
        path = Path(DATA_PATH) / "models" / path
    path = path.absolute()
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".json":
        with path.open(encoding="utf-8") as source:
            definition = json.load(source)
            return [
                (
                    parameter["name"],
                    (parameter.get("texname"), parameter.get("typstname")),
                )
                for parameter in definition["parameters"]
            ]

    # Import the UFO declarations, without constructing symbolic expressions.
    # The upstream loader reuses this module through Python's import cache.
    previous_path = sys.path[:]
    try:
        sys.path[:0] = [str(path), str(path.parent)]
        model = importlib.import_module(path.name)
    finally:
        sys.path[:] = previous_path
    return [
        (
            parameter.name,
            (
                getattr(parameter, "texname", None),
                getattr(parameter, "typstname", None),
            ),
        )
        for parameter in [
            *model.all_parameters,
            *getattr(model, "all_CTparameters", ()),
        ]
    ]
