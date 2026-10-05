"""Read parameter labels and retain the requested UFO source metadata."""

import json
from pathlib import Path

from ufo_model_loader.common import DATA_PATH


def model_metadata(input_model_path):
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
            return (
                [
                    (
                        parameter["name"],
                        (parameter.get("texname"), parameter.get("typstname")),
                    )
                    for parameter in definition["parameters"]
                ],
                None,
            )

    from ufo_model_loader.raw import imported_ufo_module

    # The loader isolates same-name UFO imports. Retain this source module so
    # metadata never comes from a different module restored in sys.modules.
    # Context cleanup runs before the Rust caller registers symbolic printers.
    with imported_ufo_module(path) as model:
        return (
            [
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
            ],
            model,
        )
