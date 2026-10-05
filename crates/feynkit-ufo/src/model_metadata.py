"""Read parameter labels and retain the requested UFO source metadata."""

import json
from fractions import Fraction
from pathlib import Path

from ufo_model_loader.common import DATA_PATH


def quantum_token(particle, field):
    value = particle.get(field)
    if isinstance(value, str):
        # The upstream JSON-file importer accepts Python's exact fraction strings.
        try:
            value = str(Fraction(value))
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(
                f"invalid model JSON field '{field}' for particle '{particle['name']}': {exc}"
            ) from exc
    return json.dumps(value)


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
                {
                    particle["name"]: tuple(
                        quantum_token(particle, field)
                        for field in ("charge", "y_charge", "y_charge_right")
                    )
                    for particle in definition["particles"]
                },
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
            None,
        )
