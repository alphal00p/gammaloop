#import "../../shared.typ": callout

#let quickstart-python = [
= Using GammaLoop from Python

This route generates the same scalar bubble as the command-line example, but keeps the session in a
Python object and returns a structured evaluation result. Use Python 3.11 or newer.

== Install the Python package

Create an isolated environment. The ordinary wheel command targets the automated `0.3.4`
release and newer:

// docs-example: syntax
```sh
python3.11 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install "gammaloop>=0.3.4"
```

Until that release is visible on PyPI, build the current package from a source checkout instead.
The public build selector is not a private license key:

// docs-example: syntax
```sh
git clone https://github.com/alphal00p/gammaloop.git
cd gammaloop
python3.11 -m venv .venv
. .venv/bin/activate
SYMBOLICA_OEM_LICENSE=SYMBOLICA_OEM_GAMMALOOP \
  python -m pip install .
```

== Generate and evaluate

Save this as `first_gammaloop.py`:

// docs-example: compile gammaloop-python-quickstart
```python
from pathlib import Path
from tempfile import TemporaryDirectory

from gammaloop import GammaLoopAPI

with TemporaryDirectory() as temporary:
    api = GammaLoopAPI(state_folder=Path(temporary) / "state")
    api.run("import model scalars-default.json")
    api.run(
        "set global kv "
        "global.generation.uv.generate_integrated=false"
    )
    api.run(
        "generate amp scalar_1 > scalar_1 [{1}] "
        "--allowed-vertex-interactions V_3_SCALAR_122 "
        "-p bubble -i one_loop"
    )

    result = api.evaluate_sample(
        [0.1, 0.2, 0.3],
        process_id=0,
        integrand_name="one_loop",
    )
    print(result.integrand_result)
```

Run `python first_gammaloop.py`. A successful execution prints a finite complex value. The
temporary state is deleted when the script exits; use a normal directory when the state should
be persisted and explicitly execute a save command before closing the session.

#callout("One session owns one mutable state", [
  `GammaLoopAPI.run` acts on the same in-memory model, settings, processes, and history as the
  structured methods. The object does not save automatically. Use separate instances for
  independent calculations, and inspect the #link("reference/interfaces/")[interface guide]
  before building a long-running workflow around mutable state.
])

The documentation harness syntax-checks this script. Runtime execution additionally needs the
native wheel or source build and a Symbolica license mode appropriate for the user.
]
