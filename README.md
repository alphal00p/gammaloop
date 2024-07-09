# gammaLoop
![tests status](https://github.com/alphal00p/gammaloop/actions/workflows/gamma_loop_tests.yml/badge.svg?event=push)

*Computation of differential cross-sections using Local Unitarity.*

See [www.alphaloop.ch](www.alphaloop.ch) for details on the theoretical framework and related literature.

See the [wiki](https://wiki.alphaloop.ch/) for more information on the project.

## TL;DR: I want to be running γLoop already!

We salute the eagerness of our users.

If you want to jump right in, run the following to immediately start integrating the scalar three-loop mercedes diagram!

```
git clone https://github.com/alphal00p/gammaloop.git && cd gammaloop
./bin/compile.sh
./bin/gammaloop example/cards/scalar_mercedes.gL
```

## Installation

### > Requirements

* `Rust`: v1.77+ You can easily install Rust with this [one-liner](https://www.rust-lang.org/tools/install)

* `Python3`: v3.11+ (equipped with `pip`, and the `python-devel`)

* `GNU gcc`: v10+ (*not* `clang`!)

* `git`

*Note*: The project has been tested on Linux and MacOS. However, on MacOS, the default `clang` compiler is not supported due to lack of `libquadmath` support, and `GNU gcc` must be installed and setup as default compiler.

Windows users are encouraged to use [WSL 2](https://learn.microsoft.com/en-us/windows/wsl/).

### > Installation using `pip`
```
pip install gammaloop
gammaloop --build_dependencies
```

### > Installation from sources
```
git clone https://github.com/alphal00p/gammaloop.git
cd gammaloop
./bin/compile.sh
```
The relevant binaries will then be in `./bin/` and the gammaloop python module is located at `./python/gammaloop`.

*Note:* Alternatively, the dependencies can be built within a python virtual environment as follows:

```
./bin/build_dependencies.sh clean
./bin/build_dependencies.sh with_venv
source `./bin/gammaloop -venv`
```

and when using the `fish` shell, the last command should be replaced by: `. (./bin/gammaloop -venv).fish`

## Tests

### > Testing an installation from `pip`

Testing your installation can be done by running
```
python -m pytest
```
from within the installation directory of the gammaloop module, which you can `cd` into with e.g.:
```
bash -c 'cd `python -c "import os; import gammaloop; print(os.path.dirname(gammaloop.__file__))"`; pwd'
```

### > Testing an installation from sources

Run:
```
/bin/run_tests.sh python
/bin/run_tests.sh rust
```

## Usage

There are three entry points to the `GammaLoop` functionalities:

1. Preferred method is through the Python command-line interface `gammaloop`.

2. Alternatively, the same functionalities can be accessed programmatically, e.g. in a Jupyter notebook, through the Python API, by importing the `gammaloop` library.

3. Finally, expert users may also find it useful to steer some of functionalities directly from the rust binary `gammaloop_rust_cli`.

Both executables `gammaloop` and `gammaloop_rust_cli` are made available as scripts in the `bin` directory.
The `gammaloop` Python module is also exposed after installation and ready to be imported in user custom Python scripts.

### 1. Usage from the Python command-line interface: ./gammaloop

`GammaLoop` is typically used through the python command-line interface `gammaloop`.
Place your list of commands in a file named e.g. `cmd.gL`, for instance:

```
# Content of file 'cmd.gL'
import_model sm --format ufo
export_model ./sm.yaml --format yaml
```
and run it with:
```
./gammaloop cmd.gL
```
You can find more information on the syntax of the available commands in the [wiki](https://wiki.alphaloop.ch/) and by running:
```
./gammaloop --help
```
to get an overview of available commands and:
```
./bin/gammaloop -c "help import_model"
```
to get help on any specific command.

You can find example of command files in the `<MODULE_PATH>/data/run_cards/` directory.

### 2. Usage from within a Jupyter notebook: the Python API

Follow the example jupyter notebook given in example to get started with the Python API.
```
cd <GAMMALOOP_INSTALLATION_DIRECTORY>/examples/jupyter
jupyter notebook steer_gammaloop.ipynb
``` 

### 3. Usage from the rust binary executable: ./gammaloop_rust_cli

All typical usecases of `GammaLoop` are available through the Python command-line interface mentioned earlier.
However, expert users might want to steer the Monte-Carlo integration directly using the `gammaloop_rust_cli` binary.
This is possible throught the `gammaloop_rust_cli` binary, which can be used as follows:
```
./bin/gammaloop_rust_cli --config <MODULE_PATH>/gammaloop/data/run_cards/rust_run_config.yaml
```
*Note*: You have to manually define your desired external momenta in this default `rust_run_config.yaml` file.

You will find more information on the content of the run configuration file and steering options in the [wiki](https://wiki.alphaloop.ch/) and by running:
```
./bin/gammaloop_rust_cli --help
```
