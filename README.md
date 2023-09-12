# GammaLoop

*Computation of differential cross-sections using Local Unitarity.*

See [www.alphaloop.ch](www.alphaloop.ch) for details on the theoretical framework and related literature.

See the [wiki](https://wiki.alphaloop.ch/) for more information on the project.

## Installation

### > Installation using `pip`
```
pip install gammaloop
```

### > Installation from sources
```
git clone https://github.com/alphal00p/gammaloop.git
cd gammaloop
python -m pip install -r ./python/gammaloop/requirements.txt
./bin/build_dependencies.sh
./bin/compile_bin.sh --release
./bin/compile_lib.sh --release
```
The relevant binaries will then be in `./bin/` and the python module in `./python/gammalop`.

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

Simply run:
```
/bin/run_tests.sh python
/bin/run_tests.sh rust
```

## Usage

`GammaLoop` can be used as a python command line interface (`gammaloop`) for scattering process code generation or also directly as a binary (`gammaloop_rust_cli`) for steering the Monte-Carlo integration. Both programs are installed as scripts in the `bin` directly when installing `GammaLoop` using `pip` or when installing it from sources.

### > Generating scattering process code with ./gammaloop

Place your list of commands in a file named e.g. `cmd.gL`, for instance:
```
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

### > Steering Monte-Carlo integration with ./gammaloop_rust_cli

Steering a Monte-Carlo integration with `gammaloop_rust_cli` can be done by running:
```
./bin/gammaloop_rust_cli --config <MODULE_PATH>/gammaloop/data/run_cards/rust_run_config.yaml
```
You will find more information on the content of the run configuration file and steering options in the [wiki](https://wiki.alphaloop.ch/) and by running:
```
./bin/gammaloop_rust_cli --help
```