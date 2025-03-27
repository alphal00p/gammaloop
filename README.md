<div align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/alphal00p/gammaloop/main/python/gammaloop/data/assets/gammalooplogo-dark.svg">
  <img src="https://raw.githubusercontent.com/alphal00p/gammaloop/main/python/gammaloop/data/assets/gammalooplogo-light.svg", width="300">
</picture>
</div> 

![tests status](https://github.com/alphal00p/gammaloop/actions/workflows/gamma_loop_tests.yml/badge.svg?event=push)

*Computation of differential cross-sections using Local Unitarity.*

See [www.alphaloop.ch](https://www.alphaloop.ch) for details on the theoretical framework and related literature.

See the [wiki](https://wiki.alphaloop.ch/) for more information on the project.

## TL;DR: I want to be running γLoop already!

We salute the eagerness of our users!

If you want to jump right in, run the following to immediately start integrating the scalar three-loop mercedes diagram,

```bash
git clone https://github.com/alphal00p/gammaloop.git && cd gammaloop
./bin/compile.sh
./bin/gammaloop example/cards/scalar_mercedes.gL
```

or a one-loop hexagon diagram from the scattering process $\gamma \gamma \rightarrow \gamma \gamma \gamma \gamma$ (through a top-quark loop) with:

```bash
./bin/gammaloop examples/cards/physical_1L_AA_AAAA.gL
```

You can also generate diagrams for the computation of cross-sections or amplitudes with:
```bash
> ./bin/gammaloop -c "import_model sm;\
generate_amplitude d d~ > d d~ e+ e- | u d g ghg e- a QED=2 [{2} QCD=2];\
output MyProcess -exp -mr -ef mathematica --yaml_only;"
[...]
INFO    : A total of 2214 graphs have been generated in 00:00:49.142.
```
You can then visualize all diagrams with: 
```bash
> cd ./MyProcess/sources/amplitudes/GL_ddx_ddxepem/drawings/dot
> make -j &> /dev/null
> ls feynman_diagrams.pdf
```
and import their corresponding symbolic expressions from:
```bash
> ls ./MyProcess/sources/amplitudes/GL_ddx_ddxepem/expressions
[...] GL1016_exp.json  GL1223_exp.json  GL1426_exp.json  GL1816_exp.json [...]
> head -2 ./MyProcess/sources/amplitudes/GL_ddx_ddxepem/expressions/GL84_exp.json
[
  "-2/3 𝑖 ee^2 G^6 Metric[mink[4,6],mink[4,13]] Metric[mink[4,7],mink[4,12]] Metric[mink[4,8],mink[4,11]] Metric[mink[4,9],[mink[4,10]] T[coad[8,10],cof[3,9],dind[cof[3,8]]] T[coad[8,13],cof[3,12],dind[cof[3,11]]] [...]",
```
which can be further processed to your liking using e.g. [Mathematica](https://www.wolfram.com/mathematica), [Symbolica](https://github.com/benruijl/symbolica-community) and / or [spenso](https://github.com/alphal00p/spenso).

Process generation syntax is detailed on γLoop's [wiki](https://wiki.alphaloop.ch/en/gammaLoop/ProcessGeneration).

## Installation

### > Requirements

The installation may be successful with older versions than the ones indicated below, but it was not tested. 

* `Rust`: v1.81+ You can easily install Rust with this [one-liner](https://www.rust-lang.org/tools/install)

* `Python3`: v3.12+ (equipped with `pip`, and the `python-devel` dependency)

* `GNU gcc`: v14+

* `git`

Windows users are encouraged to use [WSL 2](https://learn.microsoft.com/en-us/windows/wsl/).

### > Installation using `pip`
```bash
pip install gammaloop
gammaloop --build_dependencies
```

### > Installation from sources
```bash
git clone https://github.com/alphal00p/gammaloop.git
cd gammaloop
./bin/compile.sh --release
```
The relevant binaries will then be in `./bin/` and the gammaloop python module is located at `./python/gammaloop`.

*Note:* Alternatively, the dependencies can be built within a python virtual environment as follows:

```bash
./bin/build_dependencies.sh clean
./bin/build_dependencies.sh with_venv
# From within a `bash` shell, use the command below
source `./bin/gammaloop -venv`
# From within a `fish` shell, use the command below
. (./bin/gammaloop -venv).fish
```

## Tests

### > Testing an installation from `pip`

From within the installation directory of the gammaloop module, which you can `cd` into with e.g.:
```bash
bash -c 'cd `python -c "import os; import gammaloop; print(os.path.dirname(gammaloop.__file__))"`; pwd'
```

You can test your `gammaLoop` installation by running (incl. only tests with a maximum runtime of 15 seconds):
```bash
python -m pytest --max-runtime 15.0
```

### > Testing an installation from sources

From within the installation directory, run:
```bash
/bin/run_tests.sh python
/bin/run_tests.sh rust
```

## Usage

There are three entry points to the `GammaLoop` functionalities:

1. Preferred method is through the Python command-line interface `gammaloop`.

2. Alternatively, the same functionalities can be accessed programmatically, e.g. in a Jupyter notebook, through the Python API, by importing the `gammaloop` library.

3. Finally, expert users may also find it useful to steer some of functionalities directly from the rust binary `gammaloop_rust_cli`.

Both executables `gammaloop` and `gammaloop_rust_cli` are made available as scripts in the `bin` directory.
The `gammaloop` Python module is also exposed after installation and ready to be imported by custom Python scripts.

### 1. Usage from the Python command-line interface: ./gammaloop

`GammaLoop` is typically used through the python command-line interface `gammaloop`.
Place your list of commands in a file named e.g. `cmd.gL`, for instance:

```bash
# Content of file 'cmd.gL'
import_model sm --format ufo
export_model ./sm.yaml --format yaml
```
and run it with:
```bash
./bin/gammaloop cmd.gL
```
You can find more information on the syntax of the available commands in the [wiki](https://wiki.alphaloop.ch/) and by running:
```bash
./bin/gammaloop --help
```
to get an overview of available commands and:
```bash
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
```bash
./bin/gammaloop_rust_cli --config <MODULE_PATH>/gammaloop/data/run_cards/rust_run_config.yaml
```
*Note*: You have to manually define your desired external momenta in this default `rust_run_config.yaml` file.

You will find more information on the content of the run configuration file and steering options in the [wiki](https://wiki.alphaloop.ch/) and by running:
```bash
./bin/gammaloop_rust_cli --help
```
