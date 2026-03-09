<div align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://github.com/alphal00p/gammaloop/blob/2ee2ec575fa575c26bdaf89a3e7df41428b879dc/assets/gammalooplogo-dark.svg">
  <img src="https://github.com/alphal00p/gammaloop/blob/2ee2ec575fa575c26bdaf89a3e7df41428b879dc/assets/gammalooplogo-light.svg", width="300">
</picture>
</div> 


*Computation of differential cross-sections using Local Unitarity.*

See [www.alphaloop.ch](https://www.alphaloop.ch) for details on the theoretical framework and related literature.

See the [wiki](https://wiki.alphaloop.ch/) for more information on the project.
Architecture documentation is available under [docs/architecture](docs/architecture/architecture.md).

## TL;DR: I want to be running γLoop already!

We salute the eagerness of our users!

If you want to jump right in, run the following to immediately start integrating the scalar three-loop mercedes diagram,

```bash
git clone https://github.com/alphal00p/gammaloop.git && cd gammaloop
just build-cli
./bin/gammaloop examples/command_cards/simple_workflow.toml
```

or a one-loop hexagon diagram from the scattering process $\gamma \gamma \rightarrow \gamma \gamma \gamma \gamma$ (through a top-quark loop) with:

```bash
./bin/gammaloop examples/command_cards/advanced_integration.toml
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

* `just`

* `maturin` (required for `just build-api`)

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
just build-cli
just build-api
```
`just build-cli` builds the CLI, and `just build-api` builds/installs the Python API via `maturin develop` (in your active Python environment).
The Python package source is located in `./crates/gammaloop-api/python/gammaloop`.

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
Place your list of commands in a file named e.g. `cmd.toml`, for instance:

```bash
# Content of file 'cmd.toml'
commands = [
    "import model sm --format ufo",
    "export model ./sm.yaml --format yaml"
]
```
and run it with:
```bash
./bin/gammaloop cmd.toml
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

You can find example command files in the `examples/command_cards/` and `tests/resources/run_cards/` directories.

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

## CLI Usage and State File Management

### Command-Line Interface Overview

The `gammaloop` CLI provides a comprehensive interface for managing physics computations with persistent state management. The CLI automatically maintains state between sessions through a dedicated state folder.

### State Folder Structure

By default, GammaLoop creates a `./gammaloop_state/` directory containing:

```
gammaloop_state/
├── cli_settings.toml              # CLI configuration and global settings
├── default_runtime_settings.toml  # Default integration and runtime parameters
├── model.json                     # Current physics model definition
├── model_parameters.json          # Model parameter values
├── run.toml                       # Command history and execution log
├── symbolica_state.bin            # Binary state from computer algebra system
├── processes/                     # Generated process definitions and integrands
└── logs/                         # Session log files
```

### CLI Command Structure

The CLI supports the following main command categories:

#### Model and Process Management
- **`import model <model_name>`** - Load a physics model (e.g., Standard Model)
- **`import graphs <file.dot>`** - Import Feynman diagrams from DOT files
- **`generate <process>`** - Generate amplitudes or cross-sections for physics processes
- **`display processes`** - Show currently loaded processes
- **`display model`** - Show current model information

#### State Management
- **`save state <path>`** - Save current state to a specified location
- **`save dot <name>`** - Export Feynman diagrams as DOT files
- **`reset all`** - Clear all processes and reset state
- **`set <parameter> <value>`** - Modify configuration parameters

#### Integration and Evaluation
- **`integrate <options>`** - Run Monte Carlo integration
- **`inspect <options>`** - Evaluate at specific phase-space points
- **`evaluate <options>`** - Batch evaluation of integrands
- **`bench <options>`** - Benchmark integrand evaluation performance

#### Process References
Commands that target a specific process accept `--process <ref>` (short `-p`), where `<ref>` can be:
- `#<id>` - Explicit numeric id
- `name:<process_name>` - Explicit process name
- `<id>` or `<process_name>` - Implicit reference (if ambiguous, use a prefix)

Examples:
```bash
./bin/gammaloop -c "display processes"
./bin/gammaloop -c "display integrands --process #0"
./bin/gammaloop -c "integrate --process name:pp_ttbar --integrand amp_0"
```

### CLI Options and State Control

#### State Folder Management
```bash
# Use custom state folder
./bin/gammaloop -s /path/to/custom_state

# Recommended for local experimentation (keeps repo root clean)
./bin/gammaloop -s .local/scratch/my_run/gammaloop_state

# Override existing state files without prompting
./bin/gammaloop -o

# Don't save state on exit
./bin/gammaloop -n

# Load specific model file
./bin/gammaloop -m /path/to/model.json
```

#### Running Commands
```bash
# Interactive mode
./bin/gammaloop

# Execute single command
./bin/gammaloop -c "display processes"

# Run commands from file
./bin/gammaloop commands.toml

# Execute run card file (positional argument)
./bin/gammaloop saved_session.toml

# Dry run (parse but don't execute)
./bin/gammaloop -d commands.toml
```

#### Logging and Debugging
```bash
# Set log level
./bin/gammaloop -l debug

# Custom log file
./bin/gammaloop -t my_session.log

# Environment variables for fine-grained logging
GL_DISPLAY_FILTER=info GL_LOGFILE_FILTER=debug ./bin/gammaloop
```

### State Persistence and Session Management

#### Automatic State Saving

GammaLoop automatically saves state after each command execution, preserving:
- Loaded physics models and parameters
- Generated processes and integrands
- Configuration settings
- Command history

#### Command History and Replay
The `run.toml` file maintains a complete history of executed commands:

```gammaloop/gammaloop_state/run.toml#L1-10
#:schema https://raw.githubusercontent.com/alphal00p/gammaloop/refs/heads/HEAD/assets/schemas/runhistory.json
commands = [
    "import model sm-default",
    "generate xs e+ e- > d d~ / z QED^2==4 [{{1}} QCD]",
    "save dot my_diagrams",
]

[cli_settings]
state_folder = "./gammaloop_state"
```

Sessions can be replayed or continued:
```bash
# Continue previous session
./bin/gammaloop -s existing_state_folder

# Replay specific run card
./bin/gammaloop previous_session.toml
```

#### Configuration Management
Settings are hierarchically managed through:
1. **Global defaults** - Built-in default values
2. **State folder settings** - Persistent per-project configuration
3. **Command-line options** - Session-specific overrides
4. **Environment variables** - Runtime logging and display control

#### State Migration and Sharing
```bash
# Export complete state for sharing
./bin/gammaloop -c "save state /path/to/shared_project"

# Load shared state
./bin/gammaloop -s /path/to/shared_project

# Load and execute run card
./bin/gammaloop base.toml
```

### Example Workflow

```bash
# Start with a simple example
./bin/gammaloop examples/command_cards/simple_workflow.toml

# Run advanced integration example
./bin/gammaloop examples/command_cards/advanced_integration.toml

# Start new project with custom state folder
./bin/gammaloop -s my_project

# Start an isolated local scratch run
./bin/gammaloop -s .local/scratch/my_project/gammaloop_state

# Create and run custom command file
cat > my_commands.toml << EOF
commands = [
    "import model sm",
    "generate amplitude e+ e- > mu+ mu- QED=2",
    "save dot electron_muon_scattering",
    "integrate --n-cores 4 --target-accuracy 1e-6"
]
EOF
./bin/gammaloop -s my_project my_commands.toml

# Continue work later
./bin/gammaloop -s my_project -c "display processes"

# Share project
cp -r my_project shared_project
./bin/gammaloop -s shared_project
```

This state management system ensures reproducibility, enables collaboration, and provides robust session persistence for complex physics computations.
