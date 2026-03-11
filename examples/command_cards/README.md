# GammaLoop Command File Examples

This directory contains example command files that demonstrate various GammaLoop workflows and capabilities.

## File Formats

GammaLoop command files use the TOML format with a specific schema. The basic structure is:

```toml
#:schema https://raw.githubusercontent.com/alphal00p/gammaloop/refs/heads/HEAD/assets/schemas/runhistory.json

commands = [
    "command1",
    "command2",
    # ... more commands
]

[cli_settings]
# CLI configuration options

[default_runtime_settings]
# Default integration and runtime parameters
```

## Example Files

### `simple_workflow.toml`
Demonstrates a basic GammaLoop workflow:
- Loading the Standard Model
- Generating simple QED processes (e⁺e⁻ → μ⁺μ⁻)
- Generating processes with QCD corrections
- Saving Feynman diagrams as DOT files

**Usage:**
```bash
./bin/gammaloop examples/command_cards/simple_workflow.toml
```

**Key Features:**
- Basic model import and process generation
- Display commands to monitor progress
- DOT file export for diagram visualization
- Simple configuration settings

### `advanced_integration.toml`
Showcases advanced integration and evaluation capabilities:
- Monte Carlo integration with multiple cores
- Phase-space point evaluation and inspection
- Performance benchmarking
- High-precision integration settings
- Advanced UV regularization and subtraction

**Usage:**
```bash
./bin/gammaloop examples/command_cards/advanced_integration.toml
```

**Key Features:**
- Multi-core integration setup
- Stability checks and validation
- Custom kinematics configuration
- Performance benchmarking
- Advanced compiler optimizations

## Running Examples

### Direct Execution
```bash
# Run a complete workflow from file (positional argument)
./bin/gammaloop examples/command_cards/simple_workflow.toml

# Use custom state folder with run card
./bin/gammaloop -s ./my_test_state examples/command_cards/simple_workflow.toml

# Dry run (parse but don't execute)
./bin/gammaloop -d examples/command_cards/advanced_integration.toml
```

### Interactive Mode
```bash
# Start interactive session
./bin/gammaloop

# Then manually execute commands from the examples:
> import model sm-default
> generate amp e+ e- > mu+ mu- QED==2
> display processes
> save dot ee_to_mumu
```

### Command-Line Execution
```bash
# Execute single commands directly
./bin/gammaloop run -c "import model sm-default"
./bin/gammaloop run -c "generate amp e+ e- > mu+ mu- QED==2; display processes"
```

## Creating Custom Command Files

### Basic Template
```toml
#:schema https://raw.githubusercontent.com/alphal00p/gammaloop/refs/heads/HEAD/assets/schemas/runhistory.json

commands = [
    "import model sm-default",
    "generate amp <process_definition>",
    "display processes"
]

[cli_settings]
state_folder = "./my_custom_state"

[cli_settings.global]
logfile_directive = "info"
display_directive = "info"
```

### Common Command Patterns

#### Model and Process Management
```toml
commands = [
    "import model sm-default",                    # Load Standard Model
    "display model",                             # Show model details
    "generate amp e+ e- > mu+ mu- QED==2",  # Generate QED process
    "generate xs e+ e- > d d~ QCD=1",            # Generate QCD cross-section
    "display processes",                         # List all processes
    "save dot my_diagrams"                       # Export Feynman diagrams
]
```

#### Integration Workflow
```toml
commands = [
    "import model sm-default",
    "generate amp e+ e- > mu+ mu- QED==2",
    "set kinematics.e_cm 91.1876",              # Set center-of-mass energy
    "set integrator.n_max 1000000",             # Set max evaluations
    "integrate -p 0 --n-cores 4",               # Run integration
    "inspect -p 0 -x 0.1 0.2 0.3"               # Inspect specific point
]
```

## Configuration Options

### CLI Settings
- `state_folder`: Directory for persistent state storage
- `override_state`: Automatically overwrite existing files
- `try_strings`: Use string serialization when possible

### Runtime Settings
- `general.m_uv`: UV mass parameter
- `general.mu_r_sq`: Renormalization scale squared
- `kinematics.e_cm`: Center-of-mass energy
- `integrator.n_max`: Maximum number of evaluations
- `integrator.max_prob_ratio`: Maximum probability ratio for sampling

### Global Settings
- `logfile_directive`: Log level for file output
- `display_directive`: Log level for console output
- `generation.enable_thresholds`: Enable threshold computations
- `generation.compile.optimization_level`: Compiler optimization level

## Output and Results

### State Files
Running these examples will create state directories containing:
- `model.json`: Physics model definition
- `processes/`: Generated Feynman diagrams and integrands
- `run.toml`: Command history for reproducibility
- `*.log`: Session log files

### Visualization
DOT files can be converted to PDF diagrams:
```bash
cd <state_folder>/processes/<process_name>/drawings/dot
make -j
ls *.pdf  # View generated diagram PDFs
```

### Integration Results
Integration results are saved in specified output files and can be analyzed further with external tools.

## Troubleshooting

### Common Issues
1. **Model not found**: Ensure the model file exists or use `sm-default` for the built-in Standard Model
2. **State folder conflicts**: Use `-o` flag to override existing state or specify a new folder with `-s`
3. **Integration failures**: Check log files for detailed error messages
4. **Memory issues**: Reduce `n_max` or use fewer cores for large processes

### Debug Mode
```bash
# Enable verbose logging
GL_DISPLAY_FILTER=debug GL_LOGFILE_FILTER=debug ./bin/gammaloop -t debug.log examples/command_cards/simple_workflow.toml

# Check the generated log file
cat debug.log
```

For more detailed information, consult the [GammaLoop wiki](https://wiki.alphaloop.ch/) and the main project documentation.
