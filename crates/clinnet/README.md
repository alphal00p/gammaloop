# clinnet

`clinnet` is the command-line companion to the [`linnet`](https://crates.io/crates/linnet) graph library.
It scans a directory for `.dot` files, incrementally renders each figure through
[Typst](https://typst.app), and assembles the resulting PDFs into a single grid document.

## Installation

```bash
cargo install clinnet
```

This installs the `linnet` binary (the package is named `clinnet` to avoid a clash with the library
crate). Make sure the `typst` CLI is available on your `$PATH`, as the tool shells out to it for
every compilation step.

## Quick start

```bash
# Render every .dot file underneath ./examples into build/grid.pdf
linnet examples
```

On first run, clinnet writes a small set of default Typst templates and a compiled `linnest.wasm`
plugin into `build/templates/`. You can provide your own templates via `--figure-template`,
`--grid-template`, or `--layout-template` flags. Pass `--style` multiple times to add extra files
whose contents should invalidate the incremental build cache when they change.

You can pass custom variables to Typst templates using `--input key=value`, which makes the
variables available through `sys.inputs` in your templates. This works exactly like the Typst CLI's
`--input` option.

## Typical workflow

1. Generate or edit `.dot` graphs in your project directory.
2. Run `linnet <graphs-root>` to refresh the figure PDFs and the combined grid.
3. Open `build/grid.pdf` (or the path supplied via `--grid-output`) to browse the results.

The cache lives in `build/.cache/figures.json`; deleting it forces a full rebuild. See
`linnet --help` for the complete flag list.

## Passing variables to templates

You can pass custom variables to your Typst templates:

```bash
# Pass variables that will be available in templates via sys.inputs
linnet examples --input theme=dark --input scale=1.5

# Multiple inputs can be specified
linnet examples --input author="John Doe" --input version=1.0 --input debug=true
```

These variables are available in both figure and grid templates through Typst's `sys.inputs`
dictionary.

## Partial rebuilds

After a full run, you can re-render individual figures or the final grid without rescanning DOT
files:

```bash
# Rebuild a single figure that was part of the previous run
linnet --build-dir build --rebuild-figure path/to/file.dot

# Rebuild the grid PDF based on the cached figure metadata
linnet --build-dir build --rebuild-grid

# Partial rebuilds also support input variables
linnet --build-dir build --input theme=light --rebuild-grid
```
