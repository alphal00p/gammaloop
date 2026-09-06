# Kurvst Coil Benchmarks

These Typst documents compare compile time for CeTZ's native
`cetz.decorations.coil` against Kurvst's `coil` path pattern rendered back
through CeTZ.

Run from the package root:

```sh
bash benches/run.sh
```

The runner accepts optional environment overrides:

```sh
COUNT=400 REPEATS=8 bash benches/run.sh
```

- `COUNT` controls how many coils each document draws.
- `REPEATS` controls how many benchmark runs are collected.

The runner uses `hyperfine` if it is installed. Otherwise it falls back to a
small shell `time` loop.
