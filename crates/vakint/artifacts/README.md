# Vakint artifacts

Each loop directory contains `ibps/` for native RustRed reduction programs and
plans, `cache/` for explicitly loaded prepared reductions, and `masters/` for
independently computed master coefficients and scalar graphs.

One-, two- and three-loop `ibps/rules.rrbin` files contain their complete sealed
reduction programs. Four-loop candidate programs, normalization plans, ordered
family descriptors and existing FORM terminal catalogs are in `4L/ibps/`.
The historical codec receipt retains its original source paths; the relocation
manifest records the corresponding current paths and unchanged payload hashes.

Master values are archived here for inspection; Vakint does not load these JSON
files. Existing internal master tables and evaluator behavior remain unchanged.
Use `manifest.json` and `masters_summary.md` for stable identities and the
validated precision of each exported result.

Prepared caches remain opt-in: `engine.load_cache(n, path, &Default::default())`.
Pass the relevant `nL/cache` directory directly. No evaluation automatically reads
or modifies those packs. Five-loop cache preparation remains unavailable; its saved IBP application is
loaded only through the existing explicit application API.
See [the cache interface](../cache_preparation.typ).
