# Embedded Standard Model

The complete `sm` model in `sm.json.zlib` comes from
`assets/models/json/sm/sm.json` in the workspace. It preserves all parameters, evaluated values,
particles, propagators, Lorentz structures, couplings, and vertices.

The zlib stream is embedded in the library and decompressed on
demand. The loader needs no filesystem or Python/UFO runtime. A regression
test compares the entire decoded definition with the source fixture and
limits the compressed payload to 8 KiB. The test fixture links to the workspace
asset. Cargo includes the file contents when packaging the crate.

Regenerate from the workspace root using Python's standard library:

```sh
python3 - <<'PY'
import json
from pathlib import Path
import zlib

root = Path("crates/feynkit-model")
model = json.loads(Path("assets/models/json/sm/sm.json").read_text())
compact = json.dumps(model, separators=(",", ":")).encode()
(root / "data/sm.json.zlib").write_bytes(zlib.compress(compact, level=9))
PY
```
