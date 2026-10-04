# Embedded Standard Model

`sm.json.zlib` is the complete `sm` model from the existing
`tests/fixtures/sm.json`. It preserves all parameters, evaluated values,
particles, propagators, Lorentz structures, couplings, and vertices.

The zlib stream is embedded in the library and decompressed on
demand. The loader needs no filesystem or Python/UFO runtime. A regression
test compares the entire decoded definition with the source fixture and
limits the compressed payload to 8 KiB.

Regenerate from the workspace root using Python's standard library:

```sh
python3 - <<'PY'
import json
from pathlib import Path
import zlib

root = Path("crates/feynkit-model")
model = json.loads((root / "tests/fixtures/sm.json").read_text())
compact = json.dumps(model, separators=(",", ":")).encode()
(root / "data/sm.json.zlib").write_bytes(zlib.compress(compact, level=9))
PY
```
