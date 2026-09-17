import hashlib
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
p = root / "workspace/manifest.json"
old = json.loads(p.read_text())
new = json.loads((root / "fresh_same_process/manifest.json").read_text())
assert old["slots"] == new["slots"]
assert old["effective_model_parameters"] == new["effective_model_parameters"]
backup = root / "original_manifest.json"
if not backup.exists():
    backup.write_text(p.read_text())
change = {
    "old": old["integrand_fingerprints"],
    "current_process": new["integrand_fingerprints"],
    "checkpoint_sha256": hashlib.sha256(
        (root / "workspace/state/integration_state.bin").read_bytes()
    ).hexdigest(),
}
(root / "diagnostic_fingerprint_override.json").write_text(json.dumps(change, indent=2))
old["integrand_fingerprints"] = new["integrand_fingerprints"]
p.write_text(json.dumps(old, indent=2))
print(
    "Diagnostic clone only: matched current process fingerprint; original checkpoint preserved."
)
