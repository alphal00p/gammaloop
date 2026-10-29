set -euo pipefail
source docs/architecture/evidence/generation-integration/mark-geometry/round1-correction/environment.sh
python3 target/mark-geometry-round3/check-docs.py > target/mark-geometry-round3/registry.log 2>&1
