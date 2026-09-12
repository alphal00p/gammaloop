# Replay the GL638 numerical counterexamples

[replay-points.json](replay-points.json) contains five recorded raw points and
13 comparisons: the ordinary double-soft precision failure, exact/nearby A/P
roots, and two points approaching a hard local-CT sliver boundary. Coordinates
were copied from the actual full-state manifests, without reconstructing them
from rounded geometry. Source hashes and historical results are included; the
private research archive and PDFs are not required to supply the inputs.

**Historical validation:** the exact Python block below executed all 13 comparisons
against the existing compatible full-sum state, using the `530930…` binary.
All totals, final precision/status/error records, and per-cut weights match the
historical values exactly. Eleven evaluations are finite without a NaN flag
(seven Stable, four Unstable); both exact A/P cases still flag `is_nan=true`.
Their sanitized zeros are reproduced failures, not physical successes. The run
finished in 64.56 s with 5.37 GB peak process-tree memory under a 30 GB guard.
The JSON records full binary and executed-script hashes.

After the nonfinite-status correction, release `57ccd168…` completed a further
19-case replay in 93.27 s: these 13 controls, point 5586 with normal/Arb Euler
checks, and both integration maxima with normal/Arb checks. The source state was
unchanged. All 17 finite cases retain bit-identical numerical payloads and
precision flags. Both exact A/P cases now finish Arb **Unstable** with no
accuracy estimate; their null cut-1 weights, `is_nan=true` and sanitized zero
totals are unchanged. Normal evaluation retains its finite intermediate Quad
rotation error of 0.0624291. The coincident-root failure and nearby-A/P false
acceptance remain. The JSON's historical observations are preserved, including
the old incorrect Stable labels; they are not the expected current labels for
the two nonfinite cases.

The first attempt failed before evaluation because loading the generation card
conflicted with saved command blocks in the read-only state; that attempt is
preserved separately. The corrected replay uses a runtime-only card. It did not
rerun the generation command below. Historical observations are diagnostic
references, not strict cross-platform assertions: regeneration can change
cancellation roundoff and the precision level accepted by the rotation test.

Run from the repository root with a built CLI and the normal project environment.
Generate a separate state using the committed GL638 card and its DOT metadata:

```bash
export GL638_BIN="$PWD/target/release/gammaloop"
export GL638_CARD="$PWD/examples/cli/epem_a_ttxh/NNLO/epem_a_tth_NNLO_test_GL638.toml"
export GL638_STATE="$(mktemp -d "${TMPDIR:-/tmp}/gl638-replay-state.XXXXXX")"
export GL638_REPLAY_DIR="$(mktemp -d "${TMPDIR:-/tmp}/gl638-replay-output.XXXXXX")"
"$GL638_BIN" --state-folder "$GL638_STATE" --no-save-state "$GL638_CARD" run generate
```

Generation imports `sm-default.json`, applies the model parameters recorded in
the JSON, generates all 936 source orientations as an explicit sum, and saves the
state. It includes direct 3D local UV and integrated UV. Generation is expensive;
an existing compatible full-sum state can instead be assigned to `GL638_STATE`.
Do not use a selected-orientation state. Regenerate after changing DOT metadata
or generated exact-π coefficients; local-CT profiles and stability policies are
runtime settings.

The following block loads the state once, explicitly reapplies the model and
runtime settings, and evaluates the three ordinary-soft precision controls.
Set `GL638_CASES=all` to run all 13 comparisons, or give comma-separated case IDs
from the JSON. It uses existing `set process` and raw `inspect -m --graph-id 0`
commands. All writes go to the replay output directory; the state is read-only.
It creates a minimal card containing only the recorded `default_runtime_settings`,
with no command blocks or generation settings. Loading that card replaces the
session's entire default runtime object; omitted fields take code defaults.
`set process ... defaults` then replaces the integrand's runtime settings, so
saved process overrides, selectors, observables, or model overrides cannot leak
into the replay. Each case resets all profile and rotation switches explicitly.

```bash
export GL638_CASES=soft22_z,soft22_euler,soft22_arb
python3 - <<'PY'
import hashlib
import json
import os
from decimal import Decimal
from pathlib import Path
import shlex
import subprocess

data = json.loads(Path("docs/research/gl638/replay-points.json").read_text())
requested = os.environ.get("GL638_CASES", "soft22_z,soft22_euler,soft22_arb")
known = {case["id"] for case in data["cases"]}
selected = known if requested == "all" else set(requested.split(","))
assert selected and selected <= known, f"Unknown case IDs: {selected - known}"
output = Path(os.environ["GL638_REPLAY_DIR"]).resolve()
output.mkdir(parents=True, exist_ok=True)
runtime_card = output / "runtime-only.toml"
runtime_card.write_text(data["runtime_card_toml"])
commands = [data["physics"]["model_command"],
            "set process -p epem_a_tth -i NNLO defaults"]
for case in data["cases"]:
    if case["id"] not in selected:
        continue
    settings = {**data["common_runtime_kv"], **data["profiles"][case["profile"]],
                "stability.rotation_axis": data["rotations"][case["rotation"]]}
    commands.append("set process -p epem_a_tth -i NNLO kv " + " ".join(
        shlex.quote(key + "=" + json.dumps(value, separators=(",", ":")))
        for key, value in settings.items()))
    # Round-trip every original f64; expand exponents for negative CLI arguments.
    momenta = data["points"][case["point"]]["parent_momenta_GeV"]
    point = [format(Decimal(format(value, ".17g")), "f")
             for vector in momenta for value in vector]
    inspect = ["inspect", "-p", "epem_a_tth", "-i", "NNLO", "-m",
               "--graph-id", "0", "--json-output", str(output / (case["id"] + ".json"))]
    if case["precision"] == "arb":
        inspect.append("-f")
    commands.append(shlex.join(inspect + ["--point"] + point))
binary = Path(os.environ["GL638_BIN"]).resolve()
argv = [str(binary), "--state-folder", str(Path(os.environ["GL638_STATE"]).resolve()),
        "--read-only-state", "--no-save-state", "-l", "info", "-L", "off", "-p", "none",
        str(runtime_card), "run", "-c", "; ".join(commands)]
(output / "invocation.json").write_text(json.dumps({
    "argv": argv, "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest()
}, indent=2) + "\n")
subprocess.run(argv, check=True)
PY
```

The momentum arrays contain four spatial vectors in generation-LMB edge order
`[3,4,7,10]`, in GeV **before** each cut's own LU rescaling. No orientation selector,
sampling Jacobian, or LMB channel weight is supplied. The A/P threshold's native
parent `[3,6,7,10]` is a different chart; do not substitute those coordinates here.
`-f` selects the configured Arb level, currently 1000-bit arithmetic, despite the
legacy CLI help describing this flag as f128.

| Cases | Recorded observation and limit of interpretation |
|---|---|
| `soft22_z`, `soft22_euler`, `soft22_arb` | CTs off, double-soft radius λ=0.01 GeV. Historical z rotation accepts Double Im ≈+3.098e−26; Euler rejects Double and gives Quad Im ≈−4.799e−31, agreeing with Arb here. The current card already uses Euler; the explicit z override recreates the old diagnostic policy. |
| `ap_equal_normal`, `ap_equal_arb` | Equal A/P roots: `is_nan=true` and a sanitized zero total in both modes. Historical records incorrectly say Stable; the current nonfinite check reports Unstable with no accuracy estimate. The evaluation failure remains unresolved; zero is not a physical result. |
| `ap_nearby_normal`, `ap_nearby_arb` | ηP=10⁻⁶ GeV and ηA/ηP=1.0001: historical z-based normal evaluation accepts Quad, but Im differs from Arb by a factor ≈1120. ηP is a threshold energy offset, not a soft radius. |
| `boundary_001_*`, `boundary_0001_*` | ν=0, λ=0.01/0.001 GeV. Baseline Im grows roughly 10×, hard-sliver Im roughly 1000×. The opt-in smooth sliver restores roughly 10× growth on these boundary rays. Baseline/smooth Arb values remain finite but miss the requested rotation tolerance. |

The `baseline`, `hard_sliver`, and `smooth_sliver` profiles use Gaussian width
`1.0`, fixed-Q widths, and sliver widths `10.0`, `0.1`, and `0.1`, respectively;
only the last enables `smooth_sliver`. They leave the integrated-CT profile at
its default. The local envelope is distinct from the LU `h_function` (power 3)
and the integrated-CT radial profile. The JSON preserves observations from both
the earlier `ce79…` binary and the later `530930…` binary, with complete hashes.

Inspect `evaluation.evaluation_metadata.is_nan`, the entire `stability_results`
history, the final status, and per-cut event weights. A successful command exit
or Stable label alone does not establish a finite or accurate result. Compare Re
and Im separately; norm agreement can hide errors in the smaller component.
Finite Unstable outputs must remain flagged. Two boundary points do not establish
a global power law, integrability, finite variance, or a complete GL638 cure.
