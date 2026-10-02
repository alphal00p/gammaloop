= Massless Figure B1: matched complexity inventory

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<massless-figure-b1-matched-complexity-inventory>
The previous campaign does not isolate a soft-versus-ordinary-UV generation ratio for this graph. The short ordinary Figure B1 control uses a different, smaller graph, while the massless fixture compares one localized orientation with a complete explicit source sum. The matched baseline must keep the massless three-loop graph, subtraction forest, construction route, and selected orientation inventory fixed while changing its IR prescriptions to MUV.

== Exact input and controls
<exact-input-and-controls>
- Card: `/common/dev/gammaloop/lcnbr/tests/resources/run_cards/paper_figure_b1_double_triangle_soft_ir.toml`.
- Graph: `tests/resources/graphs/paper_figure_b1_double_triangle_soft_ir.dot`: one graph, three loops, eight internal edges and two external photon edges. The off-shell photon indices are saturated by the explicitly supplied metric projector. Loop coordinates are attached to e6, e5, e8, in that order.
- Model `sm-default`; `e_cm=300`; external momenta `(300,0,0,0)` and dependent; helicities `[1,-1]`; `m_uv=mu_r=20`; `SingleParametric` evaluator.
- Local-only ThreeD integrand, `hedge_poset`, UV subtraction enabled, integrated addbacks disabled, thresholds disabled, tropical generation disabled, `inner_products=true`, `add_marker=false`, iterative orientation optimization disabled, one generation core, evaluator compilation disabled in captured defaults.
- IR override 1: external PDGs `[-21,21]`, internal `[1]`, selecting the central massless-quark gluon bubble, degree two.
- IR override 2: external `[-22,22]`, internal `[1,21]`, selecting its positive-degree containing photon component. All other default schemes remain MUV.
- The card selects `orientation_delta(0,0,-1,-1,-1,1,-1,-1,1,-1)` via pattern `(0,0,-,-,-,+,-,-,+,-)`. `local_ct_cli` clears that filter when explicit summation is requested (`tests/tests/uv.rs:2532`). Thus the original localized/explcit comparison changes the amount of requested orientation work as well as its representation.
- The two-loop ordinary card `paper_figure_b1_double_triangle.toml` has five internal edges and no inserted massless bubble. Its timings are not a soft-disabled baseline for the three-loop graph.

== Exact supported U baseline
<exact-supported-u-baseline>
The ordinary baseline retains the same three-loop graph and sets all three defaults to MUV while clearing both IR overrides. This is precisely the existing scalar-spectacles ordinary control pattern at `tests/tests/uv.rs:4911`:

```json
{
  "global": {
    "generation": {
      "uv": {
        "renormalization_prescription": {
          "log_divergent": "MUV",
          "massive_power_divergent": "MUV",
          "massless_power_divergent": "MUV",
          "overrides": []
        }
      }
    }
  }
}
```

The existing CLI accepts this as `set global file /tmp/<case>/muv.json; run generate`. Keep the same explicit/localized flag, local4D flag, and generation orientation pattern for the paired H and U runs. The `uv.softct` boolean is not a reliable switch for this purpose: bounded source search found its declarations/defaults/test assignments but no production read controlling H; actual dispatch uses `ApproximationType::IR`.

For isolation, use an absolute graph-import path in a scratch copy of the card, a fresh `/tmp` state, `-n` to suppress automatic saving, and an explicit `save state --path /tmp/<case>/state`. After successful generation the existing `save uv-forest /tmp/<case>/forest -p paper_figure_b1_double_triangle_soft_ir -i paper_figure_b1_double_triangle_soft_ir` writes a structure-only DOT when `--computed` is omitted. It requires generated state and generation settings history (`uv/export.rs:131`, `processes/process.rs:908`); it cannot recover the missing forest merely from the original test\'s initial state files. The execution agent owns the bounded run plan under `/tmp/soft-ct-complexity-2026-09-14`; this inventory launches no generation.

== H versus U forest topology
<h-versus-u-forest-topology>
With identical graph, cuts, `subtract_uv=true`, and neither scheme Unsubtracted, changing IR to MUV retains the same classified spinney filters and DODs:

+ `uv/uv_graph.rs:145–161` drops a connected region only when DOD is negative or its selected scheme is Unsubtracted.
+ `uv/spinney.rs:66–121` uses scheme as a stored label and diagnostic, not as an input to compatible-sub-LMB selection or component filtering.
+ Spinney equality and inclusion order depend solely on the subgraph (`spinney.rs:128–150`).
+ Disconnected retention requires the same independently retained factors (`uv_graph.rs:173–211`); disconnected aggregate labels are neutral MUV in both cases.
+ Wood construction and trace unfolding use those subgraph filters, inclusion edges and factorized joins (`hedge_poset.rs:200–361`). Scheme does not alter this topology.

Consequently the controlled pair should have identical structural forest nodes/edges, with H versus U labels and different attached local expressions. Internal U/S/US sectors or provenance histories are different from structural forest-node counts. A same-forest comparison remains necessary to measure how much those expressions grow; the algebraic formula alone does not supply a measured multiplicative factor.

There is no automatic promotion of further ancestors to IR. `uv/orchestrator.rs:109–139` rejects a positive-degree MUV parent containing a positive-degree IR child and instructs the caller to assign IR to that parent. The card makes that containing-component assignment explicitly. The local4D boundary also rejects positive-degree MUV with soft ancestry (`local_4d.rs:2381`). Clearing both overrides removes this condition rather than silently changing other schemes.

== Available and absent campaign evidence
<available-and-absent-campaign-evidence>
The original massless testcase retained five initial-state directories: localized H, localized bare, failed explicit H, and two parent-orientation probes. None contains `amp.bin`, `integrand.bin`, or `generation_summary.json`; they are not saved generated evaluators. The successful localized H/bare directories retain their UV profile JSON. Both resolve and compare all 196 summed and 196 single-orientation records. No explicit profile or projected directory exists.

The shared JSONL log is under the first orientation-probe directory. It contains 33 `Computed global numerator` messages and two rounded generation tables. These repeated projector summaries are not forest counts or final expression sizes. Localized H reports 0.600 s expression, 1.85 s Spenso build, 3.20 s Symbolica build and 522.41 MiB sampled RAM; localized bare reports 0.128/0.014/0.003 s. There is no ordinary-U three-loop result in this campaign. The explicit H failure occurred during generation after 105.445 minutes and 313.962 GiB peak RSS.

The fixture\'s independent cheap parent-CFF probe proves more than one valid nonzero direction signature, audits its deterministic first signature, preserves all source-map multiplicities, and verifies the complete source inventory after localized generation (`uv.rs:5191–5307`, `5369–5433`). It does not persist the numerical source-map count. Its structure-only forest assertion checks at least four nodes and noncomputed metadata, but does not write an exact forest inventory. No exact source-map/forest-node/input-atom byte count can be recovered from the saved compact logs.

An existing compact source-map diagnostic is available at `cff/generation.rs:1191`: tags `#generation,#cff,#profile`, stage `raw_cff_generation`, fields `native_source_maps`, `graph`, `elapsed_ms`, `source_reconstruction_ms`, `preparation_ms`, `native_generation_ms`, `postprocessing_ms`. `native_source_maps` counts generated source entries, including repeated direction signatures. Earlier proposal events at line 899 also report map counts and must not be added to final selected counts. The original campaign did not enable this CFF tag.

A graph-only bound is at most 256 fully directed sign strings for eight internal edges, or 6561 if zero/Undirected is allowed, with external slots fixed zero. This loose bound is not a measured count and does not bound source-map multiplicity, atom size, or work per map; it is unsuitable as an explanation for the observed memory cost.

== History and limits
<history-and-limits>
Read-only jj history shows both the exact three-loop graph and its card first introduced by `yssunplynlpy` (`be7c95063af1`, implement local soft counterterms). Neither exists at pre-soft base `c08b0d9f1605`. The card\'s one-orientation filter, explicit containing-IR rationale, and atom-size-ceiling diagnostic comments were present at that first introduction. There is no identical-graph pre-soft timing in this source history or the selected campaign. This does not prove nobody ran it elsewhere; it means the available record cannot establish that the same requested graph/forest/full-sum workload previously generated cheaply.

The current first differing boundary in an H/U pair is the local subtraction operation attached to otherwise identical forest topology. For the original localized/explicit pair, orientation selection differs earlier and must be controlled first. No graph numerator was expanded, no source/test code changed, and no generation was launched for this inventory.
