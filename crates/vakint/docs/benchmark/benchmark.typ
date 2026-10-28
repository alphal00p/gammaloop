// Render from the checked measurement record; never edit timing values here.
// Visual reference: https://github.com/alphal00p/ninjars/blob/main/docs/benchmark/benchmark.pdf
#let report = json("data/results.json")
#let delivery-id = report.at("selected_delivery_checkpoint", default: none)
#let selected = report.at("selected_reduction_run", default: none)
#let data = if selected == none { report } else { report.reduction_runs.at(selected) }
#assert(data.validation_passed, message: "Only validated benchmark results may be published")
#let navy = rgb("#1f3b63")
#let pale = rgb("#eef3f9")

#set document(title: "Vakint reduction benchmark")
#set page(paper: "a4", margin: (x: 1.6cm, y: 1.7cm), numbering: "1 / 1")
#set text(font: "Libertinus Serif", size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.62em)
#show raw: set text(font: "DejaVu Sans Mono", size: 8pt)
#show link: set text(fill: rgb("#1f5fa8"))

#let fixed(x, d) = {
  let parts = str(calc.round(x, digits: d)).split(".")
  if d == 0 { return parts.first() }
  let fraction = if parts.len() > 1 { parts.at(1) } else { "" }
  parts.first() + "." + fraction + "0" * (d - fraction.len())
}
#let timing(x) = fixed(x, if x >= 10 { 1 } else if x >= 1 { 2 } else if x >= 0.1 { 3 } else { 4 })
#let multiplier(x) = fixed(x, if x < 0.1 { 3 } else { 2 })

// Same logarithmic green–white–red scale as the NinjaRS reference report.
#let palette = gradient.linear(
  rgb("#5aae61"), rgb("#d9f0d3"), rgb("#ffffff"), rgb("#fddbc7"), rgb("#d6604d"),
)
#let ratio-fill(r) = {
  let t = calc.clamp(0.5 + calc.log(r, base: 2) / 4, 0, 1)
  palette.sample(t * 100%)
}
#let comparison(family, kind) = {
  let row = data.summary.find(row => row.family == family and row.kind == kind)
  let native = row.cached.median_seconds
  let reference = row.fmft.median_seconds
  let ratio = native / reference
  table.cell(fill: ratio-fill(ratio), {
    set par(leading: 0.35em)
    text(size: 17pt, weight: "bold")[#timing(native) #text(size: 10pt)[s]]
    linebreak()
    text(size: 9pt, fill: luma(65))[FMFT #timing(reference) s]
    h(7pt)
    text(size: 10pt, weight: "semibold")[×#multiplier(ratio)]
  })
}

#align(center)[
  #v(0.4em)
  #text(23pt, weight: "bold", fill: navy)[Vakint reduction benchmark]
  #v(-0.2em)
  #text(12pt, fill: luma(80))[Cached RustRed versus the FORM / FMFT reduction route]
  #v(0.2em)
  #text(9.5pt, fill: luma(110))[
    #data.completed_utc.split("T").first() · #data.protocol.cpu_model · CPU affinity #str(data.protocol.cpu_affinity)
  ]
]

#v(0.6em)
#block(fill: pale, inset: 10pt, radius: 4pt, width: 100%)[
  #set text(9.5pt)
  *Reading guide.* Bold values are RustRed times; FMFT times appear below.
  Multipliers are _RustRed ÷ FMFT_: green means faster, white parity, red slower.
  Both receive identical four-loop scalar inputs. The first physical denominator
  has power 2 (base) or 4 (held out); all others have power 1, auxiliaries 0.
]

#v(0.5em)
#text(14pt, weight: "bold", fill: navy)[Reduction time]
#v(0.3em)
Median seconds over #data.protocol.repetitions fresh processes per backend and case.

#v(0.35em)
#{
  set text(font: "DejaVu Sans", size: 9pt)
  set par(justify: false)
  show table.cell.where(y: 0): set text(fill: white, weight: "bold")
  table(
    columns: (0.7fr, 1.5fr, 1.5fr),
    align: (x, y) => if x == 0 { left + horizon } else { center + horizon },
    inset: (x: 12pt, y: 12pt),
    stroke: (x, y) => (
      top: if y <= 1 { none } else { 0.4pt + luma(225) },
      left: if x >= 2 and y > 0 { 0.4pt + luma(225) } else { none },
    ),
    fill: (x, y) => if y == 0 { navy } else { none },
    table.header([Topology], [Base · power 2], [Held out · power 4]),
    ..for family in ("H", "FG", "BMW", "X") {
      (
        text(size: 13pt, weight: "semibold", family),
        comparison(family, "base"),
        comparison(family, "held-out"),
      )
    },
  )
}
#v(0.45em)
#align(center)[
  #set text(8pt)
  #box(width: 9cm, height: 0.3cm, fill: palette, radius: 2pt)
  #linebreak()
  #box(width: 9cm)[≤ ×1/4 #h(1fr) ×1/2 #h(1fr) ×1 #h(1fr) ×2 #h(1fr) ≥ ×4]
]

#v(0.6em)
#block(stroke: 0.6pt + luma(200), inset: 10pt, radius: 4pt, width: 100%)[
  #set text(9.5pt)
  *Method and correctness.*
  #if data.protocol.at("configured_profile", default: none) != none [
    Release build: optimization level #data.protocol.configured_profile.opt_level,
    #data.protocol.configured_profile.lto LTO;
  ] else [
    Build profile: #data.protocol.build_profile (caller-declared);
    optimization/LTO settings are not inferred.
  ]
  #(data.protocol.rustc.split(" (").first() + "; " +
    data.protocol.form_version.split("\n").first().split(" (").first() + ".")
  Both retain masters in FMFT/MATAD normalization, with identical symbolic mass.
  Exact rational parity passed all eight cases at unit mass and renormalization
  scale, without truncating epsilon or using numerical tolerances.

  Base cases hit prepared roots without applying rules. Held-out cases reuse
  cached descendants after root misses. Timing covers first public reduction
  calls, including FORM version checks, launch, I/O and parsing. Data loading,
  cache attachment, master evaluation and validation are excluded. Fresh
  processes run sequentially, alternating backends, with the recorded CPU
  affinity on a shared host.
  Filesystem pages are warm; this does not measure cold startup.
]

#v(0.5em)
#set text(8.5pt, fill: luma(100))
#if selected == none [Historical measurements, unchanged.] else [Selected validated run: #raw(selected).]
Source: #raw(data.protocol.source_commit.slice(0, 12)).
#if data.protocol.source_dirty [Modified source tree.]
Data: #raw(data.protocol.data_version).
Samples and provenance: #raw("data/results.json").
Refresh instructions: #raw("../../AGENTS.md").

#pagebreak()
#let ingress = report.five_loop_dot_ingress_acceptance
#let pending = report.five_loop_performance_status
#let five-id = report.at("selected_five_loop_run", default: none)
#let tier-selections = report.at("selected_five_loop_runs", default: none)
#let five-selections = if tier-selections == none { ((id: five-id, tier: none),) } else {
  ("heavy", "light").filter(tier => tier in tier-selections).map(tier => (
    id: tier-selections.at(tier), tier: tier,
  ))
}
#let five-run(selection) = if selection.id == none { none } else { report.five_loop_runs.at(selection.id) }
#let cache-schema(run) = if run == none or run.cache == none { none } else { run.cache.raw_report.at("schema", default: none) }
#let tier-label(tier) = if tier == "heavy" { [HEAVY · native cache] } else { [LIGHT · selected gzip] }
#for selection in five-selections.filter(item => item.tier != none) {
  let run = five-run(selection)
  assert(run.cache_tier == selection.tier)
  assert(cache-schema(run) == if selection.tier == "heavy" {
    "vakint.cumulative-native-cache-report.v1"
  } else { "vakint.strict-light-cache-export-validation.v1" })
}
#let cumulative-cache = five-selections.any(selection => cache-schema(five-run(selection)) in (
  "vakint.cumulative-native-cache-report.v1", "vakint.strict-light-cache-export-validation.v1",
))
#let physical-id = report.at("selected_physical_capture_run", default: none)
#let physical = if physical-id == none { none } else { report.physical_capture_runs.at(physical-id) }
#let physical-deferred(id) = physical != none and (
  id in physical.at("physical_scope", default: (deferred: ())).deferred
)
#let primary(run, id) = if run == none { none } else { run.primary.at(id) }
#let physical-scope-note(run) = {
  let authorization = run.at("physical_scope_authorization", default: none)
  if authorization != none {
    let scope = run.physical_scope
    let labels(ids) = ids.map(id => pending.cases.find(case => case.id == id).label).join(", ")
    block(fill: pale, inset: 8pt, radius: 4pt, width: 100%)[
      #set text(8.5pt)
      *Required physical scope:* #labels(scope.required).
      #if scope.deferred.len() > 0 [*Deferred, unmeasured:* #labels(scope.deferred).]
      #authorization.raw_receipt.reason
      Scope authorization: #raw(authorization.evidence.sha256.slice(0, 12)).
      Historical records remain unchanged.
    ]
  }
}
#let summary-time(value) = if value == none { [—] } else {
  [#timing(value.median_seconds) #text(7.5pt)[(#timing(value.min_seconds)–#timing(value.max_seconds))]]
}
#let primary-cell(run, id, variant) = {
  let row = primary(run, id)
  if row == none or row.status != "measured" {
    if physical-deferred(id) or (row != none and row.status == "deferred") { [Deferred] } else { [Pending] }
    if row != none and row.functional_records.len() > 0 {
      linebreak()
      text(8pt)[Functional evidence only]
    }
  } else {
    let metric = row.at(variant).metrics.reduction_seconds
    let reference = row.uncached.metrics.reduction_seconds.median_seconds
    let ratio = if reference > 0 { row.cached.metrics.reduction_seconds.median_seconds / reference } else { none }
    table.cell(fill: if variant == "cached" and ratio != none and ratio > 0 { ratio-fill(ratio) } else { pale }, {
      text(15pt, weight: "bold")[#timing(metric.median_seconds) s]
      linebreak()
      text(8pt)[#timing(metric.min_seconds)–#timing(metric.max_seconds) s]
      linebreak()
      text(7.5pt)[#if row.at("input_mode", default: "fresh") == "prepared" [Prepared input] else [Fresh preparation]]
      if variant == "cached" and ratio != none {
        linebreak()
        text(9pt, weight: "bold")[×#multiplier(ratio)]
      }
    })
  }
}
#assert(ingress.validation_passed, message: "Five-loop ingress checks must pass")
#set text(font: "Libertinus Serif", size: 10pt, fill: black)

#align(center)[
  #text(22pt, weight: "bold", fill: navy)[Five-loop validation and benchmarks]
  #v(0.2em)
  #text(11.5pt, fill: luma(80))[
    #if five-selections.all(selection => selection.id == none) [Historical acceptance · performance measurements pending] else [Retained-master comparisons · functional evidence reported separately]
  ]
]

#v(0.55em)
#text(13pt, weight: "bold", fill: navy)[Historical DOT ingress acceptance]
#v(0.2em)
#ingress.completed_utc.split("T").first() · one sequential debug-build pass.
This checks graph definitions; it does not time IBP reduction, master evaluation,
package loading or cache preparation.

#v(0.35em)
#{
  set text(font: "DejaVu Sans", size: 9pt)
  set par(justify: false)
  table(
    columns: (2.1fr, 1fr),
    inset: (x: 11pt, y: 8pt),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y == 0 { navy } else if calc.rem(y, 2) == 0 { pale } else { none },
    table.header(
      text(fill: white, weight: "bold")[Validation workload],
      text(fill: white, weight: "bold")[Recorded result],
    ),
    [Master graphs / dotted / with numerator],
    [#ingress.graphs.masters_checked / #ingress.graphs.dotted_masters / #ingress.graphs.numerator_masters],
    [Exact routing, full numerators, native ingress and DOT round trips], [Passed],
    [Complete graph-validation pass],
    text(size: 16pt, weight: "bold", fill: navy)[#timing(ingress.graphs.elapsed_seconds) s],
    [Per graph: minimum / median / maximum],
    [#timing(ingress.graphs.minimum_seconds) / #timing(ingress.graphs.median_seconds) / #timing(ingress.graphs.maximum_seconds) s],
    [Numerical master values / prepared 5L cache in this run], [Absent / absent],
  )
}

#v(0.3em)
#set text(9pt)
Physical powers, signed indices, mass terms, full numerators, loop-basis slots
(including zero) and integration Jacobians were checked exactly. These checks
establish the stored graph definitions, without claiming that the catalogue is
an independent master basis. One observed run supports no statistical speedup
claim. Earlier interrupted checks are not a timing baseline.

#v(0.5em)
#text(13pt, weight: "bold", fill: navy)[Five-loop reduction comparison]
#v(0.25em)
#{
  set text(font: "DejaVu Sans", size: 9pt)
  set par(justify: false)
  table(
    columns: (1.8fr, ..five-selections.map(_ => (1fr, 1fr)).flatten()),
    inset: (x: 10pt, y: if tier-selections == none { 11pt } else { 8pt }),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y == 0 { navy } else if x > 0 { rgb("#fff3dc") } else { pale },
    table.header(
      text(fill: white, weight: "bold")[Retained-master workload],
      ..for selection in five-selections {
        (
          text(fill: white, weight: "bold")[#if selection.tier == none [No persistent cache] else [#upper(selection.tier)#linebreak()Uncached]],
          text(fill: white, weight: "bold")[#if selection.tier == none [With persistent cache] else [#upper(selection.tier)#linebreak()Cached]],
        )
      },
    ),
    ..for case in pending.cases {
      (case.label, ..for selection in five-selections {
        let five = five-run(selection)
        if case.id == "held-out-suite" {
          (if five == none or five.held_out == none { [Pending] } else { [Individual rows below] },
           if five == none or five.held_out == none { [Pending] } else { [Individual rows below] })
        } else {
          (primary-cell(five, case.id, "uncached"), primary-cell(five, case.id, "cached"))
        }
      })
    },
  )
}

#{
  let shown-scopes = ()
  for selection in five-selections.filter(item => item.id != none) {
    let five = five-run(selection)
    let authorization = five.at("physical_scope_authorization", default: none)
    let digest = if authorization == none { none } else { authorization.evidence.sha256 }
    if digest == none or digest not in shown-scopes {
      v(0.3em)
      physical-scope-note(five)
      shown-scopes.push(digest)
    }
  }
}

#v(0.35em)
#block(fill: pale, inset: 9pt, radius: 4pt, width: 100%)[
  #set text(9pt)
  *Eligibility and reading guide.* Primary cells require catalogue closure and
  exact result/guard agreement,
  #if cumulative-cache [a validated snapshot of strictly completed records,] else [a completed strict cache,]
  and at least three matched fresh processes per
  mode in one accepted serialized window, with fixed workers and affinity.
  Bold values are reduction seconds; small values give the full observed range.
  Multipliers are _cached ÷ uncached_, using the same colour scale as page 1.
  Functional runs, including the exploratory q0 cohort, never supply these medians.
  Cache preparation and held-out performance retain separate evidence.
  #if delivery-id != none [The user waived new uncached q0 repetitions; its
  historical reference and completed functional checks appear on the delivery page.]
  #if tier-selections != none [Each tier retains its own matched uncached baseline;
  no baseline or sample is inferred from the other tier.]
]

#v(0.35em)
#set text(8pt, fill: luma(100))
Historical ingress: source #raw(ingress.protocol.source_commit.slice(0, 12))
(modified); RustRed #raw(ingress.protocol.rustred_commit.slice(0, 12));
Symbolica #raw(ingress.protocol.symbolica_commit.slice(0, 12)).
Exact patch, executable, data and per-graph digests remain in JSON.

#v(0.2em)
Historical compatibility record, separate from selected performance evidence:
RustRed #raw(pending.validated_compatibility_cohort.rustred_revision.slice(0, 12));
Symbolica #raw(pending.validated_compatibility_cohort.symbolica_numerica_graphica_revision.slice(0, 12)).
Published upgrade: Symbolica
#raw(pending.pending_upgrade.symbolica_numerica_graphica_revision.slice(0, 12)),
FeynKit #raw(pending.pending_upgrade.feynkit_owner_source_revision.slice(0, 12));
acceptance was pending in that historical record. Full revisions: #raw("data/results.json").


#for selection in five-selections.filter(item => item.id != none) {
  let five-id = selection.id
  let five = five-run(selection)
  let cumulative-cache = cache-schema(five) == "vakint.cumulative-native-cache-report.v1"
  let light-cache = cache-schema(five) == "vakint.strict-light-cache-export-validation.v1"
  pagebreak()
  set text(font: "Libertinus Serif", size: 9.5pt, fill: black)
  align(center)[#text(20pt, weight: "bold", fill: navy)[Five-loop evidence and timing ledger]]
  if selection.tier != none { align(center)[#tier-label(selection.tier)] }
  v(0.35em)
  [Selected report: #raw(five-id). Raw reports and their hashed launcher/compiler
  receipts remain in JSON. Qualified performance additionally binds complete
  result/guard files, unknown-master inventories and accepted window receipts.
  A successful retained expression alone does not establish catalogue closure
  or numerical evaluability.]
  v(0.4em)
  text(13pt, weight: "bold", fill: navy)[Primary workloads]
  for case in pending.cases.filter(case => case.id != "held-out-suite") {
    let row = primary(five, case.id)
    v(0.3em)
    [*#case.label:* #if row.status == "deferred" or (row.status != "measured" and physical-deferred(case.id)) [Deferred; unmeasured in this scope.] else [#row.status.]]
    if row.functional_records.len() > 0 {
      let functional = row.functional_records.filter(record => not record.performance_eligible)
      if functional.len() > 0 {
        [ #functional.len() functional record(s); no performance aggregate.]
        for record in functional {
          let original = record.raw_report
          let source = original.at("source_commit", default: none)
          let native = record.artifact.at("native", default: none)
          parbreak()
          text(8pt)[#record.status · report #raw(record.evidence.record.sha256.slice(0, 12))
            #if source != none [· source #raw(source.slice(0, 12))]
            #if native != none [· original Native #raw(native.base_revision.slice(0, 12))
              · Symbolica #raw(record.artifact.symbolica_revision.slice(0, 12))]
            #if record.catalogue_closure != none [· unknown masters: #record.catalogue_closure.unknown_count]
          ]
        }
      }
    }
    if row.status == "measured" {
      parbreak()
      let records = row.functional_records.filter(record => record.performance_eligible)
      let evidence = records.first().performance_evidence
      [#row.uncached.metrics.reduction_seconds.samples_seconds.len() samples per mode;
      exact result/guard parity and catalogue closure passed.
      CPU #evidence.cpu_affinity · source #raw(evidence.stamps.VAKINT_CACHE_BUILD_GIT_REVISION.slice(0, 12)).]
      let reused = row.at("input_mode", default: "fresh") == "prepared"
      if reused {
        parbreak()
        text(8.5pt)[*Prepared-input scope.* Both modes load the same authenticated scalar fixture
        #raw(row.prepared_input.manifest_sha256.slice(0, 12)). Original generator preparation:
        #timing(row.prepared_input.manifest.preparation_seconds) s (historical, excluded from these timings).
        Generator acceptance: #raw(row.generator_acceptance.sha256.slice(0, 12)).]
      }
      v(0.2em)
      let ledger = (
        if reused { ("Prepared input loading", "metrics", "input_loading_seconds") }
          else { ("Fresh graph preparation", "metrics", "preparation_seconds") },
        ("Data/program loading", "metrics", "data_loading_seconds"),
        ("Manual 5L cache import", "metrics", "cache_loading_seconds"),
        ("Parent transport / scalar lowering", "stages", "parent_transport_and_scalar_lowering"),
        ("Deferred owner preparation (included)", "metrics", "owner_preparation_seconds"),
        ("Saved application (inclusive)", "stages", "saved_rule_application"),
        ("Complete public reduction", "metrics", "reduction_seconds"),
        (if reused { "Input load → result serialization" } else { "End to end (before report exit)" }, "metrics", "end_to_end_seconds"),
      )
      table(columns: (2.2fr, 1fr, 1fr), inset: 4pt, stroke: 0.4pt + luma(215),
        fill: (x, y) => if y == 0 { pale } else { none },
        table.header([*Seconds: median (range)*], [*Uncached*], [*Cached*]),
        ..ledger.map(item => (item.at(0),
          summary-time(row.uncached.at(item.at(1)).at(item.at(2))),
          summary-time(row.cached.at(item.at(1)).at(item.at(2))))).flatten(),
      )
    }
  }
  v(0.4em)
  text(13pt, weight: "bold", fill: navy)[#if light-cache [LIGHT selected-cache export] else if cumulative-cache [Cumulative native cache] else [Finite selected-cache preparation]]
  if five.cache == none {
    [ Pending: no completed final strict-cache receipt selected.]
  } else if light-cache {
    let cache = five.cache
    let exported = cache.raw_report
    let source = cache.source_producer
    let validation = cache.export_validation
    parbreak()
    [Strict records selected from the preserved HEAVY collection. The export,
    native reimport and source-selection admission retain their own evidence.]
    v(0.25em)
    table(columns: (1fr, 1.3fr, 1fr), inset: (x: 8pt, y: 10pt),
      stroke: 0.4pt + luma(215),
      fill: (x, y) => if y == 0 { navy } else { pale },
      table.header(
        text(fill: white, weight: "bold")[Selected unique records],
        text(fill: white, weight: "bold")[Compressed bytes],
        text(fill: white, weight: "bold")[Publication limit],
      ),
      text(14pt, weight: "bold")[#exported.selection.unique_keys],
      [#exported.output.cache.bytes (#fixed(exported.output.cache.bytes / 1048576, 3) MiB)],
      [95,000,000 bytes],
    )
    v(0.25em)
    table(columns: (1fr, 1fr, 1fr, 1fr, 1fr), inset: 5pt, stroke: 0.4pt + luma(215),
      fill: (x, y) => if y == 0 { pale } else { none },
      table.header([*Selection s*], [*Preparation s*], [*Export s*], [*Reimport s*], [*Comparison s*]),
      ..("selection", "preparation", "export", "reimport", "comparison").map(phase => {
        let value = cache.production_seconds.at(phase, default: none)
        if value == none { [—] } else { [#timing(value)] }
      }),
    )
    v(0.2em)
    text(8.5pt)[A dash means unmeasured. Preparation loads the producer's saved application;
    export sums encoding and writing only when both are measured. Comparison checks
    the reimported native records. Selection outside the helper remains unmeasured.
    Costs belong to the native-validation producer cohort, separately from consumer
    loading and reduction. File size and record count establish neither duration
    nor completed physical work.]
    parbreak()
    text(8pt)[Cache SHA-256 #raw(exported.output.cache.sha256).
    Export #raw(validation.export_validation.sha256.slice(0, 12));
    strict native reimport #raw(validation.native_validation.sha256.slice(0, 12));
    selection admission #raw(validation.selection_admission.sha256.slice(0, 12));
    metadata #raw(exported.output.metadata.sha256.slice(0, 12)).
    Exporter source #raw(exported.compiled_cohort.main_revision.slice(0, 12));
    Native #raw(exported.compiled_cohort.native_revision.slice(0, 12));
    Symbolica #raw(exported.compiled_cohort.kernel_revision.slice(0, 12)).
    Original HEAVY producer #raw(source.evidence.sha256.slice(0, 12));
    source index #raw(exported.source.cache_native_index_identity.sha256.slice(0, 12));
    source dataset #raw(exported.source.data_manifest_sha256.slice(0, 12)).
    These source identities remain distinct from each fresh consumer qualification.]
  } else if cumulative-cache {
    let cache = five.cache
    let inventory = cache.raw_report.inventory
    let deferred = cache.at("physical_scope", default: (deferred: ())).deferred
    assert(cache.status in ("preserved-committed-prefix", "validated-committed-snapshot"))
    assert(not cache.performance_eligible)
    parbreak()
    [*#if cache.status == "preserved-committed-prefix" [Committed prefixes preserved.] else [Committed native snapshot validated.]*
    Preservation and native validation are separate from completion of the source
    reductions. Each source keeps its original receipt and provenance.]
    v(0.25em)
    table(columns: (1fr, 1fr, 1.2fr), inset: (x: 8pt, y: 8pt),
      stroke: 0.4pt + luma(215),
      fill: (x, y) => if y == 0 { navy } else { pale },
      table.header(
        text(fill: white, weight: "bold")[Stored references],
        text(fill: white, weight: "bold")[Native chunks],
        text(fill: white, weight: "bold")[Compressed bytes],
      ),
      text(12pt, weight: "bold")[#inventory.stored_record_references],
      text(12pt, weight: "bold")[#inventory.physical_chunks],
      [#inventory.compressed_bytes (#fixed(inventory.compressed_bytes / 1048576, 3) MiB)],
    )
    v(0.2em)
    text(8.5pt)[Stored references may overlap or be shadowed; they are not a unique-key count.
    No preparation duration or performance ratio is inferred from this inventory.]
    v(0.25em)
    table(columns: (2fr, 1fr), inset: 5pt, stroke: 0.4pt + luma(215),
      fill: (x, y) => if x == 0 { pale } else { none },
      ..pending.cases.filter(case => case.id != "held-out-suite").map(case => (
        case.label,
        if case.id in deferred { [Deferred] }
          else if case.id in cache.completed_physical_examples { [Completion validated] } else { [Pending] },
      )).flatten(),
      [All retained source runs complete],
      if cache.all_source_runs_completed { [Yes] } else { [Not established] },
    )
    v(0.2em)
    text(8pt)[#cache.sources.len() original source receipt(s).
    Required physical scope complete: #if cache.physical_examples_complete [yes.] else [not established.]
    Snapshot identity #raw(cache.raw_report.cache_native_index_identity.sha256.slice(0, 12));
    cache report #raw(cache.evidence.sha256.slice(0, 12));
    validation #raw(cache.validation.evidence.sha256.slice(0, 12)).
    Prefix-only sources remain prefix-only. Qualified timing comparisons require
    their own accepted consumer, exact-parity and execution evidence.]
  } else {
    let build = five.cache.raw_report
    let file = build.preparation.options.single_file
    parbreak()
    [Producer #raw(build.build_identity.build.git_revision.slice(0, 12)) · original phase: #raw(build.phase).
    #if five.cache.at("final_selection", default: none) != none [
      Retained completed candidate; selection recorded separately as
      #raw(five.cache.final_selection.evidence.sha256.slice(0, 12)).
    ]
    #box[#build.preparation.completed_roots / #build.preparation.total_roots] selected roots complete;
    #build.preparation.newly_stored_intermediates intermediate records; strict Reject policy.
    Compressed size: #fixed(build.compressed_bytes / 1048576, 3) MiB;
    decoded: #file.identity.decoded_bytes bytes (#fixed(file.identity.decoded_bytes / 1048576, 3) MiB).
    #if build.phase == "pilot" [Original selected-pilot producer-driver elapsed] else [Producer-driver elapsed]:
    #raw(str(build.elapsed_seconds)) s,
    including loading and owner setup, not pure IBP or overall campaign time.
    Cache SHA-256: #raw(build.cache_sha256). Training roots, capture diagnostics and
    source/proof bindings remain in the unmodified build receipt.]
    if five.cache.at("final_selection", default: none) != none {
      let selection = five.cache.final_selection
      let bridge = selection.cache_requalification
      parbreak()
      text(8pt)[Accepted native compatibility bridge #raw(bridge.evidence.sha256.slice(0, 12)) ·
      consumer cohort receipt #raw(bridge.raw_receipt.consumer.compiled_cohort.sha256.slice(0, 12)) ·
      consumer workload #raw(bridge.raw_receipt.consumer.final_workload.sha256.slice(0, 12)).]
      if bridge.raw_receipt.at("scope", default: (:)).at("native_compatibility_only", default: false) {
        text(8pt)[ This does not establish a fresh target-suite pass, later-cohort compatibility or performance.]
      }
      let trial = selection.raw_receipt.at("discarded_trial", default: none)
      if trial != none {
        parbreak()
        text(8pt)[Separate rejected trial: #trial.strict_complete_roots / #trial.requested_roots roots complete;
        #raw(str(trial.whole_attempt_seconds)) s whole attempt; no final gzip candidate.]
      }
    }
  }
  if not cumulative-cache and not light-cache {
    parbreak()
    text(8.5pt)[This finite selection retains its original producer and acceptance.
    The cumulative cache from physical-example reductions is a separate record;
    #if physical == none [accepted producer/reuse completion is pending.] else [see the physical capture supplement.]]
  }
  v(0.4em)
  block(fill: pale, inset: 9pt, radius: 4pt)[
    *Timing boundaries.* Stages overlap; never sum them. Owner preparation is
    included in saved-application/reduction times. #raw("ibp_excluding_owner_seconds")
    means saved application minus owner preparation, including transport,
    combination and terminal substitution. Pure IBP, terminal-only, isolated gzip decode and primary
    owner-subtracted times remain unmeasured. Prepared-input E2E includes current
    input/data/cache loading and public reduction; generator preparation stays
    historical. Fresh/prepared scopes never mix. Initial loading uses Sparse;
    representation reload belongs to reduction. Missing events are “—”, never zero.
    Full records remain in JSON.
  ]
  let reuse = five.at("warm_reuse", default: none)
  if reuse != none and reuse.functional_records.len() > 0 {
    if reuse.status == "measured-paired-reuse" { pagebreak() } else { v(0.6em) }
    text(if reuse.status == "measured-paired-reuse" { 20pt } else { 13pt }, weight: "bold", fill: navy)[Genuinely connected five-loop control]
    v(0.4em)
    [*Supplement: same-engine warm reuse.* This scalar ladder-pinch control is
    separate from the physical quadratic/q0 workloads and from held-out
    persistent-cache comparisons. Each fresh process makes one first public API
    call and then repeats it on the same engine, settings and input.]
    v(0.4em)
    if reuse.status == "measured-paired-reuse" {
      let proof = reuse.functional_records.filter(record => record.performance_eligible).first()
      [Three fresh process pairs; full values, ordered guards and family identities
      agree. Median seconds and observed ranges below cover whole public API calls.]
      v(0.3em)
      table(columns: (1fr, 1fr), inset: 12pt, stroke: 0.4pt + luma(215),
        table.header(
          table.cell(fill: navy, text(fill: white)[First public call]),
          table.cell(fill: navy, text(fill: white)[Repeated call · same engine]),
        ),
        table.cell(fill: pale, align(center, summary-time(reuse.first_call))),
        table.cell(fill: if reuse.warm_over_first != none and reuse.warm_over_first > 0 { ratio-fill(reuse.warm_over_first) } else { pale },
          align(center, summary-time(reuse.warm_call))),
      )
      if reuse.warm_over_first != none {
        parbreak()
        [Warm median ÷ first-call median: *×#multiplier(reuse.warm_over_first)*.]
      }
      v(0.4em)
      text(8pt)[Source #raw(proof.performance_evidence.stamps.VAKINT_CACHE_BUILD_GIT_REVISION.slice(0, 12));
      executable #raw(proof.performance_evidence.binary_sha256.slice(0, 12));
      data #raw(proof.input.data_manifests.at("manifest.json").sha256.slice(0, 12)).
      CPU #proof.performance_evidence.cpu_affinity; one fixed worker.]
    } else {
      [#reuse.functional_records.len() completed functional pair(s).
      These observations are non-isolated; no median, range or ratio is reported.
      Three qualified fresh pairs in one accepted window are required for timing aggregates.]
      for record in reuse.functional_records {
        v(0.3em)
        table(columns: (1fr, 1fr), inset: 6pt, stroke: 0.4pt + luma(215),
          table.header(
            table.cell(fill: navy, text(fill: white)[First public call · observed seconds]),
            table.cell(fill: navy, text(fill: white)[Repeated call · observed seconds]),
          ),
          table.cell(fill: pale, align(center, raw(str(record.calls.cold.whole_public_application_seconds)))),
          table.cell(fill: pale, align(center, raw(str(record.calls.warm.whole_public_application_seconds)))),
        )
        v(0.2em)
        text(8pt)[Functional source #raw(record.input.build_git_revision.slice(0, 12)) ·
        executable #raw(record.input.executable.sha256.slice(0, 12)) ·
        data #raw(record.input.data_manifests.at("manifest.json").sha256.slice(0, 12)).]
      }
    }
    v(0.4em)
    block(fill: pale, inset: 10pt, radius: 4pt)[
      *Scope.* Default prepared 1L–4L caches are attached; no persistent 5L cache
      is attached before or after either call. On-demand owner preparation is
      included inside the API timer. Loading, input preparation, serialization,
      cache snapshots and parity checks are excluded. Reuse may involve several
      kinds of in-memory state; this is not a selected persistent-cache speedup.
      Pure IBP, owner-only and terminal-only durations remain unmeasured.
      These values never enter the cold-process or physical-workload medians.
      Family identity parity establishes neither catalogue closure nor numerical
      master availability.
    ]
  }
  for (kind, label) in (("held_out", "Protected held-out inputs"), ("tuning", "Training-selection inputs (separate)")) {
    let suite = five.at(kind)
    if suite != none {
      pagebreak()
      text(18pt, weight: "bold", fill: navy)[#label]
      if selection.tier != none { parbreak(); tier-label(selection.tier) }
      parbreak()
      [#suite.rows.len() distinct input rows · #suite.identity.repetitions fresh processes per mode.
      Times below are saved application minus owner preparation, not pure IBP.
      Ratio is cached ÷ uncached; individual case ranges are retained.]
      if suite.at("cache_requalification", default: none) != none {
        parbreak()
        text(8pt)[Consumer source #raw(suite.identity.build_git_revision.slice(0, 12));
        #if light-cache [
          LIGHT export #raw(five.cache.export_validation.export_validation.sha256.slice(0, 12));
          original HEAVY producer #raw(five.cache.source_producer.evidence.sha256.slice(0, 12));
          consumer dataset #raw(suite.identity.data_manifest_sha256.slice(0, 12)).
        ] else if cumulative-cache [
          cumulative cache report #raw(five.cache.evidence.sha256.slice(0, 12)).
          Original producer identities remain in its #five.cache.sources.len() source receipts.
        ] else [
          original cache producer #raw(five.cache.raw_report.build_identity.build.git_revision.slice(0, 12)).
        ]
        Accepted requalification: #raw(suite.cache_requalification.evidence.sha256.slice(0, 12)).
        #if light-cache [
          The actual LIGHT file, native identity and source-selection proof are bound
          to this consumer qualification; HEAVY timings are not transferred.
        ] else if cumulative-cache [
          Cache indexes, census, data-manifest hash and native binding are checked against
          the accepted consumer qualification.
        ] else [Cache, census, D40 data and Native semantic binding are unchanged.]]
      }
      v(0.3em)
      table(columns: (1.5fr, 1fr, 1fr, 0.55fr, 0.8fr), inset: 5pt, stroke: 0.4pt + luma(215),
        fill: (x, y) => if y == 0 { pale } else { none },
        table.header([*Input*], [*Uncached s*], [*Cached s*], [*Ratio*], [*Unknowns / R + U*]),
        ..suite.rows.map(row => {
          let u = row.variants.uncached.metrics.ibp_excluding_owner_seconds
          let c = row.variants.cached.metrics.ibp_excluding_owner_seconds
          let r = if u.median_seconds > 0 { c.median_seconds / u.median_seconds } else { none }
          (text(8pt, row.id), summary-time(u), summary-time(c),
           table.cell(fill: if r != none and r > 0 { ratio-fill(r) } else { none },
             if r == none { [—] } else { [×#multiplier(r)] }),
           [#row.catalogue_closure.unknown_count / #row.requested_root_hits + #row.intermediate_hits
             #linebreak()
             #text(7pt)[#if row.requested_root_hits > 0 [Some root reuse.]
               #if row.intermediate_hits > 0 [Other scope unknown.]
               #if row.requested_root_hits == 0 and row.intermediate_hits == 0 [No backing hits.]]])
        }).flatten(),
      )
      v(0.3em)
      [R counts requested Route-record hits; U counts other hits, summed over cached
      samples. U may include requested-owner Apply hits as well as descendants;
      R = 0 does not prove intermediate-only reuse. No finer scope is inferred from
      a revision or the collector's legacy label. Unknowns are the exact combined
      final inventory. Catalogue closure is independent
      of retained-result success. No median across different physical inputs is reported.]
      v(0.4em)
      text(13pt, weight: "bold", fill: navy)[Stage ledger]
      [ Each cell gives uncached → cached median seconds. Ranges and every sample are in JSON.]
      v(0.25em)
      {
        set text(size: 7.5pt)
        let compact(row, section, field) = {
          let u = row.variants.uncached.at(section).at(field)
          let c = row.variants.cached.at(section).at(field)
          [#if u == none [—] else [#timing(u.median_seconds)] → #if c == none [—] else [#timing(c.median_seconds)]]
        }
        table(columns: (1.5fr, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr), inset: 4pt, stroke: 0.4pt + luma(215),
          fill: (x, y) => if y == 0 { pale } else { none },
          table.header([*Input*], [*Data load*], [*Cache load*], [*Lowering*], [*Owner prep*], [*Saved API*], [*E2E*]),
          ..suite.rows.map(row => (row.id,
            compact(row, "metrics", "data_loading_seconds"), compact(row, "metrics", "cache_loading_seconds"),
            compact(row, "stages", "parent_transport_and_scalar_lowering"), compact(row, "metrics", "owner_preparation_seconds"),
            compact(row, "metrics", "saved_application_including_owner_seconds"), compact(row, "metrics", "end_to_end_seconds"))).flatten(),
        )
      }
      v(0.3em)
      [Suite end-to-end timing spans worker spawn through exit, unlike the primary
      before-exit timer. Graph preparation, pure IBP and terminal-only phases are
      not separately measured here. Ratios never combine these two timing scopes.]
    }
  }
}

#if physical != none {
pagebreak()
set text(font: "Libertinus Serif", size: 10pt, fill: black)
align(center)[
  #text(21pt, weight: "bold", fill: navy)[Zero-momentum physical reduction]
  #v(0.2em)
  #text(11pt, fill: luma(80))[Completed cache capture · fresh reuse of the same input]
]
v(0.6em)
  assert(physical.status == "functional-seen-physical-reuse")
  assert(not physical.performance_eligible)
  assert(physical.seen_not_heldout)
  physical-scope-note(physical)
  v(0.3em)
  [*One completed fresh seen-input pair; functional evidence only.*
  Complete retained output, canonical public nonzero-condition file, unknown-master inventory
  and family map agree. The fresh reuse wrote no new records and preserved the
  native snapshot identity.]
  v(0.5em)
  table(columns: (1.2fr, 1fr, 1fr), inset: (x: 10pt, y: 13pt),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y == 0 { navy } else { pale },
    table.header(
      text(fill: white, weight: "bold")[Observed seconds],
      text(fill: white, weight: "bold")[Capture producer],
      text(fill: white, weight: "bold")[Fresh seen reuse],
    ),
    [Complete public reduction],
    text(16pt, weight: "bold")[#timing(physical.observed_seconds.producer_reduction) s],
    text(16pt, weight: "bold")[#timing(physical.observed_seconds.cached_reduction) s],
    [Supervised collector wall time],
    [#timing(physical.observed_seconds.producer_wall) s],
    [#timing(physical.observed_seconds.reuse_wall) s],
  )
  v(0.3em)
  text(8.5pt)[These are raw single observations, with possible concurrent historical
  workers. No median, speedup claim or performance colour is assigned. Producer
  reduction includes capture overhead; pure IBP time is unmeasured.
  Both processes load the same prepared scalar input; graph preparation is excluded.]
  v(0.4em)
  let hits = physical.cache.cached_persistent_statistics
  let native = physical.raw_reuse.report.saved_application_cache.native
  let closure = physical.raw_producer.report.catalogue_closure
  let grouped(value) = {
    let digits = str(value)
    let first = calc.rem(digits.len(), 3)
    let parts = if first == 0 { () } else { (digits.slice(0, first),) }
    (parts + range(first, digits.len(), step: 3).map(i => digits.slice(i, i + 3))).join(",")
  }
  table(columns: (1fr, 1fr, 1fr), inset: (x: 8pt, y: 10pt),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y == 0 { navy } else { pale },
    table.header(
      text(fill: white, weight: "bold")[Strict records captured],
      text(fill: white, weight: "bold")[Seen Route hits / misses],
      text(fill: white, weight: "bold")[Distinct masters / unknowns],
    ),
    text(15pt, weight: "bold")[#grouped(physical.cache.producer_accepted_records)],
    text(15pt, weight: "bold")[#grouped(hits.route_hits) / #hits.misses],
    text(15pt, weight: "bold")[#closure.master_count / #closure.unknown_count],
  )
  v(0.2em)
  text(8.5pt)[Reuse applied #native.rule_applications rules, with #hits.apply_hits Apply hits,
  #hits.intermediate_hits intermediate hits and #physical.cache.newly_captured_records new records.
  Captured records are not a deduplicated-key census or the final HEAVY file size.]
  v(0.45em)
  let captured = physical.raw_producer.artifact
  text(9pt)[Producer source #raw(captured.source_head.slice(0, 12));
  actual Native archive #raw(captured.elf.native_tree.slice(0, 12));
  ELF #raw(captured.elf.sha256.slice(0, 12));
  frozen source #raw(captured.elf.source_freeze_sha256.slice(0, 12)).
  Optimization 3 and fat LTO are bound by the compiler receipt.
  This original producer cohort is separate from the current application/dependency build.
  Selected evidence: #raw(physical-id).
  Complete source, executable, input, cache and result bindings remain in
  #raw("data/results.json").]
v(0.65em)
block(fill: pale, inset: 11pt, radius: 4pt, width: 100%)[
  *Cache production policy.* Preserve useful strict roots and completed
  descendants generated by the physical examples, together with every durable
  prefix and immutable base. The finite training driver alone does not define
  the final cumulative cache. Native packs have no 50 MiB publication target or
  95 MiB aggregate cap; report actual storage and operational limits.
]
v(0.45em)
[Seen-input reuse is separate from the independently frozen held-out suite and
from same-engine RAM repeats. Protected held-out runs read an immutable snapshot;
newly captured work goes to a separate delta. Completed descendants do not prove
completion of an unfinished parent, stronger campaign coverage or numerical
master availability.]

v(0.4em)
set text(8.5pt, fill: luma(100))
[Append accepted evidence with #raw("merge_physical_capture_results.py").
Refresh and validation contract: #raw("maintenance.typ").]
}

#let engineering-id = report.at("selected_engineering_optimization_run", default: none)
#if engineering-id != none {
  let engineering = report.engineering_optimization_runs.at(engineering-id)
  let before = engineering.observations.before
  let after = engineering.observations.after_off
  let diagnostic = engineering.observations.after_on
  let old = before.records.benchmark.raw_receipt
  let new = after.records.benchmark.raw_receipt
  let observed = diagnostic.records.benchmark.raw_receipt
  let verified = after.records.verification.raw_receipt
  assert(engineering.schema == "vakint.q0-engineering-comparison.v1")
  assert(not engineering.performance_eligible)
  assert(not engineering.matched_serialized_repetitions)
  assert(not engineering.cache_only_speedup_claim)
  pagebreak()
  set text(font: "Libertinus Serif", size: 10pt, fill: black)
  align(center)[
    #text(21pt, weight: "bold", fill: navy)[q0 coefficient assembly]
    #v(0.2em)
    #text(11pt, fill: luma(80))[Verified engineering observations · 28 October 2026]
  ]
  v(0.65em)
  [Coefficient-local retained-result assembly replaces the two whole-expression
  exact-combination passes in the earlier observation. All
  #verified.full_scalar_coefficient_count master-product coefficients agree exactly
  with the authenticated original native baseline across
  #verified.master_identity_union_count distinct masters; all
  #str(verified.guard_count) canonical public guards and the full family map agree.
  Master values and numerical samples were never evaluated.]
  v(0.55em)
  table(columns: (1.7fr, 1fr, 1fr), inset: (x: 9pt, y: 9pt),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y == 0 { navy } else { pale },
    table.header(
      text(fill: white, weight: "bold")[Observed seconds],
      text(fill: white, weight: "bold")[Earlier · #old.source_commit.slice(0, 7)],
      text(fill: white, weight: "bold")[New · #new.source_commit.slice(0, 7)\ Diagnostics OFF],
    ),
    [Complete public reduction],
    [#timing(before.public_reduction_seconds) s], [#timing(after.public_reduction_seconds) s],
    [Prepared input through result serialization],
    [#timing(before.prepared_input_to_serialization_seconds) s],
    [#timing(after.prepared_input_to_serialization_seconds) s],
    [Supervised process wall time],
    [#timing(before.supervised_wall_seconds) s], [#timing(after.supervised_wall_seconds) s],
    [Retained exact combination\*],
    [#timing(old.stages.retained_outer_exact_combination.sum() + old.stages.retained_catalogue_exact_combination.sum()) s],
    [#timing(new.stages.retained_catalogue_exact_combination.sum()) s],
    [Complete native saved call†],
    [#timing(old.stages.saved_rule_application.sum()) s],
    [#timing(new.stages.saved_rule_application.sum()) s],
  )
  v(0.35em)
  [The matching prepared-input-to-serialization observations differ by
  ×#multiplier(engineering.observed_before_over_after_off.prepared_input_to_serialization)
  (about 5.5×). *This is a single engineering comparison:* source/cache cohorts,
  instrumentation and concurrent background work differ. It supplies neither
  final three-sample medians nor a cache-only speedup claim. Both runs reuse the
  same prepared scalar input; original graph preparation is excluded.]
  v(0.35em)
  text(8.5pt)[\*Earlier: outer #timing(old.stages.retained_outer_exact_combination.sum()) s
  plus catalogue #timing(old.stages.retained_catalogue_exact_combination.sum()) s.
  New: coefficient-local catalogue #timing(new.stages.retained_catalogue_exact_combination.sum()) s;
  the outer phase is absent, not a measured zero.
  †The saved call includes native cache lookup/decode, transport, weighted exact
  assembly and terminal substitution. Its observed duration increased. Pure IBP
  and terminal-only durations remain unmeasured.]
  v(0.55em)
  [*Diagnostics ON confirms reuse.* The separate fresh instrumented process took
  #timing(diagnostic.prepared_input_to_serialization_seconds) s through serialization,
  #timing(diagnostic.supervised_wall_seconds) s supervised wall time, and
  #timing(diagnostic.public_reduction_seconds) s for the public call. It returned
  the identical #str(observed.result_bytes)-byte result and complete guards as OFF,
  with #observed.saved_application_cache.persistent.route_hits Route hits,
  #observed.saved_application_cache.persistent.misses misses,
  #observed.saved_application_cache.native.rule_applications new rule applications
  and #observed.saved_application_cache.persistent.captured_records captured records.
  OFF/ON controls hit diagnostics; both attach a persistent cache.]
  v(0.5em)
  text(8.5pt)[New source #raw(new.compiled_cohort.main_revision.slice(0, 12)) ·
  Native #raw(new.compiled_cohort.native_revision.slice(0, 12)) ·
  kernel #raw(new.compiled_cohort.kernel_revision.slice(0, 12)) ·
  #new.compiled_cohort.rustc · declared profile #new.compiled_cohort.profile.
  Original compiled cohorts and verifier-host limits remain in the raw receipts;
  these checks do not establish final Main-host acceptance.]
  v(0.35em)
  text(8pt, fill: luma(100))[Selected engineering record: #raw(engineering-id).
  Hashed source index: #raw("provenance/q0_engineering_evidence_index.json").
  Raw benchmark, execution, exact comparison and all product proofs are appended
  in #raw("data/results.json"). Earlier selections and all four-loop records remain unchanged.]
}

#let files-id = report.at("selected_cache_file_engineering_run", default: none)
#if files-id != none {
  let study = report.engineering_optimization_runs.at(files-id)
  let heavy = study.records.heavy.raw_receipt
  let new = heavy.full_actual_benchmark_metadata
  let old = report.engineering_optimization_runs.at(study.baseline_engineering_run)
    .observations.after_off.records.benchmark.raw_receipt
  let light = study.records.light.raw_receipt
  let a = light.pilots.at(0)
  let b = light.pilots.at(1)
  let export-a = study.records.candidate_a.raw_receipt
  let export-b = study.records.candidate_b.raw_receipt
  assert(study.schema == "vakint.cache-file-engineering-comparison.v1")
  assert(not study.performance_eligible and not study.final_light_selection)
  pagebreak()
  set text(font: "Libertinus Serif", size: 9pt, fill: black)
  set par(leading: 0.45em)
  align(center)[#text(19pt, weight: "bold", fill: navy)[Cache file engineering observations]]
  v(0.35em)
  [*HEAVY · one complete cold q0 process.* Same executable, source
  #raw(new.source_commit.slice(0, 7)), Native #raw(new.compiled_cohort.native_revision.slice(0, 7))
  and kernel #raw(new.compiled_cohort.kernel_revision.slice(0, 7)) as the earlier INDEX run.
  All 58 coefficients, 55 masters, 146678 canonical public guards and the full
  result agree. Capture was enabled with delta reads disabled:
  #heavy.Native_work.persistent.route_hits Route hits, #heavy.Native_work.persistent.misses misses,
  #heavy.Native_work.native.rule_applications new rules and
  #heavy.Native_work.captured_records captured records.]
  v(0.3em)
  table(columns: (1.9fr, 1fr, 1fr), inset: (x: 7pt, y: 5pt), stroke: 0.4pt + luma(215),
    fill: (x, y) => if y == 0 { navy } else { pale },
    table.header(text(fill: white)[Observed seconds], text(fill: white)[Earlier INDEX · OFF], text(fill: white)[HEAVY single FILE]),
    [Data / native cache loading],
    [#timing(old.data_loading_seconds) / #timing(old.cache_loading_seconds)],
    [#timing(new.data_loading_seconds) / #timing(new.cache_loading_seconds)],
    [Complete public call], [#timing(old.reduction_seconds)], [#timing(new.reduction_seconds)],
    [Collector end to end], [#timing(old.end_to_end_seconds)], [#timing(new.end_to_end_seconds)],
    [Inclusive controller elapsed],
    [#timing(report.engineering_optimization_runs.at(study.baseline_engineering_run).observations.after_off.supervised_wall_seconds)],
    [#timing(heavy.clock_summary_seconds.inclusive_controller_elapsed_after_preflight)],
    [Native saved call / final coefficient combination],
    [#timing(old.stages.saved_rule_application.sum()) / #timing(old.stages.retained_catalogue_exact_combination.sum())],
    [#timing(new.stages.saved_rule_application.sum()) / #timing(new.stages.retained_catalogue_exact_combination.sum())],
  )
  v(0.3em)
  [Public call #multiplier(-study.heavy_observed_change_percent.public_call)% lower;
  cold collector end to end #multiplier(study.heavy_observed_change_percent.collector_E2E)% higher.
  The lower native and CAS observations cannot isolate a cache-format benefit:
  concurrent work and capture policy differ. No causal format advantage, overall
  cold speedup, warm result or serialized three-run median is established.]
  v(0.25em)
  text(8pt)[Native FILE loading covers authentication, decode and indexing. The collector
  made two admitted SHA scans of the 96 GB file: its end-to-end timer includes the first
  scan and flush/snapshot, excludes close and final SHA; the controller includes both
  scans and cleanup. Their difference is not pure SHA time. SHA-only, decoder-only
  and pure IBP times are unmeasured. This legacy CLI used #raw("file_identity=None")
  with Native framing/binding checks and before/after SHA, distinct from typed R3 acceptance.]
  v(0.5em)
  [*LIGHT · one frozen tuning case, two utility pilots.* A contains 54 Route records;
  B adds 2048 Apply records. Both fully match the original INDEX canonical result,
  guards and five families (68 masters, zero unknowns). Neither used held-out inputs.]
  v(0.3em)
  table(columns: (1.9fr, 1fr, 1fr), inset: (x: 7pt, y: 5pt), stroke: 0.4pt + luma(215),
    fill: (x, y) => if y == 0 { navy } else { pale },
    table.header(text(fill: white)[Observed utility], text(fill: white)[A · 54 Route], text(fill: white)[B · 54 Route + 2048 Apply]),
    [Compressed bytes / native records],
    [#str(a.cache_file.bytes) / #export-a.native_counts.selected_records],
    [#str(b.cache_file.bytes) / #export-b.native_counts.selected_records],
    [Public / saved call (seconds)],
    [#timing(a.timing.evaluation_seconds) / #timing(a.timing.saved_application_including_owner_seconds)],
    [#timing(b.timing.evaluation_seconds) / #timing(b.timing.saved_application_including_owner_seconds)],
    [Data / cache loading (seconds)],
    [#timing(a.timing.data_loading_seconds) / #timing(a.timing.cache_loading_seconds)],
    [#timing(b.timing.data_loading_seconds) / #timing(b.timing.cache_loading_seconds)],
    [New rules / Route hits / Apply hits],
    [#a.actual_cache_statistics.native.rule_applications / #a.actual_cache_statistics.persistent.route_hits / #a.actual_cache_statistics.persistent.apply_hits],
    [#b.actual_cache_statistics.native.rule_applications / #b.actual_cache_statistics.persistent.route_hits / #b.actual_cache_statistics.persistent.apply_hits],
    [Closed capture records / directory bytes],
    [#a.closed_running_capture.report.committed_records / #str(a.closed_running_capture.closed_directory_bytes)],
    [#b.closed_running_capture.report.committed_records / #str(b.closed_running_capture.closed_directory_bytes)],
  )
  v(0.3em)
  [B used fewer rules but took longer. The earlier full INDEX observation for this
  same case took #timing(light.baseline.evaluation_seconds) s and
  #light.baseline.new_native_rules new rules (receipt #raw(light.prior_measurement.sha256.slice(0, 12))).
  These fresh pilots ran alongside unrelated physical work; they do not establish
  a performance winner. The later user-directed LIGHT freeze is recorded on the
  delivery page. All new completed capture nodes
  remain preserved; neither run read its delta. Public timings include capture overhead.]
  v(0.25em)
  text(8pt)[Both native exports passed full reimport/count/ordered-guard checks within
  95000000 compressed bytes. Their separate export helper measured encode+write
  #timing(export-a.phase_seconds.encode_selected_completed_cache + export-a.phase_seconds.write_selected_cache) /
  #timing(export-b.phase_seconds.encode_selected_completed_cache + export-b.phase_seconds.write_selected_cache) s,
  reimport #timing(export-a.phase_seconds.reimport_selected_cache) /
  #timing(export-b.phase_seconds.reimport_selected_cache) s and exact visitor comparison
  #timing(export-a.phase_seconds.compare_selected_native_visitor) /
  #timing(export-b.phase_seconds.compare_selected_native_visitor) s (A/B). Exporter identity
  receipts do not imply that the application supplied a typed Native file identity.
  Complete source cohorts, closed captures, raw phases and timer scopes remain in
  #raw("data/results.json"). Source hashes in
  #raw("provenance/cache_file_engineering_evidence_index.json").]
}

#if delivery-id != none {
  let delivery = report.delivery_checkpoints.at(delivery-id)
  let freeze = delivery.records.light_freeze.raw_receipt
  let monitor = delivery.records.monitor.raw_receipt
  let sources = delivery.records.heavy_sources.raw_receipt
  assert(delivery.schema == "vakint.q0-delivery-checkpoint.v1")
  assert(not delivery.performance_eligible and not delivery.cache_performance_winner)
  assert(delivery.new_uncached_repetitions_waived)
  pagebreak()
  set text(font: "Libertinus Serif", size: 10pt, fill: black)
  align(center)[
    #text(21pt, weight: "bold", fill: navy)[Zero-momentum delivery checkpoint]
    #v(0.2em)
    #text(11pt, fill: luma(80))[Frozen LIGHT cache · exact monitored reuse · 28 October 2026]
  ]
  v(0.55em)
  block(fill: pale, inset: 10pt, radius: 4pt, width: 100%)[
    *Scope and timing reference.* The required physical example is the complete
    zero-momentum four-rung fixture; the quadratic example is deferred.
    The user waived three new uncached repetitions. The completed historical
    uncached public reduction took #timing(delivery.historical_uncached_seconds) s
    (#fixed(delivery.historical_uncached_seconds / 3600, 1) h), in its original older cohort.
    It remains a historical reference, not a matched baseline for a new cache-only
    ratio. The earlier approximately 5.5× assembly observation concerns Vakint's
    coefficient-local source change; it is separate from cache-format performance.
  ]
  v(0.55em)
  text(13pt, weight: "bold", fill: navy)[LIGHT · accepted as sufficient]
  v(0.25em)
  table(columns: (1fr, 1fr, 1fr), inset: (x: 8pt, y: 11pt),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y == 0 { navy } else { rgb("#e8f1ea") },
    table.header(
      text(fill: white)[Strict native records],
      text(fill: white)[Compressed size],
      text(fill: white)[Decoded size],
    ),
    text(17pt, weight: "bold")[#freeze.selected_records],
    text(17pt, weight: "bold")[#fixed(freeze.compressed_bytes / 1000000, 2) MB],
    text(17pt, weight: "bold")[#fixed(freeze.decoded_bytes / 1000000, 2) MB],
  )
  v(0.25em)
  [Candidate B retains 54 Route and 2048 Apply records. Complete Native reimport
  verified exact coefficients and the original ordered cache-record guards.
  The user accepted this existing cache and ended selection/tuning; this does
  not declare B a performance winner. The earlier A/B observations and all
  completed capture prefixes remain preserved.]
  v(0.2em)
  block[
    #set text(8pt)
    #set par(justify: false)
    Exact lengths: #freeze.compressed_bytes compressed / #freeze.decoded_bytes decoded bytes.
    Compressed SHA-256: #raw(freeze.sha256). Publication and actual download/offline
    validation remain pending in this checkpoint.
  ]
  v(0.55em)
  text(13pt, weight: "bold", fill: navy)[Monitored HEAVY reuse · functional acceptance]
  v(0.25em)
  table(columns: (1fr, 1fr, 1fr), inset: (x: 8pt, y: 11pt),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y == 0 { navy } else { pale },
    table.header(
      text(fill: white)[Exact master identities / coefficients],
      text(fill: white)[Heartbeat freshness bound],
      text(fill: white)[Observed peak RAM],
    ),
    text(17pt, weight: "bold")[#monitor.science.master_identity_union_count / #monitor.science.exact_scalar_product_count],
    text(17pt, weight: "bold")[#fixed(monitor.monitor.maximum_conservative_freshness_bound_seconds, 3) s],
    text(17pt, weight: "bold")[#fixed(monitor.actual_OS_resources.worker_HWM_max_observed_bytes / 1000000000, 2) GB],
  )
  v(0.25em)
  [All #monitor.science.ordered_guard_Atom_count exported guard Atoms agree in order,
  with mass, couplings, normalization and family identities preserved.
  #monitor.monitor.sample_count OS/monitor samples observed no faults.
  The heartbeat covered loading, native work, exact combination, serialization
  and the final SHA scan; native scientific counters may age during opaque phases.
  All physical workers have exited. This instrumented check supplies no new
  performance median.]
  v(0.25em)
  text(8pt)[The peak is an observed OS high-water mark, not a bound on every
  instantaneous allocation. The 50 GB machine-available-memory guard stayed
  satisfied; minimum observed available RAM was
  #fixed(monitor.actual_OS_resources.machine_MemAvailable_min_observed_bytes / 1000000000, 1) GB.
  Public guard-vector equality is distinct from the original ordered guards
  inside each Native cache record.]
  v(0.5em)
  block(stroke: 0.5pt + luma(200), inset: 9pt, radius: 4pt, width: 100%)[
    *Final HEAVY export: pending.* All #sources.source_count precious sources are
    closed, with #str(sources.total_stored_record_references) stored record references
    before deduplication. Final unique count, compressed/decoded lengths,
    SHA-256 and export/reimport durations remain unmeasured. The previously
    validated 96 GB snapshot stays preserved; it is not relabeled as this final union.
  ]
  v(0.35em)
  text(8pt, fill: luma(100))[Evidence: #raw("provenance/q0_delivery_evidence_index.json").
  Exact raw receipts and user decisions are appended in #raw("data/results.json").
  No master values were evaluated: the 640 declared five-loop masters remain uncomputed.]
}
