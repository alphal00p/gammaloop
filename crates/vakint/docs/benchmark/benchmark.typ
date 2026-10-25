// Render from the checked measurement record; never edit timing values here.
// Visual reference: https://github.com/alphal00p/ninjars/blob/main/docs/benchmark/benchmark.pdf
#let report = json("data/results.json")
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
#let five = if five-id == none { none } else { report.five_loop_runs.at(five-id) }
#let primary(id) = if five == none { none } else { five.primary.at(id) }
#let summary-time(value) = if value == none { [—] } else {
  [#timing(value.median_seconds) #text(7.5pt)[(#timing(value.min_seconds)–#timing(value.max_seconds))]]
}
#let primary-cell(id, variant) = {
  let row = primary(id)
  if row == none or row.status != "measured" {
    [Pending]
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
    #if five == none [Historical acceptance · performance measurements pending] else [Retained-master comparisons · functional evidence reported separately]
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
    columns: (1.8fr, 1fr, 1fr),
    inset: (x: 10pt, y: 11pt),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y == 0 { navy } else if x > 0 { rgb("#fff3dc") } else { pale },
    table.header(
      text(fill: white, weight: "bold")[Retained-master workload],
      text(fill: white, weight: "bold")[No persistent cache],
      text(fill: white, weight: "bold")[With persistent cache],
    ),
    ..for case in pending.cases {
      if case.id == "held-out-suite" {
        (case.label, if five == none or five.held_out == none { [Pending] } else { [Individual rows below] },
         if five == none or five.held_out == none { [Pending] } else { [Individual rows below] })
      } else {
        (case.label, primary-cell(case.id, "uncached"), primary-cell(case.id, "cached"))
      }
    },
  )
}

#v(0.35em)
#block(fill: pale, inset: 9pt, radius: 4pt, width: 100%)[
  #set text(9pt)
  *Eligibility and reading guide.* Primary cells require catalogue closure and
  exact result/guard agreement, a completed strict cache, and at least three matched fresh processes per
  mode in one accepted serialized window, with fixed workers and affinity.
  Bold values are reduction seconds; small values give the full observed range.
  Multipliers are _cached ÷ uncached_, using the same colour scale as page 1.
  Functional runs, including the exploratory q0 cohort, never supply these medians.
  Cache preparation and held-out performance retain separate evidence.
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


#if five != none {
  pagebreak()
  set text(font: "Libertinus Serif", size: 9.5pt, fill: black)
  align(center)[#text(20pt, weight: "bold", fill: navy)[Five-loop evidence and timing ledger]]
  v(0.35em)
  [Selected report: #raw(five-id). Raw reports and their hashed launcher/compiler
  receipts remain in JSON. Qualified performance additionally binds complete
  result/guard files, unknown-master inventories and accepted window receipts.
  A successful retained expression alone does not establish catalogue closure
  or numerical evaluability.]
  v(0.4em)
  text(13pt, weight: "bold", fill: navy)[Primary workloads]
  for case in pending.cases.filter(case => case.id != "held-out-suite") {
    let row = primary(case.id)
    v(0.3em)
    [*#case.label:* #row.status.]
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
  text(13pt, weight: "bold", fill: navy)[Strict cache preparation]
  if five.cache == none {
    [ Pending: no completed final strict-cache receipt selected.]
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
      parbreak()
      [#suite.rows.len() distinct input rows · #suite.identity.repetitions fresh processes per mode.
      Times below are saved application minus owner preparation, not pure IBP.
      Ratio is cached ÷ uncached; individual case ranges are retained.]
      if suite.at("cache_requalification", default: none) != none {
        parbreak()
        text(8pt)[Consumer source #raw(suite.identity.build_git_revision.slice(0, 12));
        original cache producer #raw(five.cache.raw_report.build_identity.build.git_revision.slice(0, 12)).
        Accepted requalification: #raw(suite.cache_requalification.evidence.sha256.slice(0, 12)).
        Cache, census, D40 data and Native semantic binding are unchanged.]
      }
      v(0.3em)
      table(columns: (1.5fr, 1fr, 1fr, 0.55fr, 0.8fr), inset: 5pt, stroke: 0.4pt + luma(215),
        fill: (x, y) => if y == 0 { pale } else { none },
        table.header([*Input*], [*Uncached s*], [*Cached s*], [*Ratio*], [*Unknowns / hits*]),
        ..suite.rows.map(row => {
          let u = row.variants.uncached.metrics.ibp_excluding_owner_seconds
          let c = row.variants.cached.metrics.ibp_excluding_owner_seconds
          let r = if u.median_seconds > 0 { c.median_seconds / u.median_seconds } else { none }
          (text(8pt, row.id), summary-time(u), summary-time(c),
           table.cell(fill: if r != none and r > 0 { ratio-fill(r) } else { none },
             if r == none { [—] } else { [×#multiplier(r)] }),
           [#row.catalogue_closure.unknown_count / #row.requested_root_hits + #row.intermediate_hits])
        }).flatten(),
      )
      v(0.3em)
      [Hits show requested-root + intermediate hits summed over cached samples;
      unknowns are the exact combined final inventory. Catalogue closure is independent
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
