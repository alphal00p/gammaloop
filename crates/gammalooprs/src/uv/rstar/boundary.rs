//! The opt-in remainder-boundary producer (`GAMMALOOP_REMAINDER_BOUNDARY=1|check`).
//!
//! Today one root-scope R* remainder reaches vakint as the expanded
//! `Σ_{i+j≤deg} P_i D_j`: the contracted numerator times the Taylor
//! coefficients of its denominator product, expanded into monomials and
//! grouped per topology. On the new path (beta_computation DESIGN.md v2) the
//! same remainder is handed over factored, as a vakint [`RemainderBoundary`]:
//!
//! * the denominator product alone is Taylor expanded ([`denominator_orders`]);
//!   its monomials give the power patterns, their vacuum frames
//!   (`TermFrame`) and their cut (`cut_frame`), with one composite routing;
//! * the numerator `P = C(σ(N))` is contracted without the denominators and
//!   split by one pass of a recursive grammar into blocks `S_g` of grade
//!   `g = (T, s, b)` (k·p structure, power of p², power of the reference mass
//!   squared), which are free of `p` and of the reference mass ([`Grammar`]);
//! * per pattern and grade, the small coefficient `R_{π,g}` carries the
//!   structure `t_g`, the mass compensation and the denominator orders,
//!   projected onto the vacuum tensors and multiplied by the scalar root
//!   normalisation.
//!
//! Nothing here expands a product `P_i × D_j` or a sum of loop dots. Every
//! block, coefficient and routing image is checked against the vakint
//! contract before the hand-over ([`verify_vakint_scalar`]).

use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    sync::LazyLock,
    time::Instant,
};

use eyre::{Result, WrapErr, eyre};
use linnet::half_edge::{involution::EdgeIndex, subgraph::SuBitGraph};
use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    structure::{
        representation::{Minkowski, RepName},
        slot::{SlotMatch, SlotMatcher},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    function,
    id::Replacement,
    symbol,
};
use vakint::{
    TensorReductionMethod, Vakint, VakintSettings,
    p_integral::{BoundaryOptions, BoundaryPattern, MassiveClosure, RemainderBoundary},
    symbols::S,
};

use crate::{
    utils::GS,
    uv::approx::{
        integrated::{TermFrame, TermRewrites, term_frame},
        streamed,
    },
};

use super::{
    TaylorExpansion,
    expression::RStarIntegrand,
    integrate::{VacuumCut, cut_frame, tensor_reduce_by_factors},
};

// ---------------------------------------------------------------------------
// Switches
// ---------------------------------------------------------------------------

/// `GAMMALOOP_REMAINDER_BOUNDARY`: today's path (`Off`), the boundary for
/// every admitted remainder (`On`), or both compared (`Check`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoundaryMode {
    Off,
    On,
    Check,
}

impl BoundaryMode {
    pub fn as_str(&self) -> &'static str {
        match self {
            BoundaryMode::Off => "off",
            BoundaryMode::On => "on",
            BoundaryMode::Check => "check",
        }
    }
}

fn parse_mode(value: Option<&str>) -> std::result::Result<BoundaryMode, String> {
    match value {
        None | Some("" | "0" | "off") => Ok(BoundaryMode::Off),
        Some("1" | "on") => Ok(BoundaryMode::On),
        Some("check") => Ok(BoundaryMode::Check),
        Some(other) => Err(format!(
            "GAMMALOOP_REMAINDER_BOUNDARY={other:?}: expected unset, 0, off, 1, on or check"
        )),
    }
}

static MODE: LazyLock<std::result::Result<BoundaryMode, String>> =
    LazyLock::new(|| match std::env::var("GAMMALOOP_REMAINDER_BOUNDARY") {
        Ok(value) => parse_mode(Some(&value)),
        Err(std::env::VarError::NotPresent) => parse_mode(None),
        Err(std::env::VarError::NotUnicode(value)) => Err(format!(
            "GAMMALOOP_REMAINDER_BOUNDARY={value:?} is not unicode"
        )),
    });

static STRICT: LazyLock<bool> = LazyLock::new(|| {
    std::env::var("GAMMALOOP_REMAINDER_BOUNDARY_STRICT").as_deref() == Ok("1")
});

#[cfg(test)]
thread_local! {
    static TEST_MODE: std::cell::Cell<Option<(BoundaryMode, bool)>> = const { std::cell::Cell::new(None) };
}

/// The mode of `GAMMALOOP_REMAINDER_BOUNDARY`, read once; a value other than
/// unset, empty, `0`, `off`, `1`, `on` or `check` is an error at first use.
pub fn remainder_boundary_mode() -> Result<BoundaryMode> {
    #[cfg(test)]
    if let Some((mode, _)) = TEST_MODE.with(|slot| slot.get()) {
        return Ok(mode);
    }
    MODE.clone().map_err(|error| eyre!(error))
}

/// `GAMMALOOP_REMAINDER_BOUNDARY_STRICT=1`: every decline after admission is
/// an error, and vakint pattern declines are errors too.
pub(crate) fn boundary_strict() -> bool {
    #[cfg(test)]
    if let Some((_, strict)) = TEST_MODE.with(|slot| slot.get()) {
        return strict;
    }
    *STRICT
}

/// Tests: run `f` with the mode (and strictness) set on this thread, never
/// through the environment.
#[cfg(test)]
pub(crate) fn with_boundary_mode<T>(mode: BoundaryMode, strict: bool, f: impl FnOnce() -> T) -> T {
    let previous = TEST_MODE.with(|slot| slot.replace(Some((mode, strict))));
    let result = f();
    TEST_MODE.with(|slot| slot.set(previous));
    result
}

/// `GAMMALOOP_BOUNDARY_FALLBACK_MAX_SUMMANDS` (default 1,000,000 summands of
/// `P`; 0 = no limit): above it a decline is an error in every mode and vakint
/// runs strict (DESIGN 10c).
pub(crate) fn fallback_max_summands() -> usize {
    static MAX: LazyLock<usize> = LazyLock::new(|| {
        match std::env::var("GAMMALOOP_BOUNDARY_FALLBACK_MAX_SUMMANDS")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
        {
            Some(0) => usize::MAX,
            Some(max) => max,
            None => 1_000_000,
        }
    });
    *MAX
}

/// The vakint environment gates that `lower_direct_blocks` does not apply
/// itself (DESIGN A10, vakint `evaluate_direct`): `VAKINT_LOWERING` unset or
/// `direct`, the fold not off, no `VAKINT_P_FOLD_CHECK`. Read once.
pub(crate) fn vakint_gates() -> std::result::Result<(), &'static str> {
    static GATES: LazyLock<std::result::Result<(), &'static str>> = LazyLock::new(|| {
        match std::env::var("VAKINT_LOWERING").as_deref() {
            Err(_) | Ok("direct") => {}
            Ok(_) => return Err("lowering"),
        }
        if matches!(
            std::env::var("VAKINT_P_FOLD").as_deref(),
            Ok("off" | "0" | "false")
        ) {
            return Err("fold_off");
        }
        if std::env::var_os("VAKINT_P_FOLD_CHECK").is_some() {
            return Err("fold_check");
        }
        Ok(())
    });
    *GATES
}

/// The settings part of the A10 gate: RustRed present and before Forcer, as
/// the direct reader requires (`evaluate_direct`).
pub(crate) fn settings_gate(settings: &VakintSettings) -> std::result::Result<(), &'static str> {
    let first = settings.evaluation_order.0.iter().find(|method| {
        matches!(
            method,
            vakint::EvaluationMethod::RustRed(_) | vakint::EvaluationMethod::Forcer(_)
        )
    });
    match first {
        Some(vakint::EvaluationMethod::RustRed(_)) => Ok(()),
        Some(_) => Err("forcer_before_rustred"),
        None => Err("no_rustred"),
    }
}

/// `VAKINT_P_BOUNDARY_BLOCKS`: `reuse` (default) or `off`.
pub(crate) fn vakint_boundary_blocks() -> &'static str {
    match std::env::var("VAKINT_P_BOUNDARY_BLOCKS").as_deref() {
        Ok("off" | "0") => "off",
        _ => "reuse",
    }
}

/// The non-identity checkpoint field `{"mode": .., "strict": ..}` (DESIGN 0).
pub fn remainder_boundary_provenance() -> serde_json::Value {
    serde_json::json!({
        "mode": match remainder_boundary_mode() {
            Ok(mode) => mode.as_str().to_owned(),
            Err(error) => format!("error({error})"),
        },
        "strict": boundary_strict(),
    })
}

/// The start-up line of a binary that may take the new path (DESIGN 0).
pub fn startup_line() -> String {
    format!(
        "remainder_boundary={} strict={} vakint_gates={} vakint_boundary_blocks={} fallback_max_summands={}",
        match remainder_boundary_mode() {
            Ok(mode) => mode.as_str().to_owned(),
            Err(error) => format!("error({error})"),
        },
        u8::from(boundary_strict()),
        match vakint_gates() {
            Ok(()) => "ok",
            Err(reason) => reason,
        },
        vakint_boundary_blocks(),
        match fallback_max_summands() {
            usize::MAX => "unlimited".to_owned(),
            max => max.to_string(),
        }
    )
}

// ---------------------------------------------------------------------------
// The scalar root normalisation (beta)
// ---------------------------------------------------------------------------

/// The closed metric-root normalisation as data: for every closed,
/// colour-free, index-free scalar `X` in dots, the root normalisation callback
/// gives `factor · X.replace_multiple(&external_routing)` (by value).
#[derive(Clone, Debug)]
pub struct ClosedRootNormalization {
    /// The physical external edge `x0` (beta's `p`).
    pub external: EdgeIndex,
    /// `σ`: `Q(e, args) -> ±Q(x0, args)` for both external edges.
    pub external_routing: Vec<Replacement>,
    /// `c`, with the scalar dimension `dim` not yet substituted.
    pub factor: Atom,
}

// ---------------------------------------------------------------------------
// Declines and counts
// ---------------------------------------------------------------------------

/// A remainder admitted to the new path that the build cannot represent.
#[derive(Clone, Debug)]
pub(super) struct Decline {
    pub stage: &'static str,
    pub reason: String,
}

impl Decline {
    pub(super) fn new(stage: &'static str, reason: impl Into<String>) -> Self {
        Self {
            stage,
            reason: reason.into(),
        }
    }
}

/// The counts of one built remainder (`BOUNDARY_STATS`, DESIGN 8).
#[derive(Clone, Debug, Default)]
pub(super) struct BoundaryStats {
    pub loops: usize,
    pub degree: usize,
    pub p_summands: usize,
    pub p_bytes: usize,
    pub p_truncated: usize,
    pub certified: bool,
    pub notation_retries: usize,
    pub grades: usize,
    pub graded_entries: usize,
    pub d_orders: Vec<usize>,
    pub raw_patterns: usize,
    pub patterns: usize,
    pub zero_patterns: usize,
    pub classes: usize,
    pub blocks: usize,
    pub block_max: (usize, usize),
    pub blocks_total: (usize, usize),
    pub r_terms: usize,
    pub r_summands_total: usize,
    pub r_summands_max: usize,
    pub images: usize,
    pub routing_terms: usize,
    pub block_reads: usize,
    pub block_hits: usize,
    /// Whether the numerator was contracted (the certificate is meaningful).
    pub contracted: bool,
    /// Stage name and seconds, in order.
    pub seconds: Vec<(&'static str, f64)>,
    /// RSS (GiB) after the contraction, the split, the build and the evaluation.
    pub rss_gb: [f64; 4],
}

/// `BOUNDARY_STATS`: one stderr line per built remainder (DESIGN 8).
pub(super) fn log_stats(graph: &str, id: usize, stats: &BoundaryStats) {
    eprintln!(
        "BOUNDARY_STATS graph={graph} id={id} loops={} deg={} P={}/{} P_truncated={} certified={} notation_retries={} grades={} graded_entries={} D_j={} raw_patterns={} patterns={} zero_patterns={} classes={} blocks={} block_max={}/{} blocks_total={}/{} R={}/{}/{} images={} routing_terms={} block_reads={} block_hits={} boundary_blocks={} seconds={} rss_gb={}",
        stats.loops,
        stats.degree,
        stats.p_summands,
        stats.p_bytes,
        stats.p_truncated,
        stats.certified,
        stats.notation_retries,
        stats.grades,
        stats.graded_entries,
        stats
            .d_orders
            .iter()
            .map(|n| n.to_string())
            .collect::<Vec<_>>()
            .join(","),
        stats.raw_patterns,
        stats.patterns,
        stats.zero_patterns,
        stats.classes,
        stats.blocks,
        stats.block_max.0,
        stats.block_max.1,
        stats.blocks_total.0,
        stats.blocks_total.1,
        stats.r_terms,
        stats.r_summands_total,
        stats.r_summands_max,
        stats.images,
        stats.routing_terms,
        stats.block_reads,
        stats.block_hits,
        vakint_boundary_blocks(),
        stats
            .seconds
            .iter()
            .map(|(stage, seconds)| format!("{stage}:{seconds:.3}"))
            .collect::<Vec<_>>()
            .join(","),
        stats
            .rss_gb
            .iter()
            .map(|gb| format!("{gb:.2}"))
            .collect::<Vec<_>>()
            .join(","),
    );
}

/// The remainder-boundary counts of one evaluation (`BOUNDARY_SUMMARY`,
/// DESIGN 8).
#[derive(Debug, Default)]
pub(super) struct BoundarySummary {
    /// Ordinary remainders seen with the mode on or check.
    pub remainders: usize,
    /// Built boundaries (evaluated or checked).
    pub boundary: usize,
    /// Built boundaries that vanish (scaleless or zero).
    pub zero: usize,
    /// Declines after admission (today's path ran).
    pub fallback: usize,
    pub errors: usize,
    pub not_admitted: BTreeMap<String, usize>,
    pub fallback_reasons: BTreeMap<String, usize>,
    /// Contracted numerators with and without idenso's closed-scalar
    /// certificate (DESIGN Addendum A(2)).
    pub certified: usize,
    pub uncertified: usize,
    pub seconds_boundary: f64,
    /// Every today's-path outer integration (ordinary, IR, fallbacks, checks).
    pub seconds_today: f64,
}

impl BoundarySummary {
    pub(super) fn declined(&mut self, decline: &Decline) {
        self.fallback += 1;
        *self
            .fallback_reasons
            .entry(format!("{}:{}", decline.stage, decline.reason))
            .or_default() += 1;
    }

    pub(super) fn contracted(&mut self, stats: &BoundaryStats) {
        if stats.contracted {
            if stats.certified {
                self.certified += 1;
            } else {
                self.uncertified += 1;
            }
        }
    }

    pub(super) fn log(&self, graph: &str, root: &str) {
        let list = |map: &BTreeMap<String, usize>| {
            if map.is_empty() {
                "none".to_owned()
            } else {
                map.iter()
                    .map(|(reason, count)| format!("{reason}:{count}"))
                    .collect::<Vec<_>>()
                    .join(",")
            }
        };
        eprintln!(
            "BOUNDARY_SUMMARY graph={graph} root={root} remainders={} boundary={} zero={} fallback={} errors={} not_admitted={} fallback_reasons={} certified={} uncertified={} seconds_boundary={:.3} seconds_today={:.3}",
            self.remainders,
            self.boundary,
            self.zero,
            self.fallback,
            self.errors,
            list(&self.not_admitted),
            list(&self.fallback_reasons),
            self.certified,
            self.uncertified,
            self.seconds_boundary,
            self.seconds_today,
        );
    }
}

/// The process resident set size in GiB (`/proc/self/statm`).
pub(super) fn rss_gb() -> f64 {
    std::fs::read_to_string("/proc/self/statm")
        .ok()
        .and_then(|statm| statm.split_whitespace().nth(1)?.parse::<u64>().ok())
        .map_or(0.0, |pages| (pages * 4096) as f64 / (1u64 << 30) as f64)
}

/// The stage clock of one build: every stage end is logged as
/// `BOUNDARY_STAGE graph= id= stage= seconds= rss_gb=`.
pub(super) struct StageClock<'a> {
    graph: &'a str,
    id: usize,
    last: Instant,
}

impl<'a> StageClock<'a> {
    pub(super) fn new(graph: &'a str, id: usize) -> Self {
        Self {
            graph,
            id,
            last: Instant::now(),
        }
    }

    /// End `stage`: record and log its seconds and the RSS.
    pub(super) fn lap(&mut self, stats: &mut BoundaryStats, stage: &'static str) -> f64 {
        let now = Instant::now();
        let seconds = now.duration_since(self.last).as_secs_f64();
        self.last = now;
        let rss = rss_gb();
        eprintln!(
            "BOUNDARY_STAGE graph={} id={} stage={stage} seconds={seconds:.3} rss_gb={rss:.2}",
            self.graph, self.id
        );
        stats.seconds.push((stage, seconds));
        seconds
    }
}

// ---------------------------------------------------------------------------
// Grades and the recursive grammar (DESIGN 2a)
// ---------------------------------------------------------------------------

/// `g = (T, s, b)`: the sorted basis positions of the `k·p` factors, the power
/// of `p²` and the power of the reference mass squared.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Hash, Default)]
pub(super) struct Grade {
    pub structure: Vec<usize>,
    pub p_squared: usize,
    pub mass_power: usize,
}

impl Grade {
    /// The λ-degree `a(g) = |T| + 2 s`.
    pub(super) fn degree(&self) -> usize {
        self.structure.len() + 2 * self.p_squared
    }

    fn is_zero(&self) -> bool {
        self.structure.is_empty() && self.p_squared == 0 && self.mass_power == 0
    }

    fn times(&self, other: &Grade) -> Grade {
        let mut structure = self.structure.clone();
        structure.extend_from_slice(&other.structure);
        structure.sort_unstable();
        Grade {
            structure,
            p_squared: self.p_squared + other.p_squared,
            mass_power: self.mass_power + other.mass_power,
        }
    }

    /// `t_{T,s} = Π_{j∈T} dot(k(j+1), p) · dot(p, p)^s` in vakint symbols.
    pub(super) fn structure_atom(&self, p: &Atom) -> Atom {
        let mut factors = self
            .structure
            .iter()
            .map(|j| function!(S.dot, function!(S.k, *j as i64 + 1), p))
            .collect::<Vec<_>>();
        if self.p_squared > 0 {
            factors.push(function!(S.dot, p, p).pow(self.p_squared as i64));
        }
        Atom::mul_many(factors.iter())
    }

}

/// A small map grade -> coefficient (p-free, reference-mass-free, vakint
/// spelling).
pub(super) type GradedMap = BTreeMap<Grade, Atom>;

/// Why a view has no graded map: a decline, or an index slot that the
/// notation retry may still remove.
#[derive(Debug)]
enum Fail {
    Decline(Decline),
    Indexed(String),
}

fn decline<T>(reason: impl Into<String>) -> std::result::Result<T, Fail> {
    Err(Fail::Decline(Decline::new("collect", reason)))
}

/// A vector of a `dot` argument: the basis loop momentum at position `j`, or
/// the external `x0`.
#[derive(Clone, Copy, Debug)]
enum Vector {
    Loop(usize),
    External,
}

fn head_name(symbol: Symbol) -> String {
    symbol.get_name().rsplit("::").next().unwrap_or("").to_owned()
}

/// The λ symbols of the outer and hard Taylor operations.
pub(super) fn lambda_symbols() -> [Symbol; 2] {
    [
        symbol!("rstar_outer_taylor"; Scalar),
        symbol!("rstar_hard_taylor"; Scalar),
    ]
}

/// The symbols that must never reach vakint in a block, coefficient or
/// routing image (DESIGN 2d), with the reference mass squared.
fn forbidden_symbols(reference_mass_squared: Symbol) -> Vec<Symbol> {
    let minkowski = Minkowski {}
        .new_rep(GS.dim)
        .to_symbolic([])
        .get_symbol()
        .expect("the Minkowski representation is a function");
    let ags = &idenso::dirac::AGS;
    let cs = &idenso::color::CS;
    let mut symbols = vec![
        GS.emr_mom,
        GS.loop_mom,
        SPENSO_TAG.dot,
        minkowski,
        ETS.metric,
        ags.gamma,
        ags.gammaconj,
        ags.gammaadj,
        ags.gamma0,
        ags.projp,
        ags.projm,
        ags.gamma5,
        ags.charge_conjugation,
        ags.sigma,
        SPENSO_TAG.trace,
        SPENSO_TAG.chain,
        cs.cas,
        cs.idx,
        cs.gram,
        GS.den,
        GS.dim,
        reference_mass_squared,
    ];
    symbols.extend(lambda_symbols());
    symbols
}

/// The grammar of one remainder: its basis positions, the external `x0`, the
/// truncation degree and the caches of the split.
pub(super) struct Grammar {
    position: HashMap<usize, usize>,
    x0: usize,
    degree: usize,
    reference_mass_squared: Symbol,
    lambdas: [Symbol; 2],
    forbidden: Vec<Symbol>,
    minkowski: Symbol,
    regulated_dimension: Atom,
    p: Atom,
    cache: HashMap<Atom, GradedMap>,
    retry_images: HashMap<Atom, Atom>,
    matcher: SlotMatcher,
}

/// Non-sum factors up to this size are cached by atom (as `map_factors_cached`).
const CACHE_BYTES: usize = 4096;

impl Grammar {
    /// `basis`: the R* basis loop edges `b_j` in `j` order; `x0`: the external
    /// edge; `degree`: the outer Taylor degree.
    pub(super) fn new(
        basis: &[usize],
        x0: usize,
        degree: usize,
        reference_mass_squared: Symbol,
    ) -> Self {
        Self {
            position: basis.iter().enumerate().map(|(j, edge)| (*edge, j)).collect(),
            x0,
            degree,
            reference_mass_squared,
            lambdas: lambda_symbols(),
            forbidden: forbidden_symbols(reference_mass_squared),
            minkowski: Minkowski {}
                .new_rep(GS.dim)
                .to_symbolic([])
                .get_symbol()
                .expect("the Minkowski representation is a function"),
            regulated_dimension: Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon),
            p: function!(S.p, x0 as i64),
            cache: HashMap::new(),
            retry_images: HashMap::new(),
            matcher: SlotMatcher::default(),
        }
    }

    /// The vacuum label `k(j+1)` of basis position `j`.
    fn loop_vector(j: usize) -> Atom {
        function!(S.k, j as i64 + 1)
    }

    /// The external vector `p(x0)` in vakint symbols.
    pub(super) fn external(&self) -> &Atom {
        &self.p
    }

    /// `Q(e)` or `Q(e, mink(d))` (an index-free payload): the vector.
    fn vector(&self, view: AtomView<'_>) -> std::result::Result<Vector, Fail> {
        let AtomView::Fun(function) = view else {
            return decline(format!("dot_argument:{}", kind(view)));
        };
        if function.get_symbol() != GS.emr_mom {
            return decline(format!("dot_argument:{}", head_name(function.get_symbol())));
        }
        let nargs = function.get_nargs();
        if nargs == 0 || nargs > 2 {
            return decline("dot_argument:Q");
        }
        if nargs == 2 {
            match function.get(1) {
                AtomView::Fun(slot)
                    if slot.get_symbol() == self.minkowski && slot.get_nargs() == 1 => {}
                _ => return decline("dot_argument:indexed_Q"),
            }
        }
        let Ok(edge) = usize::try_from(function.get(0)) else {
            return decline("dot_argument:Q");
        };
        if let Some(j) = self.position.get(&edge) {
            Ok(Vector::Loop(*j))
        } else if edge == self.x0 {
            Ok(Vector::External)
        } else {
            decline(format!("dot_argument:Q({edge})"))
        }
    }

    /// A `dot` argument as `Σ α_e v_e` with numbers `α_e`.
    fn linear_vectors(
        &self,
        view: AtomView<'_>,
    ) -> std::result::Result<Vec<(Atom, Vector)>, Fail> {
        let terms = match view {
            AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
            other => vec![other],
        };
        let mut vectors = Vec::with_capacity(terms.len());
        for term in terms {
            let (coefficient, vector) = match term {
                AtomView::Mul(product) => {
                    let factors = product.iter().collect::<Vec<_>>();
                    match factors.as_slice() {
                        [n @ AtomView::Num(_), v] | [v, n @ AtomView::Num(_)] => (n.to_owned(), *v),
                        _ => return decline("dot_argument:nonlinear"),
                    }
                }
                other => (Atom::num(1), other),
            };
            vectors.push((coefficient, self.vector(vector)?));
        }
        Ok(vectors)
    }

    /// The bilinear expansion of a spenso `dot(u, v)` into leaves. Graded:
    /// `k·p` and `p²` become grades; ungraded (`map_scalar`): vakint dots.
    fn dot(
        &self,
        left: AtomView<'_>,
        right: AtomView<'_>,
        graded: bool,
        truncated: &mut usize,
    ) -> std::result::Result<GradedMap, Fail> {
        let left = self.linear_vectors(left)?;
        let right = self.linear_vectors(right)?;
        let mut accumulated = BTreeMap::<Grade, Vec<Atom>>::new();
        for (a, u) in &left {
            for (b, v) in &right {
                let coefficient = a * b;
                let (grade, leaf) = match (u, v) {
                    (Vector::Loop(j), Vector::Loop(l)) => (
                        Grade::default(),
                        function!(S.dot, Self::loop_vector(*j), Self::loop_vector(*l)),
                    ),
                    (Vector::Loop(j), Vector::External) | (Vector::External, Vector::Loop(j)) => {
                        if graded {
                            (
                                Grade {
                                    structure: vec![*j],
                                    ..Default::default()
                                },
                                Atom::num(1),
                            )
                        } else {
                            (
                                Grade::default(),
                                function!(S.dot, Self::loop_vector(*j), &self.p),
                            )
                        }
                    }
                    (Vector::External, Vector::External) => {
                        if graded {
                            (
                                Grade {
                                    p_squared: 1,
                                    ..Default::default()
                                },
                                Atom::num(1),
                            )
                        } else {
                            (Grade::default(), function!(S.dot, &self.p, &self.p))
                        }
                    }
                };
                if grade.degree() > self.degree {
                    *truncated += 1;
                    continue;
                }
                accumulated.entry(grade).or_default().push(coefficient * leaf);
            }
        }
        Ok(collect_map(accumulated))
    }

    /// Whether a function (outside a `dot`) carries an explicit index slot or
    /// an opaque representation, or is an indexed head (metric, gamma, trace,
    /// chain): the notation retry applies.
    fn indexed(&mut self, view: AtomView<'_>) -> Option<String> {
        let indexed_heads = [
            ETS.metric,
            idenso::dirac::AGS.gamma,
            idenso::dirac::AGS.gamma5,
            idenso::dirac::AGS.gamma0,
            SPENSO_TAG.trace,
            SPENSO_TAG.chain,
        ];
        let mut found = None;
        let matcher = &mut self.matcher;
        view.visitor(&mut |part| {
            if found.is_some() {
                return false;
            }
            let AtomView::Fun(function) = part else {
                return true;
            };
            let symbol = function.get_symbol();
            if symbol == SPENSO_TAG.dot {
                return false;
            }
            if indexed_heads.contains(&symbol) {
                found = Some(head_name(symbol));
                return false;
            }
            match matcher.classify(part) {
                SlotMatch::Explicit(_) | SlotMatch::Opaque => {
                    found = Some(head_name(symbol));
                    false
                }
                SlotMatch::Other => true,
            }
        });
        found
    }

    /// The graded map of one view (DESIGN 2a, the recursive grammar).
    fn grade(
        &mut self,
        view: AtomView<'_>,
        stats: &mut BoundaryStats,
    ) -> std::result::Result<GradedMap, Fail> {
        match view {
            AtomView::Num(_) => Ok(single(Grade::default(), view.to_owned())),
            AtomView::Var(variable) => {
                let symbol = variable.get_symbol();
                if symbol == self.reference_mass_squared {
                    Ok(single(
                        Grade {
                            mass_power: 1,
                            ..Default::default()
                        },
                        Atom::num(1),
                    ))
                } else if symbol == GS.dim {
                    Ok(single(Grade::default(), self.regulated_dimension.clone()))
                } else if self.lambdas.contains(&symbol) {
                    decline("lambda")
                } else if self.forbidden.contains(&symbol) {
                    decline(format!("symbol:{}", head_name(symbol)))
                } else {
                    Ok(single(Grade::default(), view.to_owned()))
                }
            }
            AtomView::Add(sum) => {
                let mut accumulated = BTreeMap::<Grade, Vec<Atom>>::new();
                for term in sum.iter() {
                    for (grade, coefficient) in self.grade(term, stats)? {
                        accumulated.entry(grade).or_default().push(coefficient);
                    }
                }
                Ok(collect_map(accumulated))
            }
            AtomView::Mul(product) => {
                let mut factors = Vec::with_capacity(product.get_nargs());
                for factor in product.iter() {
                    factors.push(self.factor(factor, stats)?);
                }
                Ok(self.product(factors, stats))
            }
            AtomView::Pow(_) | AtomView::Fun(_) => self.factor(view, stats),
        }
    }

    /// The graded map of a non-sum factor, cached by atom up to 4 KiB.
    fn factor(
        &mut self,
        view: AtomView<'_>,
        stats: &mut BoundaryStats,
    ) -> std::result::Result<GradedMap, Fail> {
        match view {
            AtomView::Pow(_) | AtomView::Fun(_) if view.get_byte_size() <= CACHE_BYTES => {
                let key = view.to_owned();
                if let Some(map) = self.cache.get(&key) {
                    return Ok(map.clone());
                }
                let map = self.uncached_factor(view, stats)?;
                self.cache.insert(key, map.clone());
                Ok(map)
            }
            AtomView::Pow(_) | AtomView::Fun(_) => self.uncached_factor(view, stats),
            other => self.grade(other, stats),
        }
    }

    fn uncached_factor(
        &mut self,
        view: AtomView<'_>,
        stats: &mut BoundaryStats,
    ) -> std::result::Result<GradedMap, Fail> {
        match view {
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let base_map = self.grade(base, stats)?;
                let integer = i64::try_from(exponent).ok();
                if base_map.is_empty() {
                    // A vanishing (or wholly truncated) base.
                    return match integer {
                        Some(n) if n >= 1 => Ok(GradedMap::new()),
                        _ => decline("graded_power"),
                    };
                }
                if base_map.keys().all(Grade::is_zero) {
                    let exponent = self.scalar(exponent)?;
                    let (grade, coefficient) = base_map.into_iter().next().unwrap();
                    return Ok(single(grade, coefficient.pow(&exponent)));
                }
                match integer {
                    Some(n) if n >= 1 => Ok(self.power(&base_map, n as u64, stats)),
                    Some(_) => decline("graded_power"),
                    None if base_map.keys().any(|grade| grade.mass_power > 0) => {
                        decline("fractional_reference_mass_power")
                    }
                    None => decline("graded_power"),
                }
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                if symbol == SPENSO_TAG.dot {
                    if function.get_nargs() != 2 {
                        return decline("dot_argument:arity");
                    }
                    return self.dot(function.get(0), function.get(1), true, &mut stats.p_truncated);
                }
                if let Some(head) = self.indexed(view) {
                    return Err(Fail::Indexed(head));
                }
                if self.forbidden.contains(&symbol) {
                    return decline(format!("function:{}", head_name(symbol)));
                }
                // No argument may hold a momentum, the reference mass, a λ or
                // any other forbidden head (the scalar `dim` is mapped).
                for argument in function.iter() {
                    if self
                        .forbidden
                        .iter()
                        .any(|f| *f != GS.dim && argument.contains_symbol(*f))
                    {
                        return decline(format!("function:{}", head_name(symbol)));
                    }
                }
                let mut builder = FunctionBuilder::new(symbol);
                for argument in function.iter() {
                    builder = builder.add_arg(self.scalar(argument)?);
                }
                Ok(single(Grade::default(), builder.finish()))
            }
            other => self.grade(other, stats),
        }
    }

    /// A view that must be of grade zero only (exponents, function
    /// arguments): its mapped atom.
    fn scalar(&mut self, view: AtomView<'_>) -> std::result::Result<Atom, Fail> {
        let mut stats = BoundaryStats::default();
        let map = self.grade(view, &mut stats)?;
        if map.keys().any(|grade| !grade.is_zero()) {
            return decline("graded_scalar");
        }
        Ok(map.into_values().next().unwrap_or(Atom::Zero))
    }

    /// The convolution of the factors' maps, truncating `a(g) > degree`.
    fn product(&self, factors: Vec<GradedMap>, stats: &mut BoundaryStats) -> GradedMap {
        // The common case: every factor has one grade.
        if factors.iter().all(|map| map.len() == 1) {
            let mut grade = Grade::default();
            let mut coefficients = Vec::with_capacity(factors.len());
            for map in factors {
                let (g, c) = map.into_iter().next().unwrap();
                grade = grade.times(&g);
                coefficients.push(c);
            }
            if grade.degree() > self.degree {
                stats.p_truncated += 1;
                return GradedMap::new();
            }
            let coefficient = Atom::mul_many(coefficients.iter());
            return if coefficient.is_zero() {
                GradedMap::new()
            } else {
                single(grade, coefficient)
            };
        }
        let mut accumulated = single(Grade::default(), Atom::num(1));
        for (index, factor) in factors.into_iter().enumerate() {
            if index % 4096 == 4095 {
                let _ = symbolica::utils::resources::checkpoint();
            }
            accumulated = self.convolve(&accumulated, &factor, stats);
            if accumulated.is_empty() {
                break;
            }
        }
        accumulated
    }

    fn convolve(&self, a: &GradedMap, b: &GradedMap, stats: &mut BoundaryStats) -> GradedMap {
        let mut accumulated = BTreeMap::<Grade, Vec<Atom>>::new();
        for (ga, ca) in a {
            for (gb, cb) in b {
                let grade = ga.times(gb);
                if grade.degree() > self.degree {
                    stats.p_truncated += 1;
                    continue;
                }
                accumulated.entry(grade).or_default().push(ca * cb);
            }
        }
        collect_map(accumulated)
    }

    /// `base^n` by binary powering, truncated.
    fn power(&self, base: &GradedMap, mut n: u64, stats: &mut BoundaryStats) -> GradedMap {
        let mut result = single(Grade::default(), Atom::num(1));
        let mut square = base.clone();
        while n > 0 {
            if n & 1 == 1 {
                result = self.convolve(&result, &square, stats);
            }
            n >>= 1;
            if n > 0 {
                square = self.convolve(&square, &square, stats);
            }
        }
        result
    }

    /// The graded map of one summand of `P`, with the notation retry: the
    /// summand's indexed factors through today's dot-notation pass (cook,
    /// infer, `to_dots`, uncook; once per distinct product of indexed
    /// factors), then graded again.
    pub(super) fn grade_summand(
        &mut self,
        summand: AtomView<'_>,
        stats: &mut BoundaryStats,
    ) -> std::result::Result<GradedMap, Decline> {
        let head = match self.grade(summand, stats) {
            Ok(map) => return Ok(map),
            Err(Fail::Decline(decline)) => return Err(decline),
            Err(Fail::Indexed(head)) => head,
        };
        stats.notation_retries += 1;
        let notated = match streamed::map_indexed_factors(
            summand,
            &mut self.retry_images,
            &retry_contraction,
        ) {
            Ok(Some(notated)) => notated,
            _ => return Err(Decline::new("collect", format!("indexed:{head}"))),
        };
        match self.grade(notated.as_view(), stats) {
            Ok(map) => Ok(map),
            Err(Fail::Decline(decline)) => Err(decline),
            Err(Fail::Indexed(head)) => Err(Decline::new("collect", format!("indexed:{head}"))),
        }
    }

    /// The leaf rules without grading (DESIGN 2b): `p(x0)` is allowed, dots
    /// of linear combinations are expanded bilinearly, `dim -> 4-2ε`. For the
    /// denominator coefficients and the normalisation.
    pub(super) fn map_scalar(&mut self, view: AtomView<'_>) -> std::result::Result<Atom, Decline> {
        self.map_scalar_inner(view).map_err(|fail| match fail {
            Fail::Decline(decline) => Decline::new("project", decline.reason),
            Fail::Indexed(head) => Decline::new("project", format!("indexed:{head}")),
        })
    }

    fn map_scalar_inner(&mut self, view: AtomView<'_>) -> std::result::Result<Atom, Fail> {
        match view {
            AtomView::Num(_) => Ok(view.to_owned()),
            AtomView::Var(variable) => {
                let symbol = variable.get_symbol();
                if symbol == GS.dim {
                    Ok(self.regulated_dimension.clone())
                } else if self.forbidden.contains(&symbol) {
                    decline(format!("symbol:{}", head_name(symbol)))
                } else {
                    Ok(view.to_owned())
                }
            }
            AtomView::Add(sum) => {
                let mut terms = Vec::with_capacity(sum.get_nargs());
                for term in sum.iter() {
                    terms.push(self.map_scalar_inner(term)?);
                }
                Ok(Atom::add_many(terms.iter()))
            }
            AtomView::Mul(product) => {
                let mut factors = Vec::with_capacity(product.get_nargs());
                for factor in product.iter() {
                    factors.push(self.map_scalar_inner(factor)?);
                }
                Ok(Atom::mul_many(factors.iter()))
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                Ok(self
                    .map_scalar_inner(base)?
                    .pow(&self.map_scalar_inner(exponent)?))
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                if symbol == SPENSO_TAG.dot && function.get_nargs() == 2 {
                    // Ungraded leaves all have grade zero: nothing is truncated.
                    let mut truncated = 0;
                    let map = self.dot(function.get(0), function.get(1), false, &mut truncated)?;
                    return Ok(map.into_values().next().unwrap_or(Atom::Zero));
                }
                if self.forbidden.contains(&symbol)
                    || function.iter().any(|argument| {
                        argument.contains_symbol(GS.emr_mom) || argument.contains_symbol(GS.loop_mom)
                    })
                {
                    return decline(format!("function:{}", head_name(symbol)));
                }
                let mut builder = FunctionBuilder::new(symbol);
                for argument in function.iter() {
                    builder = builder.add_arg(self.map_scalar_inner(argument)?);
                }
                Ok(builder.finish())
            }
        }
    }
}

/// The notation retry's map of a product of indexed factors: the R*
/// remainder contraction itself (cook, infer, `simplify_algebra`, `to_dots`,
/// uncook; [`RStarIntegrand::simplify_view_certified`]). Deviation from
/// DESIGN 2a, which names the dot-notation pass alone (cook, infer,
/// `to_dots`): that pass leaves a closed vector pair `Q(a,μ) Q(b,μ)`
/// unchanged, so the design's own case (c) would decline; both maps keep the
/// value, the contraction also closes such pairs.
pub(super) fn retry_contraction(numerator: AtomView<'_>) -> Result<Atom> {
    Ok(RStarIntegrand::simplify_view_certified(numerator)?.0)
}

fn kind(view: AtomView<'_>) -> &'static str {
    match view {
        AtomView::Num(_) => "number",
        AtomView::Var(_) => "symbol",
        AtomView::Add(_) => "sum",
        AtomView::Mul(_) => "product",
        AtomView::Pow(_) => "power",
        AtomView::Fun(_) => "function",
    }
}

fn single(grade: Grade, coefficient: Atom) -> GradedMap {
    let mut map = GradedMap::new();
    if !coefficient.is_zero() {
        map.insert(grade, coefficient);
    }
    map
}

fn collect_map(accumulated: BTreeMap<Grade, Vec<Atom>>) -> GradedMap {
    accumulated
        .into_iter()
        .filter_map(|(grade, mut coefficients)| {
            let coefficient = if coefficients.len() == 1 {
                coefficients.pop().unwrap()
            } else {
                Atom::add_many(coefficients.iter())
            };
            (!coefficient.is_zero()).then_some((grade, coefficient))
        })
        .collect()
}

/// The graded split of `P` (consumed): every summand's pieces pushed into the
/// block of their grade, then one sum per block. Blocks that vanish are
/// dropped. `Err(decline)` when a summand has no graded map.
pub(super) fn split_into_blocks(
    p: Atom,
    grammar: &mut Grammar,
    stats: &mut BoundaryStats,
) -> Result<std::result::Result<BTreeMap<Grade, Atom>, Decline>> {
    let mut vectors = BTreeMap::<Grade, Vec<Atom>>::new();
    {
        let mut emit = |index: usize, summand: AtomView<'_>| -> Result<Option<Decline>> {
            if index % 4096 == 4095 {
                symbolica::utils::resources::checkpoint()
                    .map_err(|error| eyre!("{error}"))
                    .wrap_err("while splitting an R* remainder into boundary blocks")?;
            }
            let map = match grammar.grade_summand(summand, stats) {
                Ok(map) => map,
                Err(decline) => return Ok(Some(decline)),
            };
            for (grade, coefficient) in map {
                stats.graded_entries += 1;
                vectors.entry(grade).or_default().push(coefficient);
            }
            Ok(None)
        };
        match p.as_view() {
            AtomView::Add(sum) => {
                for (index, summand) in sum.iter().enumerate() {
                    if let Some(decline) = emit(index, summand)? {
                        return Ok(Err(decline));
                    }
                }
            }
            single => {
                if let Some(decline) = emit(0, single)? {
                    return Ok(Err(decline));
                }
            }
        }
    }
    drop(p);
    let mut blocks = BTreeMap::new();
    while let Some((grade, coefficients)) = vectors.pop_first() {
        let block = Atom::add_many(coefficients.iter());
        drop(coefficients);
        if !block.is_zero() {
            blocks.insert(grade, block);
        }
    }
    stats.grades = blocks.len();
    Ok(Ok(blocks))
}

// ---------------------------------------------------------------------------
// Postconditions (DESIGN 2d)
// ---------------------------------------------------------------------------

/// Where a scalar goes in the boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Role {
    /// A block `S_g`.
    Block,
    /// A coefficient `R_{π,g}`.
    Coefficient,
    /// A routing image.
    Routing,
}

impl Role {
    fn as_str(&self) -> &'static str {
        match self {
            Role::Block => "block",
            Role::Coefficient => "coefficient",
            Role::Routing => "routing",
        }
    }
}

/// The vakint contract of one block, coefficient or routing image: nothing
/// forbidden anywhere; `k(n)` only as a direct `dot` argument with
/// `1 <= n <= loops` (blocks, coefficients) or as a linear term with
/// `1 <= n < loops` (routing); `p(n)`: none in blocks, only `p(x0)` inside a
/// `dot` in coefficients, only the closing head `p(head)` in routing images.
pub(super) fn verify_vakint_scalar(
    view: AtomView<'_>,
    role: Role,
    loops: usize,
    x0: usize,
    head: usize,
    reference_mass_squared: Symbol,
) -> std::result::Result<(), Decline> {
    let forbidden = forbidden_symbols(reference_mass_squared);
    let fail = |what: String| Decline::new("verify", format!("postcondition:{}:{what}", role.as_str()));
    let vector_index = |view: AtomView<'_>, symbol: Symbol| -> Option<i64> {
        match view {
            AtomView::Fun(function) if function.get_symbol() == symbol && function.get_nargs() == 1 => {
                i64::try_from(function.get(0)).ok()
            }
            _ => None,
        }
    };
    if role == Role::Routing {
        // A linear combination of cut loop momenta and the closing head.
        let terms = match view {
            AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
            other => vec![other],
        };
        for term in terms {
            let vector = match term {
                AtomView::Mul(product) => {
                    let factors = product.iter().collect::<Vec<_>>();
                    match factors.as_slice() {
                        [AtomView::Num(_), v] | [v, AtomView::Num(_)] => *v,
                        _ => return Err(fail("nonlinear".into())),
                    }
                }
                other => other,
            };
            if let Some(n) = vector_index(vector, S.k) {
                if n < 1 || n as usize >= loops {
                    return Err(fail(format!("k({n})")));
                }
            } else if let Some(n) = vector_index(vector, S.p) {
                if n as usize != head {
                    return Err(fail(format!("p({n})")));
                }
            } else {
                return Err(fail(kind(vector).into()));
            }
        }
        return Ok(());
    }
    fn walk(
        view: AtomView<'_>,
        inside_dot: bool,
        check: &mut dyn FnMut(AtomView<'_>, bool) -> std::result::Result<bool, Decline>,
    ) -> std::result::Result<(), Decline> {
        if !check(view, inside_dot)? {
            return Ok(());
        }
        match view {
            AtomView::Num(_) | AtomView::Var(_) => Ok(()),
            AtomView::Add(sum) => sum.iter().try_for_each(|v| walk(v, false, check)),
            AtomView::Mul(product) => product.iter().try_for_each(|v| walk(v, false, check)),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                walk(base, false, check)?;
                walk(exponent, false, check)
            }
            AtomView::Fun(function) => {
                let dot = function.get_symbol() == S.dot;
                function.iter().try_for_each(|v| walk(v, dot, check))
            }
        }
    }
    let mut check = |view: AtomView<'_>, inside_dot: bool| -> std::result::Result<bool, Decline> {
        match view {
            AtomView::Var(variable) => {
                let symbol = variable.get_symbol();
                if forbidden.contains(&symbol) {
                    return Err(fail(head_name(symbol)));
                }
                Ok(false)
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                if forbidden.contains(&symbol) {
                    return Err(fail(head_name(symbol)));
                }
                if symbol == S.dot && function.get_nargs() != 2 {
                    return Err(fail("dot_arity".into()));
                }
                if symbol == S.k {
                    let Some(n) = vector_index(view, S.k) else {
                        return Err(fail("k".into()));
                    };
                    if !inside_dot || n < 1 || n as usize > loops {
                        return Err(fail(format!("k({n})")));
                    }
                    return Ok(false);
                }
                if symbol == S.p {
                    let Some(n) = vector_index(view, S.p) else {
                        return Err(fail("p".into()));
                    };
                    if role == Role::Block || !inside_dot || n as usize != x0 {
                        return Err(fail(format!("p({n})")));
                    }
                    return Ok(false);
                }
                Ok(true)
            }
            _ => Ok(true),
        }
    };
    walk(view, false, &mut check)
}

// ---------------------------------------------------------------------------
// Denominator orders (DESIGN 2b)
// ---------------------------------------------------------------------------

/// Whether a factor is a `den` or a power of one.
fn is_denominator(factor: AtomView<'_>) -> bool {
    match factor {
        AtomView::Fun(function) => function.get_symbol() == GS.den,
        AtomView::Pow(power) => matches!(
            power.get_base_exp().0,
            AtomView::Fun(function) if function.get_symbol() == GS.den
        ),
        _ => false,
    }
}

/// The Taylor orders `D_0 .. D_deg` of the scaled and compensated
/// denominator product, each split into its denominator monomials `μ`:
/// `μ -> [d_{0,μ} .. d_{deg,μ}]` (coefficients unmapped), and the summand
/// count of every order.
#[allow(clippy::too_many_arguments)]
pub(super) fn denominator_orders(
    denominators: &[Atom],
    root_edges: &[EdgeIndex],
    degree: usize,
    reference_mass_squared: Symbol,
    auxiliary_mass_squared: &Atom,
) -> Result<std::result::Result<(BTreeMap<Atom, Vec<Atom>>, Vec<usize>), Decline>> {
    let parameter = symbol!("rstar_outer_taylor"; Scalar);
    let taylor = TaylorExpansion::new(parameter, degree);
    let product = Atom::mul_many(denominators.iter());
    let scaled = RStarIntegrand::scale_edges(&product, root_edges, parameter);
    let compensated = taylor.compensate_mass(&scaled, reference_mass_squared, auxiliary_mass_squared);
    let orders = match taylor.orders(&compensated) {
        Ok(orders) => orders,
        Err(error) => return Ok(Err(Decline::new("den", format!("series: {error}")))),
    };
    let mut by_monomial = BTreeMap::<Atom, Vec<Vec<Atom>>>::new();
    let mut sizes = Vec::with_capacity(orders.len());
    for (j, order) in orders.iter().enumerate() {
        let expanded = order.expand();
        let summands = match expanded.as_view() {
            AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
            AtomView::Num(_) if expanded.is_zero() => Vec::new(),
            single => vec![single],
        };
        sizes.push(summands.len());
        for summand in summands {
            let factors = match summand {
                AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
                single => vec![single],
            };
            let (dens, others): (Vec<_>, Vec<_>) =
                factors.into_iter().partition(|factor| is_denominator(*factor));
            if dens.is_empty() {
                return Ok(Err(Decline::new("den", "a denominator order term without denominators")));
            }
            if others.iter().any(|factor| factor.contains_symbol(GS.den)) {
                return Ok(Err(Decline::new("den", "a denominator inside a coefficient")));
            }
            let monomial = Atom::mul_many(dens.iter().copied());
            let coefficient = Atom::mul_many(others.iter().copied());
            by_monomial
                .entry(monomial)
                .or_insert_with(|| vec![Vec::new(); degree + 1])[j]
                .push(coefficient);
        }
    }
    let by_monomial = by_monomial
        .into_iter()
        .map(|(monomial, orders)| {
            (
                monomial,
                orders
                    .into_iter()
                    .map(|coefficients| Atom::add_many(coefficients.iter()))
                    .collect(),
            )
        })
        .collect();
    Ok(Ok((by_monomial, sizes)))
}

// ---------------------------------------------------------------------------
// The build
// ---------------------------------------------------------------------------

/// Everything the build needs from the evaluation, besides the remainder.
pub(super) struct BoundaryRequest<'b> {
    pub graph_name: &'b str,
    pub marked_edge: EdgeIndex,
    /// The integration id this remainder receives (`integration_id + 1`).
    pub integration_id: usize,
    pub root_edges: &'b [EdgeIndex],
    pub degree: usize,
    pub reference_mass_squared: Symbol,
    pub auxiliary_mass_squared: &'b Atom,
    pub vakint: &'b Vakint,
    /// `VakintSettings` with the remainder's Laurent depth.
    pub settings: &'b VakintSettings,
    /// The vakint `additional_normalization` (one power per vacuum loop).
    pub additional_normalization: &'b str,
    pub epsilon: Symbol,
    pub closed_root: &'b ClosedRootNormalization,
}

/// One kept or zero pattern, for check mode and records (the key, the zero
/// patterns and the grades are for the census record of DESIGN 6b).
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub(super) struct PatternMeta {
    /// The grouping key (today's topology key, `topo(prop(..))`).
    pub key: Atom,
    /// `I_π`: today's `VakintTerm.integral` (rebuilt, routed, vakint symbols).
    pub integral: Atom,
    pub closure_power: i64,
    pub unprojected_nonzero: bool,
}

/// A built boundary with its metadata.
#[allow(dead_code)]
pub(super) struct BuiltBoundary {
    pub boundary: RemainderBoundary,
    /// Per kept pattern, in the order of `boundary.patterns`.
    pub patterns: Vec<PatternMeta>,
    pub zero_patterns: Vec<PatternMeta>,
    /// Per block index, its grade.
    pub grades: Vec<Grade>,
    /// The vacuum loop count `L`.
    pub loops: usize,
}

/// The outcome of a build.
pub(super) enum Built {
    Boundary(BuiltBoundary, BoundaryStats),
    /// Every pattern is scaleless, or the remainder vanishes: the value is 0.
    Zero(BoundaryStats),
    Declined(Decline, BoundaryStats),
}

/// One pattern under construction.
struct Pattern {
    key: Atom,
    /// `(μ, rest)` of the raw patterns merged into this one.
    raw: Vec<(Atom, Atom)>,
    frame: TermFrame,
    cut: VacuumCut,
}

/// Build the boundary of one admitted remainder, cheap stages first (DESIGN
/// 6a): den, frames, cut, routing, contract, collect, project, assemble,
/// verify. Hard errors (a failing contraction, a failing projection) are
/// errors; anything the boundary cannot represent is a decline.
pub(super) fn build_boundary(remainder: &RStarIntegrand, request: &BoundaryRequest<'_>) -> Result<Built> {
    let loops = remainder.region.loops;
    let mut stats = BoundaryStats {
        loops,
        degree: request.degree,
        ..Default::default()
    };
    let mut clock = StageClock::new(request.graph_name, request.integration_id);
    macro_rules! declined {
        ($decline:expr) => {
            return Ok(Built::Declined($decline, stats))
        };
    }

    // den: the denominator product alone.
    if remainder.routed_is_zero() {
        clock.lap(&mut stats, "den");
        return Ok(Built::Zero(stats));
    }
    let (numerator, denominators) = match remainder.split_routed() {
        Ok(split) => split,
        Err(reason) => declined!(Decline::new("den", reason)),
    };
    let (by_monomial, sizes) = match denominator_orders(
        &denominators,
        request.root_edges,
        request.degree,
        request.reference_mass_squared,
        request.auxiliary_mass_squared,
    )? {
        Ok(orders) => orders,
        Err(decline) => declined!(decline),
    };
    stats.d_orders = sizes;
    stats.raw_patterns = by_monomial.len();
    clock.lap(&mut stats, "den");

    // frames: today's topology key per monomial, merged; the vacuum frame.
    let converter = remainder.converter_graph();
    let rewrites = TermRewrites::new(
        &converter.underlying,
        &remainder.region.remaining,
        &converter.empty_subgraph::<SuBitGraph>(),
        false,
    );
    let mut grouped = BTreeMap::<Atom, Vec<(Atom, Atom)>>::new();
    for monomial in by_monomial.keys() {
        let Some((key, rest)) = rewrites.denominator_image(monomial) else {
            declined!(Decline::new("frames", "no single topology"));
        };
        grouped.entry(key).or_default().push((monomial.clone(), rest));
    }
    let mut frames = Vec::with_capacity(grouped.len());
    for (key, raw) in grouped {
        let frame = match term_frame(&key) {
            Ok(frame) => frame,
            Err(error) => declined!(Decline::new("frames", format!("momentum solve: {error:#}"))),
        };
        frames.push((key, raw, frame));
    }
    clock.lap(&mut stats, "frames");

    // cut: once per pattern, at its massive line.
    let x0 = usize::from(request.closed_root.external);
    let probe = function!(S.p, x0 as i64);
    let mut patterns = Vec::with_capacity(frames.len());
    for (key, raw, frame) in frames {
        let cut = match cut_frame(&frame.integral, &probe) {
            Ok(Some(cut)) => cut,
            Ok(None) => declined!(Decline::new("cut", "no massive line")),
            Err(error) => declined!(Decline::new("cut", format!("{error:#}"))),
        };
        patterns.push(Pattern {
            key,
            raw,
            frame,
            cut,
        });
    }
    let scaleless = patterns.iter().filter(|pattern| pattern.cut.scaleless).count();
    if scaleless == patterns.len() {
        clock.lap(&mut stats, "cut");
        return Ok(Built::Zero(stats));
    }
    if scaleless != 0 {
        declined!(Decline::new("cut", "scaleless patterns beside others"));
    }
    if patterns
        .iter()
        .any(|pattern| !matches!(pattern.cut.integrand.as_view(), AtomView::Fun(_)))
    {
        declined!(Decline::new("cut", "a non-topo cut"));
    }
    clock.lap(&mut stats, "cut");

    // routing: the composite vacuum-label routing, the same for every pattern.
    let basis = remainder
        .loop_basis_edges()
        .iter()
        .map(|edge| usize::from(*edge))
        .collect::<Vec<_>>();
    if basis.len() != loops {
        declined!(Decline::new("routing", "the loop basis does not match the loop count"));
    }
    let mut routing: Option<Vec<(Atom, Atom)>> = None;
    for pattern in &patterns {
        let mut images = Vec::with_capacity(loops);
        for (j, edge) in basis.iter().enumerate() {
            let Some(image) = pattern.frame.loop_image(*edge) else {
                declined!(Decline::new("routing", format!("no loop image of Q({edge})")));
            };
            images.push((
                Grammar::loop_vector(j),
                image.replace_multiple(&pattern.cut.routing).expand(),
            ));
        }
        match &routing {
            None => routing = Some(images),
            Some(first) if *first == images => {}
            Some(_) => declined!(Decline::new("routing", "routing_differs")),
        }
    }
    let routing = routing.expect("at least one pattern");
    let first = &patterns[0].cut.closure;
    if patterns.iter().any(|pattern| {
        pattern.cut.closure.momentum != first.momentum
            || pattern.cut.closure.mass_squared != first.mass_squared
    }) {
        declined!(Decline::new("routing", "the closing momentum or mass differs"));
    }
    let closing = MassiveClosure {
        momentum: first.momentum.clone(),
        mass_squared: first.mass_squared.clone(),
        power: first.power,
    };
    let head = match closing.momentum.as_view() {
        AtomView::Fun(function) if function.get_symbol() == S.p && function.get_nargs() == 1 => {
            usize::try_from(function.get(0)).ok()
        }
        _ => None,
    };
    let Some(head) = head else {
        declined!(Decline::new("routing", "the closing momentum is not p(h)"));
    };
    stats.routing_terms = routing
        .iter()
        .map(|(_, image)| match image.as_view() {
            AtomView::Add(sum) => sum.get_nargs(),
            _ => 1,
        })
        .sum();
    clock.lap(&mut stats, "routing");

    // contract: P = C(σ(N)), without the denominators.
    let (p, certified) = remainder
        .numerator_contracted(numerator, &request.closed_root.external_routing)
        .wrap_err("while contracting the numerator of an R* remainder boundary")?;
    stats.certified = certified;
    stats.contracted = true;
    stats.p_summands = vakint::census::summands(p.as_view());
    stats.p_bytes = p.as_view().get_byte_size();
    remainder.assert_no_combinations()?;
    clock.lap(&mut stats, "contract");
    stats.rss_gb[0] = rss_gb();

    // collect: the graded split of P into blocks.
    let mut grammar = Grammar::new(&basis, x0, request.degree, request.reference_mass_squared);
    let blocks = match split_into_blocks(p, &mut grammar, &mut stats)? {
        Ok(blocks) => blocks,
        Err(decline) => declined!(decline),
    };
    clock.lap(&mut stats, "collect");
    stats.rss_gb[1] = rss_gb();
    if blocks.is_empty() {
        return Ok(Built::Zero(stats));
    }

    // project: R_{π,g} = c̃ · P[ρ_{π,g}], expanded.
    let normalization = match grammar.map_scalar(
        request
            .closed_root
            .factor
            .replace_multiple(&request.closed_root.external_routing)
            .as_view(),
    ) {
        Ok(factor) => factor,
        Err(decline) => declined!(decline),
    };
    let additional = symbolica::parse!(request.additional_normalization);
    let mut tensor_settings = request.settings.clone();
    tensor_settings.tensor_reduction_method = TensorReductionMethod::FeynKit;
    tensor_settings.use_dot_product_notation = true;
    let mut images = HashMap::<(Atom, Atom), Atom>::new();
    let auxiliary = request.auxiliary_mass_squared;
    let mut mapped_orders = HashMap::<Atom, Vec<Atom>>::new();
    for (monomial, orders) in &by_monomial {
        let mut mapped = Vec::with_capacity(orders.len());
        for order in orders {
            let routed = order.replace_multiple(&request.closed_root.external_routing);
            match grammar.map_scalar(routed.as_view()) {
                Ok(order) => mapped.push(order),
                Err(decline) => declined!(decline),
            }
        }
        mapped_orders.insert(monomial.clone(), mapped);
    }
    let block_grades = blocks.keys().cloned().collect::<Vec<_>>();
    let mut kept = Vec::new();
    let mut zero_patterns = Vec::new();
    for pattern in &patterns {
        // d_{j,π} = Σ_{μ∈π} rest_μ · d_{j,μ}.
        let mut orders = vec![Vec::<Atom>::new(); request.degree + 1];
        for (monomial, rest) in &pattern.raw {
            let rest = match grammar.map_scalar(rest.as_view()) {
                Ok(rest) => rest,
                Err(decline) => declined!(decline),
            };
            for (j, order) in mapped_orders[monomial].iter().enumerate() {
                if !order.is_zero() {
                    orders[j].push(&rest * order);
                }
            }
        }
        let orders = orders
            .iter()
            .map(|terms| Atom::add_many(terms.iter()))
            .collect::<Vec<_>>();
        let normalization = &normalization * additional.pow(pattern.frame.nloops as i64);
        let mut terms = Vec::new();
        let mut unprojected_nonzero = false;
        for (index, grade) in block_grades.iter().enumerate() {
            // ρ = t_g (M²)^b Σ_c C(b,c) (-1)^c Σ_{j ≤ deg-a-2c, j ≡ a (2)} d_j.
            let a = grade.degree();
            let b = grade.mass_power;
            let mut sum = Vec::new();
            for c in 0..=b {
                let Some(limit) = request.degree.checked_sub(a + 2 * c) else {
                    break;
                };
                let binomial = binomial(b, c) * if c % 2 == 0 { 1 } else { -1 };
                for (j, order) in orders.iter().enumerate().take(limit + 1) {
                    if j % 2 == a % 2 && !order.is_zero() {
                        sum.push(Atom::num(binomial) * order);
                    }
                }
            }
            let sum = Atom::add_many(sum.iter());
            if sum.is_zero() {
                continue;
            }
            unprojected_nonzero = true;
            let rho = (grade.structure_atom(grammar.external())
                * auxiliary.pow(b as i64)
                * sum)
                .expand();
            if rho.is_zero() {
                continue;
            }
            let vectors = vakint::VakintTerm::identify_vectors_in_numerator(rho.as_view())
                .map_err(|error| eyre!("{error}"))?;
            let term = vakint::VakintTerm {
                integral: patterns[0].frame.integral.clone(),
                numerator: rho,
                vectors,
            };
            let projected = tensor_reduce_by_factors(request.vakint, &tensor_settings, &term, &mut images)?
                .unwrap_or_else(|| term.numerator.clone());
            let coefficient = (&normalization * projected).expand();
            if coefficient.is_zero() {
                continue;
            }
            let summands = vakint::census::summands(coefficient.as_view());
            stats.r_terms += 1;
            stats.r_summands_total += summands;
            stats.r_summands_max = stats.r_summands_max.max(summands);
            terms.push((index, coefficient));
        }
        let meta = PatternMeta {
            key: pattern.key.clone(),
            integral: pattern.frame.integral.clone(),
            closure_power: pattern.cut.closure.power,
            unprojected_nonzero,
        };
        if terms.is_empty() {
            zero_patterns.push(meta);
        } else {
            kept.push((
                meta,
                BoundaryPattern {
                    topology: pattern.cut.integrand.clone(),
                    closure_power: pattern.cut.closure.power,
                    terms,
                },
            ));
        }
        let _ = symbolica::utils::resources::checkpoint();
    }
    stats.images = images.len();
    stats.patterns = kept.len();
    stats.zero_patterns = zero_patterns.len();
    stats.classes = kept
        .iter()
        .map(|(meta, _)| meta.closure_power)
        .collect::<BTreeSet<_>>()
        .len();
    clock.lap(&mut stats, "project");
    if kept.is_empty() {
        return Ok(Built::Zero(stats));
    }

    // assemble: the blocks in grade order (every block is used by some
    // pattern or kept for the record), the patterns in key order.
    let blocks = blocks.into_values().collect::<Vec<_>>();
    stats.blocks = blocks.len();
    for block in &blocks {
        let size = (
            vakint::census::summands(block.as_view()),
            block.as_view().get_byte_size(),
        );
        stats.block_max = (stats.block_max.0.max(size.0), stats.block_max.1.max(size.1));
        stats.blocks_total = (stats.blocks_total.0 + size.0, stats.blocks_total.1 + size.1);
    }
    let label = format!(
        "graph={} marked_edge={} integration_id={} loops={loops}",
        request.graph_name,
        usize::from(request.marked_edge),
        request.integration_id
    );
    let (metas, boundary_patterns): (Vec<_>, Vec<_>) = kept.into_iter().unzip();
    let boundary = match RemainderBoundary::new(
        label,
        closing.momentum.clone(),
        closing.mass_squared.clone(),
        routing,
        blocks,
        boundary_patterns,
    ) {
        Ok(boundary) => boundary,
        Err(error) => declined!(Decline::new("assemble", format!("{error}"))),
    };
    clock.lap(&mut stats, "assemble");

    // verify: the vakint contract of every block, coefficient and image.
    for block in &boundary.blocks {
        if let Err(decline) = verify_vakint_scalar(
            block.as_view(),
            Role::Block,
            loops,
            x0,
            head,
            request.reference_mass_squared,
        ) {
            declined!(decline);
        }
    }
    for pattern in &boundary.patterns {
        for (_, coefficient) in &pattern.terms {
            if let Err(decline) = verify_vakint_scalar(
                coefficient.as_view(),
                Role::Coefficient,
                loops,
                x0,
                head,
                request.reference_mass_squared,
            ) {
                declined!(decline);
            }
        }
        // The vakint contract of the cut and the closing line.
        let mut heads = BTreeSet::new();
        pattern.topology.visitor(&mut |part| {
            if let AtomView::Fun(function) = part
                && function.get_symbol() == S.p
            {
                heads.insert(part.to_owned());
                return false;
            }
            true
        });
        if heads.len() != 1 || !heads.contains(&closing.momentum) {
            declined!(Decline::new("verify", "vakint_contract:cut_head"));
        }
        if pattern.closure_power < 1 {
            declined!(Decline::new("verify", "vakint_contract:closure_power"));
        }
    }
    if closing.mass_squared.is_zero()
        || closing.mass_squared.contains_symbol(S.k)
        || closing.mass_squared.contains_symbol(S.p)
        || closing.mass_squared.contains_symbol(request.epsilon)
    {
        declined!(Decline::new("verify", "vakint_contract:closing_mass"));
    }
    for (_, image) in &boundary.routing {
        if let Err(decline) = verify_vakint_scalar(
            image.as_view(),
            Role::Routing,
            loops,
            x0,
            head,
            request.reference_mass_squared,
        ) {
            declined!(decline);
        }
    }
    clock.lap(&mut stats, "verify");
    stats.rss_gb[2] = rss_gb();
    Ok(Built::Boundary(
        BuiltBoundary {
            boundary,
            patterns: metas,
            zero_patterns,
            grades: block_grades,
            loops,
        },
        stats,
    ))
}

fn binomial(n: usize, k: usize) -> i64 {
    (0..k).fold(1i64, |value, i| value * (n - i) as i64 / (i as i64 + 1))
}

/// [`RemainderBoundary::evaluate`] with the block read counts: the value
/// (expanded), and `(reads, hits)`.
pub(super) fn evaluate_boundary(
    boundary: &RemainderBoundary,
    vakint: &Vakint,
    settings: &VakintSettings,
    options: BoundaryOptions,
) -> Result<(Atom, usize, usize)> {
    // As `RemainderBoundary::evaluate` (the settings' RustRed resource
    // budget around the whole lowering and reduction), keeping the counts.
    symbolica::utils::resources::ResourceBudget::new(settings.rustred_resources).run(|| {
        symbolica::utils::resources::checkpoint().map_err(|error| eyre!("{error}"))?;
        let (classes, declined, (reads, hits)) = boundary
            .lower_counted(vakint, settings, options)
            .map_err(|error| eyre!("{error}"))
            .wrap_err("while lowering an R* remainder boundary")?;
        let mut result = Atom::Zero;
        for class in classes {
            result += class
                .evaluate(vakint, settings)
                .map_err(|error| eyre!("{error}"))
                .wrap_err("while reducing an R* remainder boundary")?;
        }
        for index in declined {
            let pattern = &boundary.patterns[index];
            result += MassiveClosure::evaluate_batch(
                vakint,
                settings,
                &[(boundary.closure(pattern), boundary.direct_input(pattern).ordinary())],
            )
            .map_err(|error| eyre!("{error}"))?
            .pop()
            .unwrap_or_default();
        }
        symbolica::utils::resources::checkpoint().map_err(|error| eyre!("{error}"))?;
        Ok((result.expand(), reads, hits))
    })
}
