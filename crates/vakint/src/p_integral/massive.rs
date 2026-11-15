//! Close a massless propagator integral with one massive propagator.
//!
//! Keep its exact homogeneity through the final radial integration: expanding
//! the p-integral first would turn the necessary epsilon-dependent power of
//! the hard momentum into logarithms and lose the simple Gamma integral.

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::function;
use symbolica::transcendental::gamma;

use super::{ReducedPropagatorValue, external_head, forcer, laurent, scalar_numerator};
use crate::symbols::S;
use crate::utils::simplify_real;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{EvaluationMethod, RustRedMasterMode, Vakint, VakintError, VakintExpression};

/// The massive line used to close an already cut, scalar vacuum graph.
///
/// `momentum` is the massless graph's external `p(i)`. `evaluate` integrates
/// that momentum with an extra `(p² - mass_squared)^(-power)` and the same
/// per-loop measure as its p-integral. All other propagators must be massless.
/// The caller owns the graph cut and its momentum routing.
#[derive(Clone, Debug)]
pub struct MassiveClosure {
    pub momentum: Atom,
    pub mass_squared: Atom,
    pub power: i64,
}

/// A closure input with its cut graph, cut routing and vacuum numerator
/// apart, so that the numerator can be lowered onto the RustRed p library
/// without substituting the routing into it and expanding it
/// ([`MassiveClosure::evaluate_direct_batch`]).
#[derive(Clone, Debug)]
pub struct DirectClosureInput {
    /// `topo(prop(...))` of the cut graph, in the cut loop momenta `k(1..)`
    /// and the closing head.
    pub topology: Atom,
    /// Every vacuum loop momentum `k(i)` and its image in the cut momenta.
    pub routing: Vec<(Atom, Atom)>,
    /// The scalar vacuum numerator (after any tensor projection).
    pub numerator: Atom,
}

impl DirectClosureInput {
    /// The ordinary closure input `numerator(routing).expand() * topology`.
    pub fn ordinary(&self) -> Atom {
        let routing = self
            .routing
            .iter()
            .map(|(head, image)| symbolica::id::Replacement::new(head.to_pattern(), image.clone()))
            .collect::<Vec<_>>();
        self.numerator.replace_multiple(&routing).expand() * &self.topology
    }
}

/// The folded RustRed input of one or more closures with the same closing
/// line (mass, power, head): an additive map over canonical integral keys.
#[derive(Clone, Debug)]
pub struct LoweredClosure {
    closure: MassiveClosure,
    head: Atom,
    folded: std::collections::BTreeMap<
        crate::rustred_evaluation::p_integral::Channel,
        crate::rustred_evaluation::p_integral::FoldedInputs,
    >,
}

impl LoweredClosure {
    /// The closing line this input belongs to: (mass², power, head).
    pub fn closing_key(&self) -> (Atom, i64, Atom) {
        (self.closure.mass_squared.clone(), self.closure.power, self.head.clone())
    }

    /// Add `other` (same closing key) into this input; `Err(other)` otherwise.
    /// Values are exact rationals, summed per channel, ε-basis and key.
    pub fn absorb(&mut self, other: LoweredClosure) -> Result<(), LoweredClosure> {
        if other.closing_key() != self.closing_key() {
            return Err(other);
        }
        for (channel, inputs) in other.folded {
            let target = self.folded.entry(channel).or_default();
            target.integrals.extend(inputs.integrals);
            for (basis, values) in inputs.by_basis {
                let bucket = target.by_basis.entry(basis).or_default();
                for (key, value) in values {
                    let entry = bucket.entry(key).or_insert_with(symbolica::domains::rational::Rational::zero);
                    *entry += &value;
                }
                bucket.retain(|_, value| !value.is_zero());
            }
        }
        Ok(())
    }

    /// Reduce and finalise once (RustRed master mode of the settings).
    pub fn evaluate(self, vakint: &Vakint, settings: &crate::VakintSettings) -> Result<Atom, VakintError> {
        use crate::rustred_evaluation::p_integral as rr;
        let handle = vakint.p_library()?;
        let mode = settings
            .evaluation_order
            .0
            .iter()
            .find_map(|method| match method {
                EvaluationMethod::RustRed(options) => Some(options.master_mode),
                _ => None,
            })
            .unwrap_or(RustRedMasterMode::Symbolic);
        let reduced = rr::reduce_folded_channels(vakint, settings, &handle, self.folded, None, rr::fold_mode(), true)?;
        let mut result = Atom::Zero;
        let mut numerical = Atom::Zero;
        for ((loops, _), values) in reduced {
            let value = self.closure.finalize(&values, loops, &self.head, settings, mode)?;
            if mode == RustRedMasterMode::Numerical {
                numerical += value;
            } else {
                result += value;
            }
        }
        if !numerical.is_zero() {
            result += vakint.resolve_rustred_numerical(numerical, settings)?;
        }
        Ok(result.expand())
    }
}

/// `VAKINT_LOWERING`: `legacy` evaluates direct closure inputs through the
/// ordinary path; `check` lowers them directly and also through the ordinary
/// path, compares the values (`LOWERING_CHECK` on stderr) and returns the
/// direct value.
fn lowering_mode() -> &'static str {
    static MODE: std::sync::LazyLock<String> = std::sync::LazyLock::new(|| {
        std::env::var("VAKINT_LOWERING").unwrap_or_else(|_| "direct".to_owned())
    });
    MODE.as_str()
}

/// Highest power of loop momenta in a scalar p-integral numerator: `dot(k,k)`
/// counts two, `dot(k,p)` one (an upper bound for nested sums).
pub(crate) fn loop_rank(view: AtomView) -> usize {
    match view {
        AtomView::Add(sum) => sum.iter().map(loop_rank).max().unwrap_or(0),
        AtomView::Mul(product) => product.iter().map(loop_rank).sum(),
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            let exponent = crate::get_integer_from_atom(exponent)
                .filter(|n| *n > 0)
                .unwrap_or(0);
            exponent as usize * loop_rank(base)
        }
        AtomView::Fun(function) if function.get_symbol() == S.dot => function
            .iter()
            .filter(|argument| argument.contains_symbol(S.k))
            .count(),
        _ => 0,
    }
}

/// Experimental hybrid (`BETA_RUSTRED_MAX_RANK`): with RustRed selected, a
/// massive-closure p-integral of higher numerator rank goes to Forcer, where
/// RustRed's reduction is slower.
pub(crate) fn rustred_max_rank() -> Option<usize> {
    static RANK: std::sync::LazyLock<Option<usize>> = std::sync::LazyLock::new(|| {
        std::env::var("BETA_RUSTRED_MAX_RANK")
            .ok()
            .and_then(|value| value.parse().ok())
    });
    *RANK
}

/// The same threshold for the physical (non-closure) p-integrals,
/// `BETA_RUSTRED_MAX_RANK_PHYSICAL`; unset keeps RustRed for every rank there,
/// which measured faster on the four-loop sample.
pub(crate) fn rustred_max_rank_physical() -> Option<usize> {
    static RANK: std::sync::LazyLock<Option<usize>> = std::sync::LazyLock::new(|| {
        std::env::var("BETA_RUSTRED_MAX_RANK_PHYSICAL")
            .ok()
            .and_then(|value| value.parse().ok())
    });
    *RANK
}

/// Diagnostic (`BETA_RUSTRED_DUMP_DIR`): write the p-integrals of a RustRed
/// reduction slower than one second as one vakint expression per line, with
/// its reduction time, so they can be replayed against Forcer.
fn dump_slow_rustred_requests(requests: &[super::PropagatorRequest], elapsed: std::time::Duration) {
    static DIRECTORY: std::sync::LazyLock<Option<std::path::PathBuf>> =
        std::sync::LazyLock::new(|| std::env::var_os("BETA_RUSTRED_DUMP_DIR").map(Into::into));
    static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let Some(directory) = DIRECTORY.as_ref() else {
        return;
    };
    if elapsed.as_secs_f64() < 1.0 {
        return;
    }
    let index = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let mut text = format!("# rustred_seconds {:.3}\n", elapsed.as_secs_f64());
    for request in requests {
        let powers = request
            .powers
            .iter()
            .map(|power| power.to_string())
            .collect::<Vec<_>>()
            .join(",");
        text.push_str(&format!(
            "({})*topo({}({},{powers}))\n",
            request.canonical_numerator.to_canonical_string(),
            request.family.name,
            request.external.to_canonical_string(),
        ));
    }
    let _ = std::fs::write(directory.join(format!("closure-{index:05}.txt")), text);
}

/// `VAKINT_P_FOLD_CHECK=<file>`: run every RustRed closure job today's way
/// and folded (`VAKINT_P_FOLD_CHECK_MODE`, default `basis`), compare the
/// finalized values and append one JSON line per job to the file.
fn fold_check_log() -> Option<&'static std::path::PathBuf> {
    static LOG: std::sync::LazyLock<Option<std::path::PathBuf>> =
        std::sync::LazyLock::new(|| std::env::var_os("VAKINT_P_FOLD_CHECK").map(Into::into));
    LOG.as_ref()
}

fn write_fold_record(path: &std::path::Path, record: serde_json::Value) {
    use std::io::Write;
    static LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
    let _guard = LOCK.lock();
    if let Ok(mut file) = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
    {
        let _ = writeln!(file, "{record}");
    }
}

fn status_kib(field: &str) -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|status| {
            status.lines().find_map(|line| {
                line.strip_prefix(field)
                    .and_then(|rest| rest.trim().trim_end_matches("kB").trim().parse().ok())
            })
        })
        .unwrap_or(0)
}

/// Restart the peak (`VmHWM`) at the current RSS and return the RSS (KiB).
fn reset_peak() -> u64 {
    let _ = std::fs::write("/proc/self/clear_refs", "5");
    status_kib("VmRSS:")
}

impl MassiveClosure {
    /// Reduce `requests` in `fold` mode and finalize them: the value, and
    /// the times, counts and memory of the two stages.
    fn reduce_and_finalize(
        &self,
        vakint: &Vakint,
        settings: &crate::VakintSettings,
        requests: &[super::PropagatorRequest],
        head: &Atom,
        mode: RustRedMasterMode,
        fold: crate::rustred_evaluation::p_integral::FoldMode,
    ) -> Result<(Atom, serde_json::Value), VakintError> {
        use crate::rustred_evaluation::p_integral as rr;
        let base = reset_peak();
        let started = std::time::Instant::now();
        let reduced = match fold {
            rr::FoldMode::Off => {
                rr::reduce_requests_unfolded(vakint, settings, requests, None, true)?
            }
            fold => rr::reduce_requests_folded(vakint, settings, requests, None, fold, true)?,
        };
        let reduce_seconds = started.elapsed().as_secs_f64();
        let stats = rr::last_stats();
        let reduce_peak = status_kib("VmHWM:").saturating_sub(base);
        let base = reset_peak();
        let finalizing = std::time::Instant::now();
        let mut value = Atom::Zero;
        let mut values = 0;
        for ((loops, _), list) in reduced {
            values += list.len();
            value += self.finalize(&list, loops, head, settings, mode)?;
        }
        let finalize_seconds = finalizing.elapsed().as_secs_f64();
        let finalize_peak = status_kib("VmHWM:").saturating_sub(base);
        Ok((
            value,
            serde_json::json!({
                "fold": format!("{fold:?}"),
                "reduce_seconds": reduce_seconds,
                "split_seconds": stats.split_seconds,
                "lower_seconds": stats.lower_seconds,
                "prepare_seconds": stats.prepare_seconds,
                "library_seconds": stats.library_seconds,
                "render_seconds": stats.render_seconds,
                "channels_split": stats.channels,
                "reductions": stats.reductions,
                "pairs": stats.pairs,
                "values": values,
                "finalize_seconds": finalize_seconds,
                "reduce_peak_kib": reduce_peak,
                "finalize_peak_kib": finalize_peak,
            }),
        ))
    }

    /// [`Self::reduce_and_finalize`] today's way and folded, compared.
    fn fold_check(
        &self,
        vakint: &Vakint,
        settings: &crate::VakintSettings,
        requests: &[super::PropagatorRequest],
        head: &Atom,
        mode: RustRedMasterMode,
    ) -> Result<(Atom, serde_json::Value), VakintError> {
        use crate::rustred_evaluation::p_integral::FoldMode;
        let fold = match std::env::var("VAKINT_P_FOLD_CHECK_MODE").as_deref() {
            Ok("coefficient") => FoldMode::Coefficient,
            Ok("basis") => FoldMode::Basis,
            Ok("off") => FoldMode::Off,
            _ => FoldMode::DEFAULT,
        };
        let order = std::env::var("VAKINT_P_FOLD_CHECK_ORDER").unwrap_or_default();
        if order == "only" {
            // One mode alone (cold caches), for timing; the caller compares
            // the canonical values of separate runs.
            let (value, stats) =
                self.reduce_and_finalize(vakint, settings, requests, head, mode, fold)?;
            let record = serde_json::json!({
                "pid": std::process::id(),
                "requests": requests.len(),
                "families": requests.iter().map(|r| r.family.name).collect::<Vec<_>>(),
                "loops": requests.first().map(|r| r.family.loops),
                "epsilon_terms": settings.number_of_terms_in_epsilon_expansion,
                "only": stats,
            });
            return Ok((value, record));
        }
        let (today, folded) = if order == "fold-first" {
            let folded = self.reduce_and_finalize(vakint, settings, requests, head, mode, fold)?;
            let today =
                self.reduce_and_finalize(vakint, settings, requests, head, mode, FoldMode::Off)?;
            (today, folded)
        } else {
            let today =
                self.reduce_and_finalize(vakint, settings, requests, head, mode, FoldMode::Off)?;
            let folded = self.reduce_and_finalize(vakint, settings, requests, head, mode, fold)?;
            (today, folded)
        };
        let comparing = std::time::Instant::now();
        let difference = (&today.0 - &folded.0).expand();
        let equal = difference.is_zero() || difference.together().is_zero();
        let identical =
            today.0.expand().to_canonical_string() == folded.0.expand().to_canonical_string();
        let record = serde_json::json!({
            "pid": std::process::id(),
            "requests": requests.len(),
            "families": requests.iter().map(|r| r.family.name).collect::<Vec<_>>(),
            "loops": requests.first().map(|r| r.family.loops),
            "epsilon_terms": settings.number_of_terms_in_epsilon_expansion,
            "today": today.1,
            "folded": folded.1,
            "equal": equal,
            "identical_expanded_string": identical,
            "compare_seconds": comparing.elapsed().as_secs_f64(),
        });
        if !equal {
            return Err(super::invalid(format!(
                "VAKINT_P_FOLD_CHECK: the folded closure differs from today's: {record}"
            )));
        }
        Ok((today.0, record))
    }

    fn validated_head(&self, settings: &crate::VakintSettings) -> Result<Atom, VakintError> {
        let (head, _) = external_head(&self.momentum)?;
        if self.mass_squared.contains_symbol(S.k)
            || self.mass_squared.contains_symbol(S.p)
            || self
                .mass_squared
                .contains_symbol(vk_symbol!(&settings.epsilon_symbol))
        {
            return Err(super::invalid(
                "the closing mass must be independent of momenta and epsilon",
            ));
        }
        if self.power <= 0 || self.mass_squared.is_zero() {
            return Err(super::invalid(
                "the closing propagator needs a positive power and a nonzero mass squared",
            ));
        }
        Ok(head)
    }

    /// [`Self::evaluate_batch`] of `(closure, input.ordinary())`, with every
    /// closure whose numerator the direct lowering covers lowered without
    /// building the ordinary input ([`crate::rustred_evaluation::p_lowering`]).
    /// The values are those of the ordinary path.
    pub fn evaluate_direct_batch(
        vakint: &Vakint,
        settings: &crate::VakintSettings,
        inputs: &[(Self, DirectClosureInput)],
    ) -> Result<Vec<Atom>, VakintError> {
        let mode = lowering_mode();
        if mode == "legacy" {
            let ordinary = inputs
                .iter()
                .map(|(closure, input)| (closure.clone(), input.ordinary()))
                .collect::<Vec<_>>();
            return Self::evaluate_batch(vakint, settings, &ordinary);
        }
        let check = mode == "check";
        let mut results = Vec::with_capacity(inputs.len());
        for (closure, input) in inputs {
            let started = std::time::Instant::now();
            let mut stats = crate::rustred_evaluation::p_lowering::DirectStats::default();
            let direct = settings.with_resource_budget(|| {
                closure.evaluate_direct(vakint, settings, input, &mut stats)
            })?;
            let direct_seconds = started.elapsed().as_secs_f64();
            let (value, declined) = match direct {
                Ok(value) => (value, None),
                Err(declined) => {
                    let value = Self::evaluate_batch(
                        vakint,
                        settings,
                        &[(closure.clone(), input.ordinary())],
                    )?
                    .pop()
                    .unwrap_or_default();
                    (value, Some(declined.0))
                }
            };
            if crate::census::enabled() {
                crate::census::record(
                    "direct_closure",
                    serde_json::json!({
                        "declined": declined,
                        "numerator_summands": crate::census::summands(input.numerator.as_view()),
                        "numerator_bytes": input.numerator.as_view().get_byte_size(),
                        "numerator_terms": stats.numerator_terms,
                        "canonical_terms": stats.canonical_terms,
                        "merged_terms": stats.merged_terms,
                        "distinct_monomials": stats.distinct_monomials,
                        "spectators": stats.spectators,
                        "lowered": stats.lowered,
                        "folded": stats.folded,
                        "read_seconds": stats.read_seconds,
                        "classify_seconds": stats.classify_seconds,
                        "monomial_seconds": stats.monomial_seconds,
                        "scatter_seconds": stats.scatter_seconds,
                        "finish_seconds": stats.finish_seconds,
                        "match_seconds": stats.match_seconds,
                        "reduce_seconds": stats.reduce_seconds,
                        "finalize_seconds": stats.finalize_seconds,
                        "seconds": started.elapsed().as_secs_f64(),
                        "direct_seconds": direct_seconds,
                    }),
                );
            }
            if check && declined.is_none() {
                let checking = std::time::Instant::now();
                let ordinary = Self::evaluate_batch(
                    vakint,
                    settings,
                    &[(closure.clone(), input.ordinary())],
                )?
                .pop()
                .unwrap_or_default();
                let ordinary_seconds = checking.elapsed().as_secs_f64();
                let identical = ordinary == value;
                let equal = identical || (&ordinary - &value).expand().is_zero();
                eprintln!(
                    "LOWERING_CHECK identical={identical} equal={equal} bytes={} direct_seconds={direct_seconds:.3} ordinary_seconds={ordinary_seconds:.3}",
                    input.numerator.as_view().get_byte_size()
                );
                if !identical {
                    let dir = std::env::var("VAKINT_LOWERING_CHECK_DUMP").ok();
                    if let Some(dir) = dir {
                        let n = std::process::id() as u128 * 1_000_000
                            + std::time::SystemTime::now()
                                .duration_since(std::time::UNIX_EPOCH)
                                .map(|d| d.as_micros() % 1_000_000)
                                .unwrap_or(0);
                        let _ = std::fs::write(
                            format!("{dir}/mismatch-{n}.txt"),
                            format!(
                                "closure {:?}\ntopology {}\nrouting {:?}\nnumerator {}\ndirect {}\nordinary {}\n",
                                closure,
                                input.topology,
                                input.routing,
                                input.numerator,
                                value,
                                ordinary
                            ),
                        );
                    }
                }
            } else if check {
                eprintln!("LOWERING_CHECK declined={}", declined.unwrap_or(""));
            }
            results.push(value);
        }
        if std::env::var_os("VAKINT_LOWERING_SUMCHECK").is_some() {
            // Summing check: closures with one closing key, lowered, absorbed
            // and reduced once, against the sum of their values.
            let mut groups = std::collections::BTreeMap::<(Atom, i64, Atom), (LoweredClosure, Atom, usize)>::new();
            for ((closure, input), value) in inputs.iter().zip(&results) {
                let Some(lowered) = closure.lower_direct(vakint, settings, input)? else {
                    continue;
                };
                match groups.entry(lowered.closing_key()) {
                    std::collections::btree_map::Entry::Occupied(mut entry) => {
                        let (sum, total, count) = entry.get_mut();
                        let _ = sum.absorb(lowered);
                        *total += value;
                        *count += 1;
                    }
                    std::collections::btree_map::Entry::Vacant(entry) => {
                        entry.insert((lowered, value.clone(), 1));
                    }
                }
            }
            for (_, (sum, total, count)) in groups {
                if count < 2 {
                    continue;
                }
                let summed = sum.evaluate(vakint, settings)?;
                let total = total.expand();
                let identical = summed == total;
                let equal = identical || (&summed - &total).expand().is_zero();
                eprintln!("LOWERING_SUMCHECK closures={count} identical={identical} equal={equal}");
            }
        }
        Ok(results)
    }

    /// The direct lowering of one closure without its reduction: the folded
    /// RustRed input keyed by canonical integrals, which can be added to the
    /// input of other closures (of this or other diagrams) with the same
    /// closing mass, power and head before one reduction ([`LoweredClosure`]).
    /// `Ok(None)` when the direct lowering declines the input (use
    /// [`Self::evaluate_direct_batch`]).
    pub fn lower_direct(
        &self,
        vakint: &Vakint,
        settings: &crate::VakintSettings,
        input: &DirectClosureInput,
    ) -> Result<Option<LoweredClosure>, VakintError> {
        use crate::rustred_evaluation::p_lowering as direct;
        let head = self.validated_head(settings)?;
        if input.numerator.is_zero() {
            return Ok(Some(LoweredClosure { closure: self.clone(), head, folded: Default::default() }));
        }
        let cut = match direct::match_cut(vakint, settings, &input.topology, &head)? {
            Ok(Some(cut)) => cut,
            Ok(None) => {
                return Ok(Some(LoweredClosure { closure: self.clone(), head, folded: Default::default() }));
            }
            Err(_) => return Ok(None),
        };
        let handle = vakint.p_library()?;
        let mut stats = direct::DirectStats::default();
        match direct::lower_closure(vakint, settings, &handle, &cut, &input.routing, input.numerator.as_view(), &mut stats)? {
            Ok(folded) => Ok(Some(LoweredClosure { closure: self.clone(), head, folded })),
            Err(_) => Ok(None),
        }
    }

    /// One closure through the direct lowering; `Err(Declined)` when the
    /// ordinary path must evaluate it.
    fn evaluate_direct(
        &self,
        vakint: &Vakint,
        settings: &crate::VakintSettings,
        input: &DirectClosureInput,
        stats: &mut crate::rustred_evaluation::p_lowering::DirectStats,
    ) -> Result<Result<Atom, crate::rustred_evaluation::p_lowering::Declined>, VakintError> {
        use crate::rustred_evaluation::p_integral as rr;
        use crate::rustred_evaluation::p_lowering::{self as direct, Declined};
        let head = self.validated_head(settings)?;
        if input.numerator.is_zero() {
            return Ok(Ok(Atom::Zero));
        }
        settings.validate_execution()?;
        if fold_check_log().is_some() {
            return Ok(Err(Declined("the fold check runs the ordinary path")));
        }
        let fold = rr::fold_mode();
        if fold == rr::FoldMode::Off {
            return Ok(Err(Declined("the unfolded reduction")));
        }
        let Some(EvaluationMethod::RustRed(options)) = settings
            .evaluation_order
            .0
            .iter()
            .find(|method| matches!(method, EvaluationMethod::RustRed(_)))
            .cloned()
        else {
            return Ok(Err(Declined("no RustRed method")));
        };
        let first = settings.evaluation_order.0.iter().find_map(|method| match method {
            EvaluationMethod::RustRed(options) => Some(options.master_mode),
            EvaluationMethod::Forcer(options) => Some(options.master_mode),
            _ => None,
        });
        if first != Some(options.master_mode) {
            return Ok(Err(Declined("a Forcer method precedes RustRed")));
        }
        let matching = std::time::Instant::now();
        let cut = match direct::match_cut(vakint, settings, &input.topology, &head)? {
            Ok(Some(cut)) => cut,
            Ok(None) => return Ok(Ok(Atom::Zero)),
            Err(declined) => return Ok(Err(declined)),
        };
        stats.match_seconds += matching.elapsed().as_secs_f64();
        let handle = vakint.p_library()?;
        let folded = match direct::lower_closure(
            vakint,
            settings,
            &handle,
            &cut,
            &input.routing,
            input.numerator.as_view(),
            stats,
        )? {
            Ok(folded) => folded,
            Err(declined) => return Ok(Err(declined)),
        };
        let reducing = std::time::Instant::now();
        let reduced = rr::reduce_folded_channels(vakint, settings, &handle, folded, None, fold, true)?;
        stats.reduce_seconds += reducing.elapsed().as_secs_f64();
        let finalizing = std::time::Instant::now();
        let mode = options.master_mode;
        let mut result = Atom::Zero;
        let mut numerical = Atom::Zero;
        for ((loops, _), values) in reduced {
            let value = self.finalize(&values, loops, &head, settings, mode)?;
            if mode == RustRedMasterMode::Numerical {
                numerical += value;
            } else {
                result += value;
            }
        }
        if !numerical.is_zero() {
            result += vakint.resolve_rustred_numerical(numerical, settings)?;
        }
        let result = result.expand();
        stats.finalize_seconds += finalizing.elapsed().as_secs_f64();
        Ok(Ok(result))
    }

    /// Evaluate independent closures in bounded batches. Forcer reductions of
    /// exact scalar requests share FORCER invocations; RustRed reduces each
    /// request through its p library. Every closure keeps its own mass,
    /// power, external momentum and Laurent expansion. Spectator
    /// coefficients, including epsilon poles, remain attached to their
    /// original homogeneous channels until that expansion.
    pub fn evaluate_batch(
        vakint: &Vakint,
        settings: &crate::VakintSettings,
        inputs: &[(Self, Atom)],
    ) -> Result<Vec<Atom>, VakintError> {
        let mut results = Vec::with_capacity(inputs.len());
        for chunk in inputs.chunks(32) {
            results
                .extend(settings.with_resource_budget(|| {
                    Self::evaluate_batch_inner(vakint, settings, chunk)
                })?);
        }
        Ok(results)
    }

    fn evaluate_batch_inner(
        vakint: &Vakint,
        settings: &crate::VakintSettings,
        inputs: &[(Self, Atom)],
    ) -> Result<Vec<Atom>, VakintError> {
        let _stage = crate::stage_timing::StageSpan::new("vakint.massive.closure_batch");
        let mode = settings
            .evaluation_order
            .0
            .iter()
            .find_map(|method| match method {
                EvaluationMethod::RustRed(options) => Some(options.master_mode),
                EvaluationMethod::Forcer(options) => Some(options.master_mode),
                _ => None,
            })
            .unwrap_or(RustRedMasterMode::Symbolic);
        let mut results = vec![Atom::Zero; inputs.len()];
        let mut numerical = vec![Atom::Zero; inputs.len()];
        let mut heads = Vec::with_capacity(inputs.len());
        let mut forcer_requests = Vec::new();
        let mut forcer_owners = Vec::new();
        // RustRed requests borrow their match; keep each closure's matches so
        // that all of its requests are reduced together.
        let mut rustred_jobs: Vec<Vec<(crate::ReplacementRules, Atom)>> =
            (0..inputs.len()).map(|_| Vec::new()).collect();
        let mut rustred_mode = None;
        let census = crate::census::enabled();
        let census_started = std::time::Instant::now();
        let mut census_counts = [0_usize; 6];
        let mut census_families = std::collections::BTreeSet::new();
        for (owner, (closure, input)) in inputs.iter().enumerate() {
            let head = closure.validated_head(settings)?;
            heads.push(head.clone());
            if input.is_zero() {
                continue;
            }
            settings.validate_execution()?;
            let mut topologies = Vec::new();
            input.visitor(&mut |part| {
                if let AtomView::Fun(function) = part
                    && function.get_symbol() == S.topo
                {
                    topologies.push(part.to_owned());
                    return false;
                }
                true
            });
            topologies.sort();
            topologies.dedup();
            for (topology, numerator) in input.coefficient_list::<i16>(&topologies) {
                if topology == Atom::one() {
                    let value = closure.finalize(
                        &[ReducedPropagatorValue {
                            scale_power: 0,
                            reduced: Vakint::convert_to_dot_notation(numerator.as_view()),
                            spectator: None,
                        }],
                        0,
                        &head,
                        settings,
                        mode,
                    )?;
                    if mode == RustRedMasterMode::Numerical {
                        numerical[owner] += value;
                    } else {
                        results[owner] += value;
                    }
                    continue;
                }
                let mut expression = VakintExpression::try_from(topology * numerator)?;
                if census {
                    census_counts[0] += 1;
                    census_counts[1] += expression.0.len();
                }
                for term in &mut expression.0 {
                    if !term.normalize_propagator()? {
                        continue;
                    }
                    vakint.validate_settings(settings)?;
                    let external = term.propagator_external()?.ok_or_else(|| {
                        super::invalid("a massive closure requires a massless p-integral")
                    })?;
                    if external != head {
                        return Err(super::invalid(format!(
                            "the cut graph uses {external}, but the closing line uses {head}"
                        )));
                    }
                    term.numerator = scalar_numerator(term.numerator.as_view(), &head)?;
                    let matched = term.prepare_integral_specs(vakint, settings, None)?;
                    let method = vakint.evaluation_method(settings, &matched.canonical_topology);
                    let request = matched
                        .propagator_request(term.numerator.as_view())
                        .expect("a p-integral match has a propagator request");
                    if census {
                        census_counts[2] += 1;
                        census_counts[3] += crate::census::summands(request.canonical_numerator);
                        census_counts[4] += request.canonical_numerator.get_byte_size();
                        census_families
                            .insert((request.family.name, request.powers.to_vec()));
                    }
                    match method {
                        Some(EvaluationMethod::Forcer(options)) => {
                            if options.master_mode != mode {
                                return Err(super::invalid(
                                    "batched massive closures need one master mode",
                                ));
                            }
                            if let Some(input) = forcer::forcer_input_with_epsilon(
                                request.family,
                                request.powers,
                                request.canonical_numerator,
                                &request.external,
                                vk_symbol!(&settings.epsilon_symbol),
                            )? && !input.parts.is_empty()
                            {
                                forcer_requests.push(input);
                                forcer_owners.push((owner, request.family.loops));
                            }
                        }
                        Some(EvaluationMethod::RustRed(options))
                            if rustred_max_rank().is_some_and(|max| {
                                loop_rank(request.canonical_numerator) > max
                            }) =>
                        {
                            if options.master_mode != mode {
                                return Err(super::invalid(
                                    "batched massive closures need one master mode",
                                ));
                            }
                            if let Some(input) = forcer::forcer_input_with_epsilon(
                                request.family,
                                request.powers,
                                request.canonical_numerator,
                                &request.external,
                                vk_symbol!(&settings.epsilon_symbol),
                            )? && !input.parts.is_empty()
                            {
                                forcer_requests.push(input);
                                forcer_owners.push((owner, request.family.loops));
                            }
                        }
                        Some(EvaluationMethod::RustRed(options)) => {
                            if *rustred_mode.get_or_insert(options.master_mode)
                                != options.master_mode
                            {
                                return Err(super::invalid(
                                    "batched massive closures need one master mode",
                                ));
                            }
                            drop(request);
                            let numerator = term.numerator.clone();
                            rustred_jobs[owner].push((matched, numerator));
                        }
                        _ => {
                            return Err(super::invalid(
                                "massive closure of a p-integral requires the RustRed or Forcer \
                                 backend",
                            ));
                        }
                    }
                }
            }
        }
        let census_matching = census_started.elapsed().as_secs_f64();
        if let Some(rustred_mode) = rustred_mode {
            for (owner, jobs) in rustred_jobs.iter().enumerate() {
                if jobs.is_empty() {
                    continue;
                }
                let requests = jobs
                    .iter()
                    .map(|(matched, numerator)| {
                        matched
                            .propagator_request(numerator.as_view())
                            .expect("a p-integral match has a propagator request")
                    })
                    .collect::<Vec<_>>();
                // One closure has one external head; its values are summed
                // per loop count before the radial integration.
                if let Some(log) = fold_check_log() {
                    // Paired check: today's path, then the folded one, on the
                    // same requests; today's value is returned.
                    let (value, record) = inputs[owner].0.fold_check(
                        vakint,
                        settings,
                        &requests,
                        &heads[owner],
                        rustred_mode,
                    )?;
                    write_fold_record(log, record);
                    if rustred_mode == RustRedMasterMode::Numerical {
                        numerical[owner] += value;
                    } else {
                        results[owner] += value;
                    }
                    continue;
                }
                let started = std::time::Instant::now();
                let reduced = crate::rustred_evaluation::p_integral::reduce_requests_factored(
                    vakint, settings, &requests, None,
                )?;
                dump_slow_rustred_requests(&requests, started.elapsed());
                if census {
                    let stats = crate::rustred_evaluation::p_integral::last_stats();
                    crate::census::record(
                        "reduce_folded",
                        stats.census_json(requests.len(), started.elapsed().as_secs_f64()),
                    );
                }
                for ((loops, _), values) in reduced {
                    let value = inputs[owner].0.finalize(
                        &values,
                        loops,
                        &heads[owner],
                        settings,
                        rustred_mode,
                    )?;
                    if rustred_mode == RustRedMasterMode::Numerical {
                        numerical[owner] += value;
                    } else {
                        results[owner] += value;
                    }
                }
            }
        }
        if census {
            crate::census::record(
                "closure_batch",
                serde_json::json!({
                    "closures": inputs.len(),
                    "topology_groups": census_counts[0],
                    "vakint_terms": census_counts[1],
                    "requests": census_counts[2],
                    "request_summands": census_counts[3],
                    "request_bytes": census_counts[4],
                    "distinct_family_powers": census_families.len(),
                    "input_summands": inputs.iter().map(|(_, input)| crate::census::summands(input.as_view())).sum::<usize>(),
                    "input_bytes": inputs.iter().map(|(_, input)| input.as_view().get_byte_size()).sum::<usize>(),
                    "matching_seconds": census_matching,
                    "seconds": census_started.elapsed().as_secs_f64(),
                }),
            );
        }
        let values = forcer::reduce_inputs(&forcer_requests, settings)?;
        for ((owner, loops), values) in forcer_owners.iter().copied().zip(values) {
            let value = inputs[owner]
                .0
                .finalize(&values, loops, &heads[owner], settings, mode)?;
            if mode == RustRedMasterMode::Numerical {
                numerical[owner] += value;
            } else {
                results[owner] += value;
            }
        }
        for ((result, numerical), (_, input)) in results.iter_mut().zip(numerical).zip(inputs) {
            if input.is_zero() {
                continue;
            }
            if !numerical.is_zero() {
                *result += vakint.resolve_rustred_numerical(numerical, settings)?;
            }
            *result = result.expand();
        }
        Ok(results)
    }

    /// Evaluate scalar p-integrals and loop-free powers of `dot(p,p)`, including
    /// sums of both. The p backend the evaluation order selects for a
    /// p-integral (RustRed with a p library covering its loop count, else
    /// Forcer) reduces it; loop-free terms need no external backend. Tensor
    /// dependence on the hard momentum is rejected rather than discarded.
    ///
    /// The Laurent depth counts the closing loop as well: a cut with `L` loops
    /// is expanded through `number_of_terms_in_epsilon_expansion - L - 2`.
    /// The master mode of the selecting backend gives retained or expanded
    /// output; loop-free terms use the mode of the first RustRed or Forcer
    /// entry of the order, else symbolic expansion.
    pub fn evaluate(
        &self,
        vakint: &Vakint,
        settings: &crate::VakintSettings,
        input: AtomView,
    ) -> Result<Atom, VakintError> {
        settings.with_resource_budget(|| {
            let (head, _) = external_head(&self.momentum)?;
            if self.mass_squared.contains_symbol(S.k)
                || self.mass_squared.contains_symbol(S.p)
                || self
                    .mass_squared
                    .contains_symbol(vk_symbol!(&settings.epsilon_symbol))
            {
                return Err(super::invalid(
                    "the closing mass must be independent of momenta and epsilon",
                ));
            }
            if self.power <= 0 || self.mass_squared.is_zero() {
                return Err(super::invalid(
                    "the closing propagator needs a positive power and a nonzero mass squared",
                ));
            }
            if input.is_zero() {
                return Ok(Atom::Zero);
            }
            settings.validate_execution()?;
            let mode = settings
                .evaluation_order
                .0
                .iter()
                .find_map(|method| match method {
                    EvaluationMethod::RustRed(options) => Some(options.master_mode),
                    EvaluationMethod::Forcer(options) => Some(options.master_mode),
                    _ => None,
                })
                .unwrap_or(RustRedMasterMode::Symbolic);

            // VakintExpression intentionally requires a topology in every term.
            // Separate the local remainder here so closing a tree or a local
            // counterterm also receives its final massive loop integration.
            let mut topologies = Vec::new();
            input.visitor(&mut |part| {
                if let AtomView::Fun(function) = part
                    && function.get_symbol() == S.topo
                {
                    topologies.push(part.to_owned());
                    return false;
                }
                true
            });
            topologies.sort();
            topologies.dedup();
            // Terms finalized in `Numerical` mode get their constants
            // substituted together at the end.
            let mut result = Atom::Zero;
            let mut numerical = Atom::Zero;
            for (topology, numerator) in input.coefficient_list::<i16>(&topologies) {
                if topology == Atom::one() {
                    let value = self.finalize(
                        &[ReducedPropagatorValue {
                            scale_power: 0,
                            reduced: Vakint::convert_to_dot_notation(numerator.as_view()),
                            spectator: None,
                        }],
                        0,
                        &head,
                        settings,
                        mode,
                    )?;
                    if mode == RustRedMasterMode::Numerical {
                        numerical += value;
                    } else {
                        result += value;
                    }
                    continue;
                }
                let mut expression = VakintExpression::try_from(topology * numerator)?;
                for term in &mut expression.0 {
                    if !term.normalize_propagator()? {
                        continue;
                    }
                    vakint.validate_settings(settings)?;
                    let external = term.propagator_external()?.ok_or_else(|| {
                        super::invalid("a massive closure requires a massless p-integral")
                    })?;
                    if external != head {
                        return Err(super::invalid(format!(
                            "the cut graph uses {external}, but the closing line uses {head}"
                        )));
                    }
                    term.numerator = scalar_numerator(term.numerator.as_view(), &head)?;
                    let matched = term.prepare_integral_specs(vakint, settings, None)?;
                    let method = vakint.evaluation_method(settings, &matched.canonical_topology);
                    let term_mode = match method {
                        Some(EvaluationMethod::Forcer(options)) => options.master_mode,
                        Some(EvaluationMethod::RustRed(options)) => options.master_mode,
                        _ => {
                            return Err(super::invalid(
                                "massive closure of a p-integral requires the RustRed or Forcer \
                                 backend",
                            ));
                        }
                    };
                    let request = matched
                        .propagator_request(term.numerator.as_view())
                        .expect("a p-integral match has a propagator request");
                    let values = if matches!(method, Some(EvaluationMethod::Forcer(_))) {
                        forcer::reduce_p(
                            request.family,
                            request.powers,
                            request.canonical_numerator,
                            &request.external,
                            settings,
                        )?
                    } else {
                        crate::rustred_evaluation::p_integral::reduce_requests_factored(
                            vakint,
                            settings,
                            std::slice::from_ref(&request),
                            None,
                        )?
                        .into_values()
                        .flatten()
                        .collect()
                    };
                    let value =
                        self.finalize(&values, request.family.loops, &head, settings, term_mode)?;
                    if term_mode == RustRedMasterMode::Numerical {
                        numerical += value;
                    } else {
                        result += value;
                    }
                }
            }
            if !numerical.is_zero() {
                result += vakint.resolve_rustred_numerical(numerical, settings)?;
            }
            Ok(result.expand())
        })
    }

    /// The coefficients of the powers of the closing square `dot(head,
    /// head)` in the expanded `expression`, by power. They must be free of
    /// loop momenta and of the hard momentum and use only known master
    /// symbols.
    fn radial_terms(expression: &Atom, head: &Atom) -> Result<Vec<(i64, Atom)>, VakintError> {
        let square = function!(S.dot, head, head);
        let mut terms = Vec::new();
        for (monomial, coefficient) in
            expression.coefficient_list::<i16>(std::slice::from_ref(&square))
        {
            let radial_power = if monomial == Atom::one() {
                0
            } else if monomial == square {
                1
            } else if let AtomView::Pow(power) = monomial.as_view()
                && power.get_base_exp().0 == square.as_view()
                && let Some(power) = crate::get_integer_from_atom(power.get_base_exp().1)
            {
                power
            } else {
                return Err(super::invalid(format!(
                    "the closing momentum needs integer powers of its square: {monomial}"
                )));
            };
            Self::check_scalar(&coefficient, head)?;
            terms.push((radial_power, coefficient));
        }
        Ok(terms)
    }

    /// `coefficient` is free of loop momenta and of the hard momentum and
    /// uses only known master symbols.
    fn check_scalar(coefficient: &Atom, head: &Atom) -> Result<(), VakintError> {
        let mut hard_dependence = coefficient.contains_symbol(S.k);
        coefficient.visitor(&mut |part| {
            if let AtomView::Fun(function) = part
                && function.get_symbol() == S.p
                && let AtomView::Fun(external) = head.as_view()
                && function.iter().next() == external.iter().next()
            {
                hard_dependence = true;
            }
            !hard_dependence
        });
        if hard_dependence {
            return Err(VakintError::InvalidNumerator(format!(
                "massive closure requires scalar dependence on the hard momentum: {coefficient}"
            )));
        }
        super::compare::check_symbols(coefficient)
    }

    fn finalize(
        &self,
        values: &[ReducedPropagatorValue],
        cut_loops: usize,
        head: &Atom,
        settings: &crate::VakintSettings,
        mode: RustRedMasterMode,
    ) -> Result<Atom, VakintError> {
        let epsilon = vk_symbol!(&settings.epsilon_symbol);
        let loops = cut_loops as i64 + 1;
        let target = settings.number_of_terms_in_epsilon_expansion - loops - 1;
        let mut result = Atom::Zero;
        // The radial prefactor depends only on the scale power.
        let mut prefactors = std::collections::BTreeMap::<i64, Atom>::new();
        for value in values {
            if mode != RustRedMasterMode::Retain
                && let Some(spectator) = &value.spectator
                && let Some(finalized) = self.finalize_factored(
                    &value.reduced,
                    spectator,
                    value.scale_power,
                    cut_loops,
                    head,
                    settings,
                    &mut prefactors,
                )?
            {
                result += finalized;
                continue;
            }
            // The backend has already reduced a scalar numerator. Its PG(x,y)
            // factors are scalar master tags, not two-index vectors for the
            // generic tensor notation converter.
            let reduced = value.value().expand();
            for (radial_power, coefficient) in Self::radial_terms(&reduced, head)? {
                let scale_power = value.scale_power + radial_power;
                let prefactor =
                    self.memoised_prefactor(&mut prefactors, settings, cut_loops, scale_power)?;
                let coefficient = Atom::num(self.radial_sign(radial_power)) * coefficient;
                result += if mode == RustRedMasterMode::Retain {
                    simplify_real((prefactor * coefficient).as_view())
                } else {
                    laurent::expand(prefactor, &coefficient, epsilon, target)?
                };
            }
        }
        if mode != RustRedMasterMode::Retain {
            crate::truncation_guard::ensure_known_orders(
                &result,
                epsilon,
                Some(target),
                "Forcer massive closure",
                crate::truncation_guard::plain_name,
            )?;
        }
        let mu = Atom::var(vk_symbol!(&settings.mu_r_sq_symbol));
        Ok(simplify_real(result.as_view())
            .replace(Atom::var(vk_symbol!("logMassiveClosure")).to_pattern())
            .with(function!(Symbol::LOG, &self.mass_squared).to_pattern())
            .replace(Atom::var(vk_symbol!("log_mu_sq")).to_pattern())
            .with(function!(Symbol::LOG, mu).to_pattern()))
    }

    /// The Laurent expansion of one factored value `masters · spectator`,
    /// with the masters and every radial coefficient of the spectator
    /// expanded apart ([`laurent::expand_factored`]); the same atom as the
    /// expansion of the product. `None` when a radial coefficient is not an
    /// admitted rational spectator or an expansion declines, and the caller
    /// expands the product.
    ///
    /// Zero masters give zero before any check, as the product does. For
    /// nonzero invalid inputs, checking the factors separately can report
    /// a different invalid factor first than checking their product.
    #[allow(clippy::too_many_arguments)]
    fn finalize_factored(
        &self,
        masters: &Atom,
        spectator: &Atom,
        scale_power: i64,
        cut_loops: usize,
        head: &Atom,
        settings: &crate::VakintSettings,
        prefactors: &mut std::collections::BTreeMap<i64, Atom>,
    ) -> Result<Option<Atom>, VakintError> {
        let epsilon = vk_symbol!(&settings.epsilon_symbol);
        let loops = cut_loops as i64 + 1;
        let target = settings.number_of_terms_in_epsilon_expansion - loops - 1;
        if masters.is_zero() {
            return Ok(Some(Atom::Zero));
        }
        let terms = Self::radial_terms(&spectator.expand(), head)?;
        if terms.is_empty() {
            return Ok(Some(Atom::Zero));
        }
        Self::check_scalar(masters, head)?;
        let mut spectators = Vec::with_capacity(terms.len());
        for (radial_power, coefficient) in terms {
            let coefficient = Atom::num(self.radial_sign(radial_power)) * coefficient;
            let Some(prepared) = laurent::Spectator::new(&coefficient, epsilon) else {
                return Ok(None);
            };
            spectators.push((radial_power, prepared));
        }
        let mut result = Atom::Zero;
        for (radial_power, prepared) in spectators {
            let prefactor = self.memoised_prefactor(
                prefactors,
                settings,
                cut_loops,
                scale_power + radial_power,
            )?;
            let Some(value) =
                laurent::expand_factored(prefactor, masters, &prepared, epsilon, target)?
            else {
                return Ok(None);
            };
            result += value;
        }
        Ok(Some(result))
    }

    /// The sign `(-1)^(power + radial_power)` of a radial term.
    fn radial_sign(&self, radial_power: i64) -> i64 {
        if (self.power + radial_power) % 2 == 0 {
            1
        } else {
            -1
        }
    }

    /// [`Self::radial_prefactor`], computed once per scale power.
    fn memoised_prefactor<'p>(
        &self,
        prefactors: &'p mut std::collections::BTreeMap<i64, Atom>,
        settings: &crate::VakintSettings,
        cut_loops: usize,
        scale_power: i64,
    ) -> Result<&'p Atom, VakintError> {
        if let std::collections::btree_map::Entry::Vacant(entry) = prefactors.entry(scale_power) {
            entry.insert(self.radial_prefactor(settings, cut_loops, scale_power)?);
        }
        Ok(&prefactors[&scale_power])
    }

    fn radial_prefactor(
        &self,
        settings: &crate::VakintSettings,
        cut_loops: usize,
        scale_power: i64,
    ) -> Result<Atom, VakintError> {
        let epsilon = Atom::var(vk_symbol!(&settings.epsilon_symbol));
        let loops = cut_loops as i64 + 1;
        let normalization = settings
            .get_staged_integral_normalization_factor_atom()?
            .replace(S.n_loops.to_pattern())
            .with(Atom::num(loops).to_pattern());
        let measure = vk_parse!(&format!(
            "(𝑖*𝜋^(2-({eps})))^{loops}\
             *(Gamma(1+{eps})*Gamma(1-{eps})^2/Gamma(2-2*{eps}))^{cut_loops}\
             *exp(-logMassiveClosure)^({eps}*{loops})",
            eps = settings.epsilon_symbol,
        ))
        .map_err(VakintError::SymbolicaError)?;
        // ∫ d^d q_E/π^(d/2) (q_E²)^(w-l ε)/(q_E²+M²)^a
        // = (M²)^(2+w-a-(l+1)ε) Γ(2+w-(l+1)ε)
        //   Γ(a-2-w+(l+1)ε)/(Γ(2-ε) Γ(a)).
        let gamma_ratio = gamma().call(Atom::num(2 + scale_power) - Atom::num(loops) * &epsilon)
            * gamma().call(Atom::num(self.power - 2 - scale_power) + Atom::num(loops) * &epsilon)
            / (gamma().call(Atom::num(2) - &epsilon) * gamma().call(Atom::num(self.power)));
        Ok(Vakint::normalize_gamma_functions(
            (normalization
                * measure
                * self
                    .mass_squared
                    .pow(Atom::num(2 + scale_power - self.power))
                * gamma_ratio)
                .as_view(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use std::sync::LazyLock;

    use ahash::HashMap;
    use symbolica::domains::float::Float;

    use super::*;
    use crate::{EvaluationOrder, ForcerOptions, LoopNormalizationFactor, VakintSettings};

    static VAKINT: LazyLock<Vakint> = LazyLock::new(|| Vakint::new().unwrap());

    fn parse(text: &str) -> Atom {
        Vakint::initialize_vakint_symbols();
        vk_parse!(text).unwrap()
    }

    fn settings(terms: i64) -> VakintSettings {
        VakintSettings {
            evaluation_order: EvaluationOrder::empty(),
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
            n_cores: 1,
            number_of_terms_in_epsilon_expansion: terms,
            use_dot_product_notation: true,
            ..Default::default()
        }
    }

    fn closure(power: i64) -> MassiveClosure {
        MassiveClosure {
            momentum: parse("p(1)"),
            mass_squared: Atom::one(),
            power,
        }
    }

    fn assert_equal(left: &Atom, right: &Atom) {
        assert!(
            (left - right).expand().is_zero(),
            "left: {left}\nright: {right}\ndifference: {}",
            (left - right).expand()
        );
    }

    #[test]
    fn massive_closure_preserves_tadpole_poles_finite_and_higher_orders() {
        let settings = settings(5);
        let value = closure(2)
            .evaluate(&VAKINT, &settings, Atom::one().as_view())
            .unwrap();
        let expected = parse("1/ε+𝜋^2*ε/12+𝜋^4*ε^3/160")
            - symbolica::transcendental::zeta().call(3) * parse("ε^2/3");
        // Both m² and μ² are one; no scale logarithms survive.
        let value = value
            .replace(Atom::var(vk_symbol!("mursq")).to_pattern())
            .with(Atom::one().to_pattern());
        assert_equal(&value, &expected);
    }

    #[test]
    fn massive_closure_matches_the_independent_vacuum_backend() {
        let mut settings = settings(5);
        for power in [1, 2, 3, 4] {
            let value = closure(power)
                .evaluate(&VAKINT, &settings, Atom::one().as_view())
                .unwrap();
            settings.evaluation_order = EvaluationOrder(vec![EvaluationMethod::RustRed(
                crate::RustRedEvaluationOptions {
                    master_mode: RustRedMasterMode::Symbolic,
                },
            )]);
            let reference = VAKINT
                .evaluate(&settings, parse(&format!("topo(I1L(1,{power}))")).as_view())
                .unwrap();
            settings.evaluation_order = EvaluationOrder::empty();
            assert_equal(&value, &reference);
        }
    }

    #[test]
    fn massive_closure_handles_mass_scaling_and_radial_numerators() {
        let settings = settings(4);
        let massive = MassiveClosure {
            mass_squared: parse("mass_squared"),
            ..closure(3)
        };
        // q²/(q²-m²)^3 = 1/(q²-m²)^2 + m²/(q²-m²)^3.
        let with_numerator = massive
            .evaluate(&VAKINT, &settings, parse("dot(p(1),p(1))").as_view())
            .unwrap();
        let lower_power = MassiveClosure {
            power: 2,
            ..massive.clone()
        }
        .evaluate(&VAKINT, &settings, Atom::one().as_view())
        .unwrap();
        let scalar = massive
            .evaluate(&VAKINT, &settings, Atom::one().as_view())
            .unwrap();
        assert_equal(
            &with_numerator,
            &(lower_power + &massive.mass_squared * scalar),
        );
    }

    #[test]
    fn massive_closure_keeps_the_massless_tree_denominator() {
        let settings = settings(3);
        let value = closure(1)
            .evaluate(&VAKINT, &settings, parse("1/dot(p(1),p(1))").as_view())
            .unwrap()
            .replace(Atom::var(vk_symbol!("mursq")).to_pattern())
            .with(Atom::one().to_pattern());
        assert_equal(&value, &parse("1/ε+1+ε*(1+𝜋^2/12)"));
    }

    #[test]
    fn existing_tensor_prepass_projects_a_mixed_mass_vacuum_before_matching() {
        let settings = settings(3);
        let topology = parse("topo(prop(1,edge(1,2),k(1),1,1)*prop(2,edge(2,1),k(1),0,1))");
        let input = parse("dot(k(1),p(2))^2") * &topology;
        let projected = VAKINT.tensor_reduce(&settings, input.as_view()).unwrap();
        let expression = VakintExpression::try_from(projected).unwrap();
        assert_eq!(expression.0.len(), 1);
        assert_eq!(expression.0[0].integral, topology);
        assert!(
            (expression.0[0].numerator.clone() - parse("dot(k(1),k(1))*dot(p(2),p(2))/(4-2*ε)"))
                .together()
                .is_zero()
        );
    }

    #[test]
    fn massive_closure_rejects_open_hard_indices_and_unintegrated_loops() {
        for input in ["p(1,mu)", "dot(p(1),p(2))", "dot(k(1),k(1))"] {
            let error = closure(2)
                .evaluate(&VAKINT, &settings(3), parse(input).as_view())
                .unwrap_err();
            assert!(matches!(error, VakintError::InvalidNumerator(_)), "{error}");
        }
    }

    #[test]
    fn factored_zero_masters_do_not_validate_the_spectator() {
        let settings = settings(3);
        for spectator in ["p(1,mu)", "dot(p(1),p(2))", "dot(k(1),k(1))"] {
            let value = closure(2)
                .finalize(
                    &[ReducedPropagatorValue {
                        scale_power: 0,
                        reduced: Atom::Zero,
                        spectator: Some(parse(spectator)),
                    }],
                    1,
                    &parse("p(1)"),
                    &settings,
                    RustRedMasterMode::Symbolic,
                )
                .unwrap();
            assert!(value.is_zero());
        }
    }

    #[test]
    fn massive_closure_expands_the_exact_bubble_before_truncation() {
        let settings = settings(5);
        let value = closure(2)
            .finalize(
                &[ReducedPropagatorValue {
                    scale_power: 0,
                    reduced: parse("1/ε"),
                    spectator: None,
                }],
                1,
                &parse("p(1)"),
                &settings,
                RustRedMasterMode::Symbolic,
            )
            .unwrap();
        // Direct successive one-loop integration, with both loop measures
        // e^(γε)/(iπ^(d/2)); Γ(2-2ε) cancels before the final expansion.
        let reference =
            parse("exp(2*EulerGamma*ε)*Gamma(1+ε)*Gamma(1-ε)^2*Gamma(2*ε)/(ε*Gamma(2-ε))");
        let reference =
            super::super::finalize::gamma_series(&reference, vk_symbol!("ε"), 2).unwrap();
        let value = value
            .replace(Atom::var(vk_symbol!("mursq")).to_pattern())
            .with(Atom::one().to_pattern());
        assert_equal(&value, &reference);
    }

    #[test]
    fn massive_closure_retains_exact_gamma_dependence_when_requested() {
        let mut settings = settings(5);
        let expanded = closure(2)
            .evaluate(&VAKINT, &settings, Atom::one().as_view())
            .unwrap();
        settings.evaluation_order =
            EvaluationOrder(vec![EvaluationMethod::Forcer(ForcerOptions {
                master_mode: RustRedMasterMode::Retain,
            })]);
        let retained = closure(2)
            .evaluate(&VAKINT, &settings, Atom::one().as_view())
            .unwrap();
        assert!(retained.contains_symbol(gamma()));
        let restored = super::super::finalize::gamma_series(&retained, vk_symbol!("ε"), 3).unwrap();
        assert_equal(&restored, &expanded);
    }

    #[test]
    fn massive_closure_keeps_products_of_scalar_forcer_tags() {
        let settings = settings(5);
        let reduced = parse("PG(1,0)^2*PG(2,0)");
        let value = closure(2)
            .finalize(
                &[ReducedPropagatorValue {
                    scale_power: 0,
                    reduced,
                    spectator: None,
                }],
                2,
                &parse("p(1)"),
                &settings,
                RustRedMasterMode::Symbolic,
            )
            .unwrap();
        let gamma_value = super::super::finalize::pg_value(1, 0, vk_symbol!("ε")).pow(Atom::num(2))
            * super::super::finalize::pg_value(2, 0, vk_symbol!("ε"));
        let reference = closure(2)
            .finalize(
                &[ReducedPropagatorValue {
                    scale_power: 0,
                    reduced: gamma_value,
                    spectator: None,
                }],
                2,
                &parse("p(1)"),
                &settings,
                RustRedMasterMode::Symbolic,
            )
            .unwrap();
        assert_equal(&value, &reference);
    }

    #[test]
    fn massive_closure_checks_master_depth_after_the_extra_loop_pole() {
        let error = closure(2)
            .finalize(
                &[ReducedPropagatorValue {
                    scale_power: 0,
                    reduced: parse("PMaster(no)"),
                    spectator: None,
                }],
                3,
                &parse("p(1)"),
                &settings(13),
                RustRedMasterMode::Symbolic,
            )
            .unwrap_err();
        assert!(
            matches!(error, VakintError::MissingMasterOrders(_)),
            "{error}"
        );
    }

    /// Folded input poles and cancellations must survive the extra radial
    /// pole. Spectators outside the rational fast path must give the same
    /// Laurent series as multiplying the whole coefficient before finalizing.
    #[test]
    fn folded_closures_preserve_poles_cancellations_and_fallback_spectators() {
        use crate::rustred_evaluation::p_integral as rr;

        let mut engine = Vakint::new().unwrap();
        engine
            .inject_p_library(rustred::propagator::testing::one_loop_reference_library().unwrap())
            .unwrap();
        let family = super::super::forcer::fixtures::registry_family("P1L").unwrap();
        let massive = closure(2);
        let head = parse("p(1)");
        let cases: &[&[&str]] = &[
            &["T/ε^2+L/(ε*(3-2*ε))"],
            &["T/(3-2*ε)"],
            &["T/(1-ε)", "-T*ε/(1-ε)", "-T"],
            &["T*exp(ε)/(3-2*ε)"],
            &["T/(muvsq-ε)"],
        ];
        for epsilon_symbol in ["ε", "gammalooprs::{}::ε"] {
            let settings = VakintSettings {
                epsilon_symbol: epsilon_symbol.into(),
                integral_normalization_factor: LoopNormalizationFactor::MSbar,
                ..settings(6)
            };
            let epsilon = Atom::var(vk_symbol!(&settings.epsilon_symbol));
            for coefficients in cases {
                let numerators = coefficients
                    .iter()
                    .map(|text| {
                        parse(text)
                            .replace(parse("ε").to_pattern())
                            .with(epsilon.to_pattern())
                            .replace(parse("T").to_pattern())
                            .with(parse("beta_projection::transverse").to_pattern())
                            .replace(parse("L").to_pattern())
                            .with(parse("beta_projection::longitudinal").to_pattern())
                    })
                    .collect::<Vec<_>>();
                let requests = numerators
                    .iter()
                    .map(|numerator| super::super::PropagatorRequest {
                        family,
                        powers: &[1, 1],
                        canonical_numerator: numerator.as_view(),
                        external: head.clone(),
                    })
                    .collect::<Vec<_>>();
                let finalize = |values: std::collections::BTreeMap<
                    (usize, Atom),
                    Vec<ReducedPropagatorValue>,
                >| {
                    values
                        .into_iter()
                        .fold(Atom::Zero, |sum, ((loops, _), values)| {
                            sum + massive
                                .finalize(
                                    &values,
                                    loops,
                                    &head,
                                    &settings,
                                    RustRedMasterMode::Symbolic,
                                )
                                .unwrap()
                        })
                };
                let reference = finalize(
                    rr::reduce_requests_unfolded(&engine, &settings, &requests, None, false)
                        .unwrap(),
                );
                for mode in [rr::FoldMode::Basis, rr::FoldMode::Coefficient] {
                    let actual = finalize(
                        rr::reduce_requests_folded(&engine, &settings, &requests, None, mode, true)
                            .unwrap(),
                    );
                    assert_equal(&actual, &reference);
                }
                if coefficients.len() == 3 {
                    assert!(reference.expand().is_zero(), "cancellation: {reference}");
                }
            }
        }
    }

    /// Closure-level paired replay of dumped requests (`FOLD_REPLAY_INPUT`,
    /// comma separated `(numerator)*topo(FAMILY(ext,powers))` files, one
    /// closure per file): today's path and the folded one, finalized as a
    /// massive closure `p(1)`, `rstar_mass_squared`, power 1 with the
    /// production epsilon symbol and `FOLD_REPLAY_TERMS` epsilon terms
    /// (default 6), compared exactly. Writes JSON lines to `VAKINT_P_FOLD_CHECK`.
    #[test]
    #[ignore]
    fn fold_replay_dumps() {
        let inputs = std::env::var("FOLD_REPLAY_INPUT").unwrap();
        let log = std::path::PathBuf::from(std::env::var("VAKINT_P_FOLD_CHECK").unwrap());
        let terms = std::env::var("FOLD_REPLAY_TERMS").map_or(6, |t| t.parse().unwrap());
        let mut engine = Vakint::new().unwrap();
        let library = rustred_app::propagator::load_p_library(
            rustred_app::propagator::PLibrarySource::Directory(
                std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4")
                    .unwrap()
                    .into(),
            ),
            &rustred_app::propagator::PLibraryOptions {
                verify_descent: false,
            },
        )
        .unwrap();
        engine.inject_p_library(library).unwrap();
        let settings = VakintSettings {
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
                crate::RustRedEvaluationOptions {
                    master_mode: RustRedMasterMode::Symbolic,
                },
            )]),
            integral_normalization_factor: LoopNormalizationFactor::MSbar,
            epsilon_symbol: "gammalooprs::{}::ε".into(),
            mu_r_sq_symbol: "gammalooprs::{}::mu_r_sq".into(),
            number_of_terms_in_epsilon_expansion: terms,
            use_dot_product_notation: true,
            n_cores: 1,
            ..Default::default()
        };
        let massive = MassiveClosure {
            momentum: parse("p(1)"),
            mass_squared: parse("gammalooprs::{scalar}::rstar_mass_squared"),
            power: 1,
        };
        let head = massive.validated_head(&settings).unwrap();
        for path in inputs.split(',') {
            let text = std::fs::read_to_string(path).unwrap();
            let records = text
                .lines()
                .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
                .map(|line| {
                    let at = line.rfind("*topo(").unwrap();
                    let numerator = parse(&line[..at]);
                    let tail = &line[at + 6..line.len() - 1];
                    let open = tail.find('(').unwrap();
                    let family =
                        super::super::forcer::fixtures::registry_family(&tail[..open]).unwrap();
                    let inner = &tail[open + 1..tail.len() - 1];
                    let mut depth = 0i32;
                    let mut cut = inner.len();
                    for (position, character) in inner.char_indices() {
                        match character {
                            '(' | '{' => depth += 1,
                            ')' | '}' => depth -= 1,
                            ',' if depth == 0 => {
                                cut = position;
                                break;
                            }
                            _ => {}
                        }
                    }
                    let powers: Vec<i64> = inner[cut + 1..]
                        .split(',')
                        .map(|p| p.trim().parse().unwrap())
                        .collect();
                    (family, powers, parse(&inner[..cut]), numerator)
                })
                .collect::<Vec<_>>();
            let requests = records
                .iter()
                .map(
                    |(family, powers, external, numerator)| super::super::PropagatorRequest {
                        family,
                        powers,
                        canonical_numerator: numerator.as_view(),
                        external: external.clone(),
                    },
                )
                .collect::<Vec<_>>();
            let (value, mut record) = massive
                .fold_check(
                    &engine,
                    &settings,
                    &requests,
                    &head,
                    RustRedMasterMode::Symbolic,
                )
                .unwrap();
            record["dump"] = serde_json::json!(path);
            if let Some(directory) = std::env::var_os("FOLD_REPLAY_VALUES") {
                // The canonical expanded value, for comparing separate runs.
                let name = std::path::Path::new(path).file_name().unwrap();
                std::fs::write(
                    std::path::Path::new(&directory).join(name),
                    value.expand().to_canonical_string(),
                )
                .unwrap();
            }
            eprintln!("{record}");
            write_fold_record(&log, record);
        }
    }

    /// The default order after the built-in 1L-3L library (grant U1): with
    /// `[RustRed, Forcer]` the closures of a 1L bubble, the 2L kite and the 3L
    /// NO top go to RustRed and equal the Forcer-only closures, exactly in
    /// `Symbolic` mode and within 1e-30 in `Numerical` mode (40 digits).
    #[test]
    fn massive_closures_through_rustred_equal_the_forcer_closures() {
        let base = VakintSettings {
            run_time_decimal_precision: 40,
            ..settings(5)
        };
        if forcer::forcer_form_version(&base.form_exe_path).is_err() {
            assert_ne!(
                std::env::var("VAKINT_REQUIRE_FORM5").as_deref(),
                Ok("1"),
                "FORM 5 is unavailable, but VAKINT_REQUIRE_FORM5=1 requires it"
            );
            eprintln!("skipping the massive closures through RustRed: FORM 5 is unavailable");
            return;
        }
        let massive = closure(2);
        let mursq = HashMap::from_iter([(
            base.mu_r_sq_symbol.clone(),
            Float::with_val(base.get_binary_precision(), 1),
        )]);
        for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
            let forcer = EvaluationMethod::Forcer(ForcerOptions { master_mode: mode });
            let rustred =
                EvaluationMethod::RustRed(crate::RustRedEvaluationOptions { master_mode: mode });
            let both = VakintSettings {
                evaluation_order: EvaluationOrder(vec![rustred, forcer.clone()]),
                ..base.clone()
            };
            let forcer_only = VakintSettings {
                evaluation_order: EvaluationOrder(vec![forcer]),
                ..base.clone()
            };
            for input in [
                "topo(P1L(p(1),1,1))",
                "topo(P2L(p(1),1,1,1,1,1))",
                "topo(P3L_NO(p(1),1,1,1,1,1,1,1,1))",
            ] {
                let input = parse(input);
                let value = |settings| massive.evaluate(&VAKINT, settings, input.as_view());
                let (actual, expected) = (value(&both).unwrap(), value(&forcer_only).unwrap());
                if mode == RustRedMasterMode::Symbolic {
                    assert_equal(&actual, &expected);
                    continue;
                }
                let number = |value: &Atom| {
                    VAKINT
                        .numerical_evaluation(
                            &both,
                            value.as_view(),
                            &mursq,
                            &Default::default(),
                            None,
                        )
                        .unwrap()
                        .0
                };
                let (agree, diagnostic) =
                    number(&actual).does_approx_match(&number(&expected), None, 1e-30, 5.0);
                assert!(agree, "{input}: {diagnostic}");
            }
        }
        // RustRed reduced them: the built-in program's cache recorded steps.
        let statistics = VAKINT.ordered_cache_statistics();
        assert!(
            statistics
                .programs
                .iter()
                .any(|program| program.data_version == "built-in" && program.native.inserts > 0),
            "{statistics:?}"
        );
    }

    #[test]
    fn massive_closure_runs_forcer_and_retains_local_terms() {
        let mut settings = settings(5);
        if forcer::forcer_form_version(&settings.form_exe_path).is_err() {
            eprintln!("skipping live massive closure: FORM 5 is unavailable");
            return;
        }
        settings.evaluation_order =
            EvaluationOrder(vec![EvaluationMethod::Forcer(ForcerOptions {
                master_mode: RustRedMasterMode::Symbolic,
            })]);
        let bubble = parse("topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,1))");
        let massive = closure(2);
        let expected = massive
            .finalize(
                &[ReducedPropagatorValue {
                    scale_power: 0,
                    reduced: parse("1/ε"),
                    spectator: None,
                }],
                1,
                &massive.momentum,
                &settings,
                RustRedMasterMode::Symbolic,
            )
            .unwrap();
        let local = massive
            .evaluate(&VAKINT, &settings, Atom::one().as_view())
            .unwrap();
        let value = massive
            .evaluate(&VAKINT, &settings, (&bubble + Atom::one()).as_view())
            .unwrap();
        assert_equal(&value, &(expected + local));
        let scaleless = parse("topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,0))");
        assert!(
            massive
                .evaluate(&VAKINT, &settings, scaleless.as_view())
                .unwrap()
                .is_zero()
        );
    }
}

#[cfg(test)]
#[path = "massive_replay_tests.rs"]
mod replay_tests;
