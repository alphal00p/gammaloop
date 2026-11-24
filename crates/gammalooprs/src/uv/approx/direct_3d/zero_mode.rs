//! Matsubara zero-mode projection of completed thermal residue branches.
//!
//! A frequency class is a raised group of bosonic edges without a chemical
//! potential: their signatures agree up to sign and their masses agree, so
//! they carry one Matsubara frequency at one on-shell energy, like the two
//! ring gluons that both carry q. Its total propagator power `nu` is its
//! number of edges. Evaluating the completed thermal CFF at one auxiliary
//! energy x for every class edge, the class's zero mode is the pole
//! `T/x^(2 nu)` times the static remainder of the sector. The other modes are
//! analytic at x = 0, because numerator and contact surfaces are polynomial in
//! x, so no summed pole is deeper. In a vacuum-subtracted sector the vacuum
//! 4D-soft region adds the odd pole `x^-(2 nu - 1)`. A coefficient `[x^-a]`
//! times `OSE(rep)^-a` restores static class propagators. The zero mode's
//! temperature is written as the tag `ZERO_MODE_TEMPERATURE(rep)` of its
//! class, apart from the temperature of every other occupation.
//!
//! Several classes are projected in turn. The summed CFF is analytic in each x
//! on a punctured neighbourhood of zero, so iterated coefficients are the
//! joint Laurent coefficients and the projections commute. A static condition
//! implied by the others restores its propagators without a further `T`. The
//! class propagators are restored once every class is projected, as top-level
//! factors of each residue root; a caller that Taylor-expands the coefficients
//! first multiplies them in only afterwards.
//!
//! A summand that a class has already projected, `c OSE(rep)^-a` with `c`
//! free of class energies, is its own coefficient of order `a` and is not
//! expanded again. Its propagators need not be top-level factors: a later
//! operation's Taylor coefficients may carry them in each term of a sum,
//! through which the order is read without an expansion. Projections are
//! therefore idempotent, and a frozen projected branch never enters the
//! hard-vacuum sector, whose order such a summand only has if an earlier
//! hard-vacuum projection gave it. A summand
//! without any class energy has no zero-mode pole, for instance once an
//! ordinary operation has rewritten every class energy of the family into
//! routed, massive ones.
//!
//! Class energies must be visible as `OSE(e)`. A routed `OnShellEnergy` that
//! still names a class edge is rejected in a summand that needs its
//! expansion, and in a summand without a visible class energy while another
//! summand of the family still shows one. One already rewritten in loop
//! momenta names no edge, so callers keep class edges outside every rescaled
//! subgraph of earlier operations.

use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet};
use std::num::TryFromIntError;
use std::slice;

use color_eyre::Result;
use eyre::Report;
use linnet::half_edge::involution::EdgeIndex;
use symbolica::atom::{Atom, AtomCore, AtomOrView, AtomView, Symbol};
use symbolica::function;
use symbolica::id::Replacement;
use thiserror::Error;
use three_dimensional_reps::{MediumMode, ThermalDistributionCall, ThermalDistributionCallError};

use crate::cff::CutCFFIndex;
use crate::graph::Graph;
use crate::integrands::process::param_builder::ParamBuilderGraph;
use crate::numerator::symbolica_ext::NumeratorAtomExt;
use crate::utils::{GS, ose_atom_from_index};
use crate::uv::Integrands;
use crate::uv::approx::direct_3d::branches::DirectResidueKey;
use crate::uv::thermal_zero_mode::{ZeroModeClassSpec, ZeroModeSector};
use crate::uv::uv_graph::{NoZeroMode, UVE};

use super::{DirectResidueBranches, Local3DDeformation, SoftBranches, class_zero_mode_temperature};

#[cfg(test)]
#[path = "zero_mode_tests.rs"]
mod tests;

/// Why a zero-mode projection was refused. [`Self::TruncationReached`],
/// [`Self::NotLaurent`], [`Self::FrozenHardVacuum`] and
/// [`Self::HardVacuumOccupation`] are internal errors; every other variant
/// names an unsupported class or sector, or an error of the expression
/// machinery.
#[derive(Debug, Error)]
pub(super) enum ZeroModeError {
    #[error("a zero-mode projection needs at least one frequency class")]
    NoClass,
    #[error("a zero-mode frequency class needs at least one edge")]
    EmptyClass,
    #[error("edge {edge} is not an edge of the graph")]
    UnknownEdge { edge: EdgeIndex },
    #[error("edge {edge} {reason}")]
    NoZeroMode { edge: EdgeIndex, reason: NoZeroMode },
    #[error("edges {first} and {second} do not carry one frequency at one on-shell energy")]
    NotRaised { first: EdgeIndex, second: EdgeIndex },
    #[error("zero-mode orders {orders:?} of the class {edges:?} must lie in 1..={max}")]
    InvalidOrders {
        edges: BTreeSet<EdgeIndex>,
        orders: BTreeSet<u32>,
        max: u64,
    },
    #[error(
        "the hard-vacuum sector of the class {edges:?} only has the vacuum-soft order \
         {expected}, not {orders:?}"
    )]
    HardVacuumOrders {
        edges: BTreeSet<EdgeIndex>,
        orders: BTreeSet<u32>,
        expected: u64,
    },
    #[error("zero-mode frequency classes must be disjoint, but edge {edge} is in two")]
    OverlappingClasses { edge: EdgeIndex },
    #[error("the zero-mode class {edges:?} enters a formal function argument: {expression}")]
    FormalArgument {
        edges: BTreeSet<EdgeIndex>,
        expression: Atom,
    },
    #[error(
        "the distribution of edge {edge} depends on the zero-mode class {edges:?} other than \
         at the class energy and finite temperature: energy {energy}, temperature flag \
         {temperature_flag}"
    )]
    DeformedDistribution {
        edges: BTreeSet<EdgeIndex>,
        edge: EdgeIndex,
        energy: Atom,
        temperature_flag: Atom,
    },
    #[error("a routed on-shell energy hides the zero-mode class {edges:?}: {expression}")]
    RoutedEnergy {
        edges: BTreeSet<EdgeIndex>,
        expression: Atom,
    },
    #[error(
        "a shared numerator body is not polynomial in the zero-mode class {edges:?}: {expression}"
    )]
    NotPolynomial {
        edges: BTreeSet<EdgeIndex>,
        expression: Atom,
    },
    /// A bug: the truncation bound let an omitted Laurent term reach a
    /// requested coefficient.
    #[error(
        "internal error: a truncated distribution of the zero-mode class {edges:?} reached \
         [x^-{order}]"
    )]
    TruncationReached {
        edges: BTreeSet<EdgeIndex>,
        order: u32,
    },
    /// A bug: the series left the class energy in a coefficient or a body.
    #[error(
        "internal error: the zero-mode series of the class {edges:?} is not a Laurent \
         polynomial in its energy: {expression}"
    )]
    NotLaurent {
        edges: BTreeSet<EdgeIndex>,
        expression: Atom,
    },
    /// A bug: the vacuum projection of the hard-vacuum sector would discard
    /// the occupations a frozen soft branch retains. Only branches that the
    /// class has already projected may reach it, and those are read off
    /// without an expansion.
    #[error(
        "internal error: a frozen soft branch reached the hard-vacuum sector of the zero-mode \
         class {edges:?} before that class projected it: {expression}"
    )]
    FrozenHardVacuum {
        edges: BTreeSet<EdgeIndex>,
        expression: Atom,
    },
    /// A bug: the hard-vacuum sector reads a summand that the class has
    /// already projected at its order without a vacuum projection, so the
    /// summand must already be hard-vacuum, as the output of an earlier
    /// hard-vacuum projection is.
    #[error(
        "internal error: a summand that the zero-mode class {edges:?} projected before keeps an \
         occupation outside the hard-vacuum sector: {expression}"
    )]
    HardVacuumOccupation {
        edges: BTreeSet<EdgeIndex>,
        expression: Atom,
    },
    #[error(transparent)]
    Distribution(#[from] ThermalDistributionCallError),
    #[error(transparent)]
    Integer(#[from] TryFromIntError),
    #[error(transparent)]
    Expression(Report),
}

/// Branch transformations report through `eyre`; recover a zero-mode error
/// raised inside one of their closures instead of nesting it.
impl From<Report> for ZeroModeError {
    fn from(report: Report) -> Self {
        report.downcast::<Self>().unwrap_or_else(Self::Expression)
    }
}

/// One validated Matsubara frequency class and the Laurent orders requested
/// of it.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(super) struct ZeroModeClass {
    /// The raised edges, their representative and total propagator power
    /// `nu`, one per raised edge.
    spec: ZeroModeClassSpec,
    orders: BTreeSet<u32>,
}

/// Reads off the order of the summands a class has already projected.
struct ProjectedOrders<'c> {
    class: &'c ZeroModeClass,
    representative: Atom,
    /// The numerator families whose shared bodies name a class energy.
    dependent: Vec<Replacement>,
    /// Whether each family call names a class energy at its arguments.
    instantiated: BTreeMap<Atom, bool>,
}

impl<'c> ProjectedOrders<'c> {
    fn new(class: &'c ZeroModeClass, branches: &DirectResidueBranches) -> Self {
        // A shared body carries the energy-map coefficients of every residue
        // row; it may name a class energy only through coefficients that
        // vanish on a re-keyed row, so calls of it are checked at their
        // arguments.
        let dependent = branches
            .numerator_definitions()
            .into_values()
            .filter(|entry| class.names_energy(entry.rhs.as_view()))
            .map(|entry| entry.replacement())
            .collect();
        Self {
            class,
            representative: ose_atom_from_index(class.representative()),
            dependent,
            instantiated: BTreeMap::new(),
        }
    }

    /// The order `a` of a summand `c OSE(rep)^-a` whose `c` names no class
    /// energy, also not through the numerator families it calls at their
    /// arguments, or `None` when it still needs its Laurent expansion. A
    /// summand without any class energy has order zero. The class
    /// propagators need not be top-level factors of the summand
    /// ([`Self::propagator_degree`]).
    fn order(&mut self, summand: AtomView<'_>) -> Option<i64> {
        let degree = self.propagator_degree(summand)?;
        (!self.calls_dependent_family(summand)).then_some(-degree)
    }

    /// The degree `d` of `view` as `c OSE(rep)^d` with `c` free of class
    /// energies, or `None` when `view` depends on a class energy in any
    /// other way. A product adds the degrees of its factors and an integer
    /// power multiplies that of its base; a sum has a degree when all its
    /// terms share it.
    ///
    /// The restored class propagators are top-level factors of every
    /// summand only until a later operation collects its Taylor
    /// coefficients: the static U of a sibling loop on a completed branch
    /// can gather two summands of one coefficient into a sum that each
    /// carry the propagators, as a factor of that coefficient's numerator
    /// family. The degree is read through that structure without expanding
    /// it, so such a summand is still the projected `c OSE(rep)^d`.
    fn propagator_degree(&self, view: AtomView<'_>) -> Option<i64> {
        if view == self.representative.as_view() {
            return Some(1);
        }
        if !self.class.names_energy(view) {
            return Some(0);
        }
        match view {
            AtomView::Mul(product) => product.iter().try_fold(0, |degree, factor| {
                Some(degree + self.propagator_degree(factor)?)
            }),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                Some(i64::try_from(exponent).ok()? * self.propagator_degree(base)?)
            }
            AtomView::Add(sum) => {
                let mut degrees = sum.iter().map(|term| self.propagator_degree(term));
                let first = degrees.next()??;
                degrees.all(|degree| degree == Some(first)).then_some(first)
            }
            _ => None,
        }
    }

    /// The coefficient `c` of a summand `c OSE(rep)^d` that has a
    /// [`Self::propagator_degree`]: the summand with its class propagators
    /// removed wherever that degree found them.
    ///
    /// # Panics
    ///
    /// `view` has no propagator degree.
    fn coefficient<'a>(&self, view: AtomView<'a>) -> AtomOrView<'a> {
        if view == self.representative.as_view() {
            return Atom::num(1).into();
        }
        if !self.class.names_energy(view) {
            return view.into();
        }
        match view {
            AtomView::Mul(product) => {
                Atom::mul_many(product.iter().map(|factor| self.coefficient(factor))).into()
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                self.coefficient(base).pow(exponent).into()
            }
            AtomView::Add(sum) => {
                Atom::add_many(sum.iter().map(|term| self.coefficient(term))).into()
            }
            _ => unreachable!("a summand with a propagator degree has no other class energy"),
        }
    }

    /// Whether `view` calls a numerator family whose body names a class
    /// energy at the arguments of that call. Each call is instantiated once.
    fn calls_dependent_family(&mut self, view: AtomView<'_>) -> bool {
        let Self {
            class,
            dependent,
            instantiated,
            ..
        } = self;
        if dependent.is_empty() {
            return false;
        }
        family_calls(view).into_iter().any(|call| {
            *instantiated.entry(call).or_insert_with_key(|call| {
                let body = call.replace_multiple(dependent.as_slice());
                body != *call && class.names_energy(body.as_view())
            })
        })
    }
}

/// The summands of one residue root, by whether the class has already
/// projected them.
struct RootSummands<'a> {
    root: AtomView<'a>,
    /// The number of summands of `root`, none when it is zero.
    count: usize,
    /// Projected summands `c OSE(rep)^-a` by their order `a`. Order zero
    /// holds the summands without any class energy.
    projected: BTreeMap<i64, Vec<AtomView<'a>>>,
    /// Summands that still need their Laurent expansion.
    pending: Vec<AtomView<'a>>,
}

impl<'a> RootSummands<'a> {
    fn new(root: &'a Atom, orders: &mut ProjectedOrders<'_>) -> Self {
        let summands = match root.as_view() {
            _ if root.is_zero() => Vec::new(),
            AtomView::Add(sum) => sum.iter().collect(),
            view => vec![view],
        };
        let mut split = Self {
            root: root.as_view(),
            count: summands.len(),
            projected: BTreeMap::new(),
            pending: Vec::new(),
        };
        for summand in summands {
            match orders.order(summand) {
                Some(order) => split.projected.entry(order).or_default().push(summand),
                None => split.pending.push(summand),
            }
        }
        split
    }

    /// Whether a summand still shows a class energy: it needs its expansion
    /// or is projected at a nonzero order.
    fn shows_class(&self) -> bool {
        !self.pending.is_empty() || self.projected.keys().any(|order| *order != 0)
    }

    /// The summands that still need their expansion, as one root.
    fn pending_root(&self) -> Atom {
        if self.pending.len() == self.count {
            self.root.to_owned()
        } else {
            Atom::add_many(self.pending.iter().copied())
        }
    }

    /// The projected summands of `order` with their class propagators
    /// stripped off: their coefficient `[x^-order]`.
    fn projected_root(&self, order: i64, orders: &ProjectedOrders<'_>) -> Atom {
        self.projected.get(&order).map_or(Atom::Zero, |summands| {
            Atom::add_many(summands.iter().map(|summand| orders.coefficient(*summand)))
        })
    }
}

/// One residue row with the summands of each of its roots split.
struct SplitRow<'a> {
    key: &'a DirectResidueKey,
    roots: &'a Integrands,
    cuts: Vec<(CutCFFIndex, RootSummands<'a>)>,
}

impl<'a> SplitRow<'a> {
    fn new(
        key: &'a DirectResidueKey,
        roots: &'a Integrands,
        orders: &mut ProjectedOrders<'_>,
    ) -> Self {
        Self {
            key,
            roots,
            cuts: roots
                .iter()
                .map(|(cut, root)| (*cut, RootSummands::new(root, orders)))
                .collect(),
        }
    }

    /// The family of `rows` whose roots are `root` of their summands, with
    /// the shared definitions of each row.
    fn family(
        rows: &[Self],
        root: impl Fn(&RootSummands<'a>) -> Atom,
    ) -> Result<DirectResidueBranches> {
        DirectResidueBranches::from_keyed(
            rows.iter()
                .map(|row| {
                    Ok((
                        row.key.clone(),
                        row.cuts
                            .iter()
                            .map(|(cut, split)| (*cut, root(split)))
                            .collect::<Integrands>()
                            .with_numerators(row.roots.numerators().iter().cloned())?
                            .with_scalar_definitions(
                                row.roots.scalar_definitions().iter().cloned(),
                            )?,
                    ))
                })
                .collect::<Result<Vec<_>>>()?,
        )
    }
}

/// The distinct numerator-family calls in `view`.
fn family_calls(view: AtomView<'_>) -> BTreeSet<Atom> {
    let family = Integrands::family_symbol();
    let mut calls = BTreeSet::new();
    view.visitor(&mut |inner| match inner {
        AtomView::Fun(function) if function.get_symbol() == family => {
            calls.insert(inner.to_owned());
            false
        }
        _ => true,
    });
    calls
}

/// The first residue root of `branches`, or numerator body at the arguments
/// of a call in one, for which `holds` is true. Each shared body is tested
/// once; only a body for which it holds is instantiated at the calls, since
/// the arguments of a row may select terms of the body that vanish on it.
/// Bodies that only other summands call are never instantiated.
fn find_expression(
    branches: &DirectResidueBranches,
    mut holds: impl FnMut(AtomView<'_>) -> Result<bool, ZeroModeError>,
) -> Result<Option<Atom>, ZeroModeError> {
    let mut candidates = Vec::new();
    for entry in branches.numerator_definitions().into_values() {
        if holds(entry.rhs.as_view())? {
            candidates.push(entry.replacement());
        }
    }
    let mut calls = BTreeSet::new();
    for (_, roots) in branches.iter_keys() {
        for (_, root) in roots.iter() {
            if holds(root.as_view())? {
                return Ok(Some(root.clone()));
            }
            if !candidates.is_empty() {
                calls.append(&mut family_calls(root.as_view()));
            }
        }
    }
    for call in calls {
        let body = call.replace_multiple(&candidates);
        if body != call && holds(body.as_view())? {
            return Ok(Some(body));
        }
    }
    Ok(None)
}

/// A branch family split by whether one class has already projected its
/// summands.
struct ProjectedSplit<'a> {
    /// The coefficients `[x^-a]` of the projected summands, by order `a`.
    projected: BTreeMap<u32, DirectResidueBranches>,
    /// The summands that still need their Laurent expansion, if any.
    pending: Option<Cow<'a, DirectResidueBranches>>,
}

impl ZeroModeClass {
    /// Validate `edges` of `graph` as one frequency class. Each order `a`
    /// requests the coefficient of `x^-a`; the zero mode itself is
    /// `a = 2 nu`. The smallest edge is the class representative.
    ///
    /// # Errors
    ///
    /// The class needs an edge, and its edges must be bosons without a
    /// chemical potential that carry one frequency at one on-shell energy.
    /// Every order must lie in `1..=2 nu`, between the zero mode and the
    /// regular part.
    pub(super) fn new(
        graph: &Graph,
        edges: impl IntoIterator<Item = EdgeIndex>,
        orders: impl IntoIterator<Item = u32>,
    ) -> Result<Self, ZeroModeError> {
        let spec =
            ZeroModeClassSpec::new(edges.into_iter().collect()).ok_or(ZeroModeError::EmptyClass)?;
        let orders = orders.into_iter().collect::<BTreeSet<_>>();
        let representative = spec.representative();
        let n_edges = graph.underlying.n_edges();
        for &edge in spec.edges() {
            if usize::from(edge) >= n_edges {
                return Err(ZeroModeError::UnknownEdge { edge });
            }
            let data = &graph[edge];
            data.zero_mode()
                .map_err(|reason| ZeroModeError::NoZeroMode { edge, reason })?;
            if !graph
                .loop_momentum_basis
                .edges_are_raised(representative, edge)
                || data.mass_atom() != graph[representative].mass_atom()
            {
                return Err(ZeroModeError::NotRaised {
                    first: representative,
                    second: edge,
                });
            }
        }
        let max = 2 * u64::from(spec.power());
        if orders.is_empty()
            || orders
                .iter()
                .any(|order| !(1..=max).contains(&u64::from(*order)))
        {
            return Err(ZeroModeError::InvalidOrders {
                edges: spec.edges().clone(),
                orders,
                max,
            });
        }
        Ok(Self { spec, orders })
    }

    fn edges(&self) -> &BTreeSet<EdgeIndex> {
        self.spec.edges()
    }

    fn representative(&self) -> EdgeIndex {
        self.spec.representative()
    }

    /// The smallest requested order sets the series endpoint `x^-a_min`.
    fn lowest_order(&self) -> i64 {
        i64::from(
            *self
                .orders
                .first()
                .expect("a validated zero-mode class requests an order"),
        )
    }

    fn contains(&self, call: &ThermalDistributionCall) -> bool {
        self.edges().contains(&call.edge)
    }

    /// The static class propagators `OSE(rep)^-a` of the order `a`.
    fn restored(&self, order: u32) -> Atom {
        ose_atom_from_index(self.representative()).pow(Atom::num(-i64::from(order)))
    }

    /// Split `branches` into the summands this class has already projected
    /// and the rest. A summand `c OSE(rep)^-a` whose `c` names no class
    /// energy, also not through the numerator families it calls at their
    /// arguments, contributes `c` to `[x^-a]` and nothing to any other
    /// coefficient, without inlining, resolving or a series.
    ///
    /// # Errors
    ///
    /// A summand without a visible class energy that hides one in a routed
    /// on-shell energy is refused while another summand of the family still
    /// shows the class energy.
    fn split_projected<'a>(
        &self,
        branches: &'a DirectResidueBranches,
    ) -> Result<ProjectedSplit<'a>, ZeroModeError> {
        let mut orders = ProjectedOrders::new(self, branches);
        let rows = branches
            .iter_keys()
            .map(|(key, roots)| SplitRow::new(key, roots, &mut orders))
            .collect::<Vec<_>>();
        let splits = || {
            rows.iter()
                .flat_map(|row| row.cuts.iter().map(|(_, split)| split))
        };
        if splits().all(|split| split.projected.is_empty()) {
            return Ok(ProjectedSplit {
                projected: BTreeMap::new(),
                pending: Some(Cow::Borrowed(branches)),
            });
        }
        // Zero is exact for a summand without a class energy, also where an
        // earlier operation rewrote every class energy into a routed, massive
        // one. Beside a summand that still shows the class energy, a hidden
        // one may cancel a pole of it only in their sum.
        if splits().any(RootSummands::shows_class) {
            for summand in splits()
                .filter_map(|split| split.projected.get(&0))
                .flatten()
            {
                self.check_routed_energies(*summand)?;
            }
        }
        let pending = if splits().any(|split| !split.pending.is_empty()) {
            Some(Cow::Owned(SplitRow::family(
                &rows,
                RootSummands::pending_root,
            )?))
        } else {
            None
        };
        let mut projected = BTreeMap::new();
        for order in &self.orders {
            let exponent = i64::from(*order);
            if splits().all(|split| !split.projected.contains_key(&exponent)) {
                continue;
            }
            // Strip the class propagators; they are restored once every
            // class is projected.
            projected.insert(
                *order,
                SplitRow::family(&rows, |split| split.projected_root(exponent, &orders))?,
            );
        }
        Ok(ProjectedSplit { projected, pending })
    }

    /// The coefficients `[x^-a]` of `branches`, keyed by `a`, by a Laurent
    /// expansion in the class energy.
    fn expand(
        &self,
        graph: &Graph,
        branches: &DirectResidueBranches,
    ) -> Result<BTreeMap<u32, DirectResidueBranches>, ZeroModeError> {
        // Templates of earlier operations would carry class energies in
        // formal arguments, where neither substitution nor series may act.
        let branches = branches
            .inline_factors(|entry| self.names_energy(entry.rhs.as_view()))?
            .resolved_scalars()?;
        let x = Atom::var(GS.zero_mode_energy);
        let energies = self
            .edges()
            .iter()
            .map(|edge| Replacement::new(ose_atom_from_index(*edge).to_pattern(), x.to_pattern()))
            .collect::<Vec<_>>();
        // A truncated distribution N^(d) errs from x^(r-d) on. Products may
        // be more singular than the zero mode, e.g. through E0 + E1 = 2x of a
        // raised class, as long as their sum is not: keep enough orders that
        // the deepest product never reaches x^-a_min.
        let mut deepest = 0;
        let deformed = branches
            .map_numerators(|entry| {
                let body = entry.rhs.replace_multiple(&energies);
                self.check_numerator(body.as_view())?;
                Ok(body)
            })?
            .fallible_map(|_, root| {
                let root = self.expose_distributions(root)?.replace_multiple(&energies);
                self.check_routed_energies(root.as_view())?;
                deepest = deepest.min(self.valuation(root.as_view())?);
                Ok(root)
            })?;
        let relative_order = u32::try_from((-self.lowest_order() - deepest).max(0))?;
        let deformed = deformed.fallible_map(|_, root| {
            Ok(self.insert_laurent_series(graph, root, relative_order)?)
        })?;
        let series = deformed.series_preserving_numerators(
            GS.zero_mode_energy,
            Atom::Zero.as_view(),
            -self.lowest_order(),
            DirectResidueBranches::numerator_scope().1,
        )?;
        for (_, roots) in series.iter_keys() {
            if let Some(entry) = roots
                .numerators()
                .iter()
                .chain(roots.scalar_definitions())
                .find(|entry| entry.rhs.contains_symbol(GS.zero_mode_energy))
            {
                return Err(ZeroModeError::NotLaurent {
                    edges: self.edges().clone(),
                    expression: entry.rhs.clone(),
                });
            }
        }
        self.orders
            .iter()
            .map(|order| Ok((*order, self.coefficient(&series, *order)?)))
            .collect()
    }

    /// Validate where x enters a deformed root and bound the x-valuation of
    /// each of its products from below. x may enter polynomially, as in
    /// numerator and contact surfaces, and through energy denominators,
    /// on-shell normalizations and the class distributions, the latter also
    /// in positive powers. Any other function would hide x from the Laurent
    /// series. A class distribution `N^(d)` counts as its pole `x^-(d+1)`, x
    /// itself as `x^1`, and a linear energy sum, retained as a surface or
    /// not, as `x^1` when it vanishes at `x = 0`. A retained surface also
    /// enters polynomially where a scalar-pole cancellation brought surfaces
    /// over a common denominator.
    fn valuation(&self, view: AtomView<'_>) -> Result<i64, ZeroModeError> {
        if !view.contains_symbol(GS.zero_mode_energy) {
            return Ok(0);
        }
        match view {
            AtomView::Add(sum) => sum
                .iter()
                .try_fold(0, |deepest, term| Ok(deepest.min(self.valuation(term)?))),
            AtomView::Mul(product) => product
                .iter()
                .try_fold(0, |total, factor| Ok(total + self.valuation(factor)?)),
            // The only variable that contains x is x itself.
            AtomView::Var(_) => Ok(1),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                match i64::try_from(exponent) {
                    Ok(exponent) if exponent < 0 => self
                        .surface_valuation(base)
                        .map(|vanishes| exponent * vanishes)
                        .ok_or_else(|| self.formal_argument(view)),
                    Ok(exponent) => Ok(exponent * self.valuation(base)?),
                    Err(_) => Err(self.formal_argument(view)),
                }
            }
            AtomView::Fun(function) if function.get_symbol() == GS.thermal_distribution => {
                let call = ThermalDistributionCall::try_from(view)?;
                self.check_class_call(&call)?;
                Ok(-1 - i64::try_from(call.derivative_order)?)
            }
            AtomView::Fun(function) if function.get_symbol() == GS.esurface => self
                .surface_valuation(view)
                .ok_or_else(|| self.formal_argument(view)),
            _ => Err(self.formal_argument(view)),
        }
    }

    /// A linear energy sum: x, a retained surface, or a sum whose terms are
    /// each x-free or a number times x. Its valuation is one when it vanishes
    /// at `x = 0` and zero otherwise; any other dependence on x is not an
    /// energy sum.
    fn surface_valuation(&self, view: AtomView<'_>) -> Option<i64> {
        let x = Atom::var(GS.zero_mode_energy);
        // Whether a valid term is a number times x rather than x-free.
        let dependent = |term: AtomView<'_>| match term {
            _ if !term.contains_symbol(GS.zero_mode_energy) => Some(false),
            AtomView::Mul(product) => {
                let mut factors = product
                    .iter()
                    .filter(|factor| factor.contains_symbol(GS.zero_mode_energy));
                (factors.next().is_some_and(|factor| factor == x.as_view())
                    && factors.next().is_none())
                .then_some(true)
            }
            _ => (term == x.as_view()).then_some(true),
        };
        let terms = match view {
            AtomView::Var(_) => vec![view],
            AtomView::Add(sum) => sum.iter().collect(),
            AtomView::Fun(function) if function.get_symbol() == GS.esurface => {
                function.iter().collect()
            }
            _ => return None,
        };
        terms.into_iter().try_fold(1, |vanishes, term| {
            Some(vanishes.min(i64::from(dependent(term)?)))
        })
    }

    fn formal_argument(&self, view: AtomView<'_>) -> ZeroModeError {
        ZeroModeError::FormalArgument {
            edges: self.edges().clone(),
            expression: view.to_owned(),
        }
    }

    /// A class distribution must be taken at the class energy x and at
    /// finite temperature: any other energy or limit would hide x from the
    /// substitution or misstate its pole.
    fn check_class_call(&self, call: &ThermalDistributionCall) -> Result<(), ZeroModeError> {
        if self.contains(call)
            && call.energy == Atom::var(GS.zero_mode_energy)
            && call.temperature_flag.is_one()
        {
            Ok(())
        } else {
            Err(ZeroModeError::DeformedDistribution {
                edges: self.edges().clone(),
                edge: call.edge,
                energy: call.energy.clone(),
                temperature_flag: call.temperature_flag.clone(),
            })
        }
    }

    /// Whether `view` names a class edge, through x or through a call of one
    /// of `heads` whose first argument is a class edge.
    fn names_class(&self, view: AtomView<'_>, heads: &[Symbol]) -> bool {
        if view.contains_symbol(GS.zero_mode_energy) {
            return true;
        }
        let mut found = false;
        view.visitor(&mut |inner| {
            if found {
                return false;
            }
            if let AtomView::Fun(function) = inner
                && heads.contains(&function.get_symbol())
                && let Some(edge) = function.iter().next()
                && usize::try_from(edge).is_ok_and(|edge| self.edges().contains(&EdgeIndex(edge)))
            {
                found = true;
            }
            !found
        });
        found
    }

    /// Whether `view` depends on a class energy: x, a class `OSE` or a class
    /// distribution.
    fn names_energy(&self, view: AtomView<'_>) -> bool {
        self.names_class(view, &[GS.ose, GS.thermal_distribution])
    }

    /// A routed `OnShellEnergy` that names a class edge would keep a class
    /// energy away from the substitution.
    fn check_routed_energies(&self, view: AtomView<'_>) -> Result<(), ZeroModeError> {
        let edge_heads = [GS.ose, GS.energy, GS.emr_mom, GS.emr_vec];
        let mut routed = None;
        view.visitor(&mut |inner| {
            if routed.is_some() {
                return false;
            }
            match inner {
                AtomView::Fun(function) if function.get_symbol() == GS.on_shell_energy => {
                    if self.names_class(inner, &edge_heads) {
                        routed = Some(inner.to_owned());
                    }
                    false
                }
                _ => true,
            }
        });
        routed.map_or(Ok(()), |expression| {
            Err(ZeroModeError::RoutedEnergy {
                edges: self.edges().clone(),
                expression,
            })
        })
    }

    /// The Laurent terms `(k, [x^k])` of an expanded root.
    fn laurent_terms(&self, root: &Atom) -> Result<Vec<(i64, Atom)>, ZeroModeError> {
        let x = Atom::var(GS.zero_mode_energy);
        let not_laurent = || ZeroModeError::NotLaurent {
            edges: self.edges().clone(),
            expression: root.clone(),
        };
        root.coefficient_list_exact(slice::from_ref(&x))
            .ok_or_else(not_laurent)?
            .into_iter()
            .map(|(monomial, coefficient)| {
                let exponent = match monomial.as_view() {
                    AtomView::Num(_) if monomial.is_one() => 0,
                    AtomView::Var(_) if monomial == x => 1,
                    AtomView::Pow(power) if power.get_base_exp().0 == x.as_view() => {
                        i64::try_from(power.get_base_exp().1).map_err(|_| not_laurent())?
                    }
                    _ => return Err(not_laurent()),
                };
                Ok((exponent, coefficient))
            })
            .collect()
    }

    /// Move the class distributions out of each medium-weight wrapper. The
    /// wrapper is numerically the identity; it keeps the remaining x-free
    /// occupations together for degree certification and vacuum projection.
    fn expose_distributions(&self, root: &Atom) -> Result<Atom, ZeroModeError> {
        let mut error = None;
        let exposed = root.replace_map(|view, _, out| {
            if error.is_some() {
                return;
            }
            let AtomView::Fun(wrapper) = view else {
                return;
            };
            if wrapper.get_symbol() != GS.thermal_weight_wrapper || wrapper.get_nargs() != 1 {
                return;
            }
            let weight = wrapper.get(0);
            let mut calls = BTreeSet::new();
            weight.visitor(&mut |inner| {
                if error.is_some() {
                    return false;
                }
                let AtomView::Fun(function) = inner else {
                    return true;
                };
                if function.get_symbol() != GS.thermal_distribution {
                    return true;
                }
                match ThermalDistributionCall::try_from(inner) {
                    Ok(call) if self.contains(&call) => {
                        calls.insert(inner.to_owned());
                    }
                    Ok(_) => {}
                    Err(parse_error) => error = Some(parse_error.into()),
                }
                false
            });
            if calls.is_empty() {
                return;
            }
            let calls = calls.into_iter().collect::<Vec<_>>();
            match weight.coefficient_list_exact(&calls) {
                Some(terms) => {
                    **out = terms
                        .into_iter()
                        .map(|(monomial, coefficient)| {
                            function!(GS.thermal_weight_wrapper, coefficient) * monomial
                        })
                        .sum();
                }
                // The class distributions are not polynomial in the weight.
                None => error = Some(self.formal_argument(view)),
            }
        });
        error.map_or(Ok(exposed), Err)
    }

    /// Replace the class distributions by their Laurent polynomials through
    /// `relative_order` orders past their poles, at the class's own zero-mode
    /// temperature. Each derivative order keeps an opaque remainder at its
    /// first omitted power, which certifies that no truncated term reaches a
    /// requested coefficient.
    fn insert_laurent_series(
        &self,
        graph: &Graph,
        root: &Atom,
        relative_order: u32,
    ) -> Result<Atom, ZeroModeError> {
        let x = Atom::var(GS.zero_mode_energy);
        let inverse_temperature = class_zero_mode_temperature(self.representative())
            .pow(-1)
            .to_pattern();
        ThermalDistributionCall::replace_selected(
            root.as_view(),
            |call| {
                if !self.contains(call) {
                    return Ok(false);
                }
                // A class call with an x-free energy escaped the substitution.
                self.check_class_call(call)?;
                Ok(true)
            },
            |call| {
                let derivative_order = i64::try_from(call.derivative_order)?;
                let laurent = graph
                    .zero_mode_distribution_laurent(
                        call.edge,
                        call.derivative_order,
                        call.thermal_sign,
                        GS.zero_mode_energy,
                        usize::try_from(relative_order)?,
                    )?
                    .replace(GS.inverse_temperature)
                    .with(&inverse_temperature);
                let omitted = i64::from(relative_order) - derivative_order;
                let remainder = function!(GS.zero_mode_remainder, derivative_order);
                Ok(laurent + x.pow(Atom::num(omitted)) * remainder)
            },
        )
    }

    /// Shared bodies may depend polynomially on class energies, but never
    /// through a distribution, a medium weight, a routed energy or a
    /// template.
    fn check_numerator(&self, body: AtomView<'_>) -> Result<(), ZeroModeError> {
        let formal = [
            GS.thermal_weight_wrapper,
            Integrands::scalar_symbol(),
            Integrands::factor_symbol(),
            Integrands::family_symbol(),
        ];
        let mut invalid = None;
        body.visitor(&mut |view| {
            if invalid.is_some() {
                return false;
            }
            let dependent = view.contains_symbol(GS.zero_mode_energy);
            match view {
                AtomView::Fun(function) if function.get_symbol() == GS.thermal_distribution => {
                    invalid = match ThermalDistributionCall::try_from(view) {
                        Ok(call) if !dependent && !self.contains(&call) => None,
                        Ok(_) => Some(self.formal_argument(view)),
                        Err(parse_error) => Some(parse_error.into()),
                    };
                    false
                }
                AtomView::Fun(function) if function.get_symbol() == GS.on_shell_energy => {
                    invalid = self.check_routed_energies(view).err();
                    false
                }
                AtomView::Fun(function) if dependent && formal.contains(&function.get_symbol()) => {
                    invalid = Some(self.formal_argument(view));
                    false
                }
                AtomView::Pow(power)
                    if dependent
                        && !i64::try_from(power.get_base_exp().1)
                            .is_ok_and(|exponent| exponent >= 0) =>
                {
                    invalid = Some(ZeroModeError::NotPolynomial {
                        edges: self.edges().clone(),
                        expression: view.to_owned(),
                    });
                    false
                }
                _ => true,
            }
        });
        invalid.map_or(Ok(()), Err)
    }

    /// `[x^-a]` of an expanded family.
    fn coefficient(
        &self,
        series: &DirectResidueBranches,
        order: u32,
    ) -> Result<DirectResidueBranches, ZeroModeError> {
        let exponent = -i64::from(order);
        Ok(series.fallible_map(|_, root| {
            let coefficient = self
                .laurent_terms(root)?
                .into_iter()
                .find_map(|(power, coefficient)| (power == exponent).then_some(coefficient))
                .unwrap_or(Atom::Zero);
            if coefficient.contains_symbol(GS.zero_mode_remainder) {
                return Err(ZeroModeError::TruncationReached {
                    edges: self.edges().clone(),
                    order,
                }
                .into());
            }
            if coefficient.contains_symbol(GS.zero_mode_energy) {
                return Err(ZeroModeError::NotLaurent {
                    edges: self.edges().clone(),
                    expression: coefficient,
                }
                .into());
            }
            Ok(coefficient)
        })?)
    }
}

/// Iterated zero-mode projection onto disjoint frequency classes.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(super) struct ZeroModeProjection {
    classes: Vec<ZeroModeClass>,
    /// The edges of every class.
    edges: BTreeSet<EdgeIndex>,
    sector: ZeroModeSector,
}

impl ZeroModeProjection {
    /// Project onto `classes` in turn, reading the coefficients from `sector`.
    ///
    /// # Errors
    ///
    /// The projection needs at least one class, its classes must be disjoint,
    /// and the hard-vacuum sector only serves the order `2 nu - 1` of each
    /// class.
    pub(super) fn new(
        classes: Vec<ZeroModeClass>,
        sector: ZeroModeSector,
    ) -> Result<Self, ZeroModeError> {
        if classes.is_empty() {
            return Err(ZeroModeError::NoClass);
        }
        let mut edges = BTreeSet::new();
        for class in &classes {
            if let Some(edge) = class.edges().iter().find(|edge| !edges.insert(**edge)) {
                return Err(ZeroModeError::OverlappingClasses { edge: *edge });
            }
            let expected = 2 * u64::from(class.spec.power()) - 1;
            if sector == ZeroModeSector::HardVacuum
                && !class
                    .orders
                    .iter()
                    .all(|order| u64::from(*order) == expected)
            {
                return Err(ZeroModeError::HardVacuumOrders {
                    edges: class.edges().clone(),
                    orders: class.orders.clone(),
                    expected,
                });
            }
        }
        Ok(Self {
            classes,
            edges,
            sector,
        })
    }

    /// The coefficients `prod_i OSE(rep_i)^-a_i [prod_i x_i^-a_i]` of every
    /// requested order tuple: [`Self::stripped_coefficients`] times
    /// [`Self::class_propagators`], which are top-level factors of every
    /// residue root.
    ///
    /// # Errors
    ///
    /// As for [`Self::stripped_coefficients`].
    #[cfg(test)]
    pub(super) fn coefficients(
        &self,
        graph: &Graph,
        branches: &DirectResidueBranches,
    ) -> Result<BTreeMap<Vec<u32>, DirectResidueBranches>, ZeroModeError> {
        Ok(self
            .stripped_coefficients(graph, branches)?
            .into_iter()
            .map(|(orders, coefficient)| {
                let propagators = self.class_propagators(&orders);
                (orders, coefficient.map(|root| root * &propagators))
            })
            .collect())
    }

    /// The Laurent coefficients `[prod_i x_i^-a_i]` of every requested order
    /// tuple, keyed by `(a_1, a_2, ...)` in class order, without their class
    /// propagators [`Self::class_propagators`]. A caller that transforms them
    /// further multiplies the propagators in only afterwards, so that no
    /// Taylor series wraps them into a scalar coefficient, where a later
    /// projection could not read them off.
    ///
    /// # Errors
    ///
    /// Fails for class energies in a formal function or a nonlinear
    /// denominator, for class distributions away from the class energy, for
    /// routed class energies, and for bodies that are not polynomial in a
    /// class energy. A frozen soft branch that reaches the hard-vacuum sector
    /// unprojected, and a projected summand of the hard-vacuum order that
    /// keeps an occupation outside the classes, are internal errors.
    pub(super) fn stripped_coefficients(
        &self,
        graph: &Graph,
        branches: &DirectResidueBranches,
    ) -> Result<BTreeMap<Vec<u32>, DirectResidueBranches>, ZeroModeError> {
        // Every completed soft branch is frozen relative to the whole graph.
        let whole = graph.full_filter();
        let soft_branches = SoftBranches::new(graph, &whole, Local3DDeformation::Ordinary);
        let mut projected = BTreeMap::from([(Vec::new(), Cow::Borrowed(branches))]);
        for class in &self.classes {
            let mut next = BTreeMap::new();
            for (prefix, branches) in projected {
                for (order, coefficient) in
                    self.class_coefficients(graph, &soft_branches, class, &branches)?
                {
                    let mut orders = prefix.clone();
                    orders.push(order);
                    next.insert(orders, Cow::Owned(coefficient));
                }
            }
            projected = next;
        }
        // A later projection reads restored propagators off without a new
        // expansion only while no series has wrapped them into a scalar
        // coefficient, so they are restored after the last class.
        Ok(projected
            .into_iter()
            .map(|(orders, coefficient)| (orders, coefficient.into_owned()))
            .collect())
    }

    /// The static class propagators `prod_i OSE(rep_i)^-a_i` of the order
    /// tuple `orders`, given in class order.
    pub(super) fn class_propagators(&self, orders: &[u32]) -> Atom {
        self.classes
            .iter()
            .zip(orders)
            .fold(Atom::one(), |restored, (class, order)| {
                restored * class.restored(*order)
            })
    }

    /// The coefficients `[x^-a]` of one class, keyed by `a`: summands the
    /// class has already projected are read off, the others are expanded in
    /// this projection's sector.
    fn class_coefficients(
        &self,
        graph: &Graph,
        soft_branches: &SoftBranches<'_>,
        class: &ZeroModeClass,
        branches: &DirectResidueBranches,
    ) -> Result<Vec<(u32, DirectResidueBranches)>, ZeroModeError> {
        let ProjectedSplit {
            mut projected,
            pending,
        } = class.split_projected(branches)?;
        if self.sector == ZeroModeSector::HardVacuum {
            // Projected summands are read off without the sector's vacuum
            // projection, so they must already be hard-vacuum.
            for family in projected.values() {
                if let Some(expression) =
                    find_expression(family, |view| self.keeps_hard_occupation(view))?
                {
                    return Err(ZeroModeError::HardVacuumOccupation {
                        edges: class.edges().clone(),
                        expression,
                    });
                }
            }
        }
        let mut expanded = match pending {
            Some(pending) => {
                let pending = self.read_in_sector(graph, soft_branches, class, pending)?;
                class.expand(graph, &pending)?
            }
            None => BTreeMap::new(),
        };
        class
            .orders
            .iter()
            .map(|order| {
                let coefficient = match (expanded.remove(order), projected.remove(order)) {
                    (Some(expanded), Some(projected)) => expanded.zip_add(&projected)?,
                    (Some(coefficient), None) | (None, Some(coefficient)) => coefficient,
                    (None, None) => branches.fallible_map(|_, _| Ok(Atom::Zero))?,
                };
                Ok((*order, coefficient))
            })
            .collect()
    }

    /// The summands `class` still has to expand, as this projection's sector
    /// reads them. In the hard-vacuum sector every occupation outside the
    /// classes takes its vacuum value, which a frozen soft branch must keep.
    fn read_in_sector<'a>(
        &self,
        graph: &Graph,
        soft_branches: &SoftBranches<'_>,
        class: &ZeroModeClass,
        pending: Cow<'a, DirectResidueBranches>,
    ) -> Result<Cow<'a, DirectResidueBranches>, ZeroModeError> {
        if self.sector == ZeroModeSector::Complete {
            return Ok(pending);
        }
        if let Some(expression) =
            find_expression(&pending, |view| Ok(soft_branches.contains_frozen(view)?))?
        {
            return Err(ZeroModeError::FrozenHardVacuum {
                edges: class.edges().clone(),
                expression,
            });
        }
        Ok(Cow::Owned(pending.map_expressions(|atom| {
            graph.make_selected_thermal_distributions_explicit(atom, MediumMode::Vacuum, |call| {
                Ok(!self.edges.contains(&call.edge))
            })
        })?))
    }

    /// Whether `view` keeps the occupation of a line outside every class,
    /// which the hard-vacuum sector takes at its vacuum value.
    fn keeps_hard_occupation(&self, view: AtomView<'_>) -> Result<bool, ZeroModeError> {
        let mut found = None;
        view.visitor(&mut |inner| {
            if found.is_some() {
                return false;
            }
            match inner {
                AtomView::Fun(function) if function.get_symbol() == GS.thermal_distribution => {
                    found = match ThermalDistributionCall::try_from(inner) {
                        Ok(call) if self.edges.contains(&call.edge) => None,
                        Ok(_) => Some(Ok(())),
                        Err(parse_error) => Some(Err(parse_error)),
                    };
                    false
                }
                _ => true,
            }
        });
        Ok(found.transpose()?.is_some())
    }
}
