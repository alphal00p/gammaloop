//! IRH-DEPTH round 5: the quotient bound `B` ([`super::depth`]) is trusted
//! only on factors that are PROVABLY regular, or whose pole order is
//! provably read.
//!
//! The child route (`GAMMALOOP_IR_HARD=child`) truncates the hard part of
//! an IR correction at `ε^(extra-1)` before it knows the pole order of the
//! quotient factors it multiplies; `B` is that pole order's bound. `B` adds
//! the excess (`depth::excess_poles`) of every factor of the product to the
//! quotient's loops. An excess read from the expression as written is only
//! a proof when the expression's poles are visible to it: a pole can hide
//! in a closed index trace (`-2/(g_ab g^ab - 4)` is `1/ε`, IRH-DEPTH review
//! of 2026-11-21, finding 1), in an opaque function or in a symbol whose
//! value comes later, or in a factor the bound never reads (the graph's
//! `numerator_prefactor`, finding 2).
//!
//! # The rule
//!
//! Every factor that enters the product is named in [`FACTORS`] (a closed
//! list, with the place where each is multiplied). Before the hard part is
//! integrated, each runtime factor is NORMALISED ([`normalise`]): every
//! tensor structure that sits in a denominator (or in a function argument)
//! must be a closed trace (no open index; fresh dummy indices included),
//! and is contracted, with Lorentz slots in `d` dimensions (conventional
//! dimensional regularisation, as the integration itself); then `d = 4-2ε`
//! is substituted in every scalar. Its excess is read on the normalised
//! form. A factor is REFUSED, and with it the WHOLE correction goes to
//! today's (generic) path before the hard truncation is used, loudly
//! (`IR_HARD_DEPTH_REFUSED reason=<class>`), counted per class in
//! `IR_HARD_DEPTH_SUMMARY`, when after normalisation it still holds
//! * an unresolved tensor structure in a denominator or a function argument
//!   (an open index, or a contraction the algebra does not resolve to a
//!   scalar): `unresolved_tensor`;
//! * a function that is neither a tensor structure (spenso tensor,
//!   representation or index symbol) nor a known constant (`zeta`, `log`,
//!   `exp`, `sqrt`, `sin`, `cos`, `abs`, `conj` of pure numbers), or a power
//!   with a symbolic exponent or a non-integer power of an `ε`-dependent
//!   base: `opaque_function`;
//! * a symbol other than `ε` and `d` in a denominator (its value is unknown,
//!   it could be `ε`): `symbol_in_denominator`;
//! * a factor this list does not name: `outside_list` (raised when a caller
//!   asks to certify a factor kind the child route does not know).
//!
//! The hard numerator ([`Factor::HardNumerator`]) is stricter: `H(γ)` is
//! integrated and TRUNCATED at `ε^(extra-1)` with `d` still symbolic in its
//! scalar coefficients, so a `1/(d-4)` (or a closed trace in a denominator,
//! which contracts to one) is a pole of `h_j` the truncation cannot see
//! (+117/32 against −117/64 on L3_006 with v0·ε⁴·(−2/(d−4))). There a
//! `d`-denominator is admitted only when it does not vanish at `d = 4`
//! (`1/(d-1)`, `1/(d+2)`), and a tensor in a denominator is refused without
//! contraction.
//!
//! A symbol in a numerator is a polynomial factor: it cannot lower a
//! valuation, so it needs no refusal (substituting a pole for it later is
//! outside the API: every listed place multiplies the factor as written).
//!
//! # Moment insertions and the tensor reduction
//!
//! The review's release ruling asks that every POSSIBLE moment insertion
//! `E_j` (also one the moment extraction drops, which no row check sees)
//! belongs to a certified regular tensor class, and notes that an
//! intermediate rank-ten metric projector has a `d-4` denominator. The
//! certificate is in two parts.
//! * The class: `E_j` is a leaf of the hard numerator's tensor structure
//!   (`moments::hard_moments` collects complete tensor leaves; scalars go to
//!   `h_j`). The hard numerator ([`Factor::HardNumerator`]) is normalised
//!   like the quotient's: a tensor in a denominator or a function argument
//!   is a refusal, so every leaf is a monomial in momenta, metrics, gamma
//!   matrices and colour tensors with positive integer powers. The scalar
//!   coefficients of the hard numerator may carry poles (they belong to
//!   `h_j`, which is exact through `ε^(extra-1)` whatever its poles).
//! * The reduction: the quotient integrals of such a monomial insertion are
//!   vacuum (and massless propagator) tensor integrals with a POLYNOMIAL
//!   numerator. Their complete tensor coefficient is regular at `d = 4`
//!   whatever intermediate projector computes it: by Gaussian integration in
//!   Schwinger parameters, `∫ Π_i k_{a_i}^{μ_i} exp(-k·A·k)` is a sum over
//!   metric pairings of products of `(A^{-1})_{ab}`, covariances independent
//!   of `d`, times the scalar integral; the reduction is an exact identity in
//!   `d` (rational arithmetic), so the `1/(d-4)` of a rank-ten projector
//!   cancels in the complete coefficient (an evanescent combination of the
//!   metric basis meets its own zero). What remains are the scalar
//!   integrals' poles, one per loop, which is the `l` of each `K` in `B`.
//!   Contracted traces of such monomials are polynomials in `d`. So every
//!   possible insertion is in the regular class once the hard numerator is
//!   normalised; [`Factor::MomentReduction`] names this step and its witness
//!   checks a closed trace and the rank-two reduction coefficient.
//!
//! # Production
//!
//! Production vertex and edge numerators are polynomials in tensor
//! structures with rational coefficients (0 ε, 0 d in all 80,681 L5 DOTs);
//! propagator-insertion surrogates carry `ε` polynomials, `ζ` values and
//! rational functions of `d` regular at `d = 4`; the graph prefactor and
//! projector are `1`. None is refused, and the normalised forms have the
//! excess the expressions as written have, so `B` and every record are
//! unchanged.

use idenso::{
    CookMode, CookSettings, IndexTooling,
    color::ColorSimplifySettings,
    dirac::GammaSimplifySettings,
    tensor::{AlgebraContraction, AlgebraSettings, SymbolicTensor},
};
use spenso::{
    network::tags::SPENSO_TAG,
    structure::{
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName},
    },
};
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};

use crate::{utils::GS, uv::approx::integrated::Integrated};

/// Every factor that enters the child route's product `h_j · z_j · closure
/// · compensation · root projection` (the closed list of the module doc).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Factor {
    QuotientVertexNumerator,
    QuotientEdgeNumerator,
    HardNumerator,
    PropagatorPower,
    InsertionCoefficient,
    GraphNumeratorPrefactor,
    GraphProjector,
    OverallFactor,
    RootClosure,
    RootNormalisation,
    RootProjector,
    RowInsertion,
    Compensation,
    MomentReduction,
}

/// One entry of [`FACTORS`].
#[derive(Clone, Copy, Debug)]
pub struct Entry {
    pub factor: Factor,
    /// The log name.
    pub name: &'static str,
    /// Where the factor is multiplied.
    pub multiplied: &'static str,
    /// How it is certified (or refused).
    pub certified: &'static str,
    /// The test that witnesses a pole in it.
    pub witness: &'static str,
}

/// THE closed list: every factor of the product, where it is multiplied and
/// how it is certified. A factor not named here is refused
/// ([`Refusal::OutsideList`]).
pub const FACTORS: [Entry; 14] = [
    Entry {
        factor: Factor::QuotientVertexNumerator,
        name: "quotient_vertex_numerator",
        multiplied: "Graph::numerator(support, hard) (numerator/uninit.rs fill_in_reduced: graph[node].get_num()), the integrand of every forest step of Z̃_F(G/γ; E_j) (Evaluation::ir_forest, quotient)",
        certified: "normalise + excess at run time, per correction, before H(γ) is integrated (Evaluation::ir_hard_child); its excess is in B",
        witness: "depth_tests::round5::the_quotient_vertex_pole_is_read_in_every_spelling",
    },
    Entry {
        factor: Factor::QuotientEdgeNumerator,
        name: "quotient_edge_numerator",
        multiplied: "Graph::numerator(support, hard) (fill_in_reduced: e.data.num), as the vertex numerators",
        certified: "as quotient_vertex_numerator (the same product)",
        witness: "depth_tests::round5::the_quotient_edge_pole_is_read_in_every_spelling",
    },
    Entry {
        factor: Factor::HardNumerator,
        name: "hard_numerator",
        multiplied: "Graph::numerator(hard, ∅): the integrand of H(γ) (Evaluation::integrate_ir_hard_term); its tensor leaves are the insertions E_j",
        certified: "normalise at run time before H(γ) is integrated (grammar of every possible E_j); its explicit scalar poles stay in h_j, but h_j is truncated before d = 4-2ε is substituted in it, so a d-denominator that vanishes at d = 4 (symbol_in_denominator) or any tensor in a denominator (unresolved_tensor, not contracted) is refused",
        witness: "depth_tests::round5::the_hard_vertex_spellings_agree",
    },
    Entry {
        factor: Factor::PropagatorPower,
        name: "propagator_power",
        multiplied: "the denominators 1/p² of every edge (uv_graph.rs denominator)",
        certified: "by type: feynkit's propagator_power is an i32 and the runtime graph rejects any power but 1 (graph/parse/mod.rs from_feynkit_diagram)",
        witness: "depth_tests::round5::a_propagator_power_cannot_carry_a_pole",
    },
    Entry {
        factor: Factor::InsertionCoefficient,
        name: "insertion_coefficient",
        multiplied: "a propagator-insertion surrogate's coefficient vertex (insertions.rs with_coefficient), a vertex numerator of the quotient or of the hard part",
        certified: "as quotient_vertex_numerator / hard_numerator (it IS a vertex numerator); insertion weights are rejected with ε (PropagatorInsertionEvaluator::validate_weight)",
        witness: "depth::tests::a_surrogate_coefficient_has_no_excess",
    },
    Entry {
        factor: Factor::GraphNumeratorPrefactor,
        name: "numerator_prefactor",
        multiplied: "not by R*: the caller's root projector (beta src/integrate.rs project_physical_root_mode: tensor · σ(global_prefactor.num), and closed_root_scalar_normalization), i.e. onto h_j z_j after the outer Taylor step",
        certified: "normalise + excess at run time (Graph::global_prefactor.num); its excess is in B; beta refuses a non-regular one at its public boundary (no route carries finite orders through a callback applied after an inner pole projection)",
        witness: "depth_tests::round5::a_numerator_prefactor_pole_is_read_in_every_spelling",
    },
    Entry {
        factor: Factor::GraphProjector,
        name: "projector",
        multiplied: "not by R* nor by beta's R* path (beta only multiplies it in its direct integration, src/inputs.rs)",
        certified: "normalise + excess at run time (Graph::global_prefactor.projector), in B as the prefactor; beta refuses a non-regular one",
        witness: "depth_tests::round5::a_projector_pole_is_read",
    },
    Entry {
        factor: Factor::OverallFactor,
        name: "overall_factor",
        multiplied: "not by R*: the caller after the counterterm's pole projection (beta finish_counterterm weight)",
        certified: "excluded from B (R*'s value does not contain it); beta refuses a non-regular one at its public boundary (a pole there needs finite orders no route keeps)",
        witness: "depth_tests::round5::the_graph_weights_are_certified",
    },
    Entry {
        factor: Factor::RootClosure,
        name: "root_closure",
        multiplied: "z_j · closure (Evaluation::ir_outer_integrated); a metric or Dirac closure",
        certified: "normalise + excess at run time; its excess is in B",
        witness: "depth::tests::the_bound_adds_the_quotient_loops_and_every_excess",
    },
    Entry {
        factor: Factor::RootNormalisation,
        name: "root_normalisation",
        multiplied: "the caller's closed-root normalisation c = transverse·σ(prefactor)/(i q²)/(d−1) (beta closed_root_scalar_normalization) after Taylor",
        certified: "σ(prefactor) as numerator_prefactor; i q² is the external invariant (ε-free), d−1 is regular at d = 4; beta builds it from these alone",
        witness: "depth_tests::round5::the_fixed_normalisations_are_regular",
    },
    Entry {
        factor: Factor::RootProjector,
        name: "root_projector",
        multiplied: "the caller's root projector callback on the complete outer polynomial (Evaluation::ir_outer_integrated)",
        certified: "a callback is opaque: by contract (documented on RStarPlan::projected_counterterm) regular in ε; beta's multiplies σ(numerator_prefactor) (certified above, refused by beta when poleful) and the fixed normalisations i p², d−1, d_A",
        witness: "depth_tests::round5::a_numerator_prefactor_pole_is_read_in_every_spelling",
    },
    Entry {
        factor: Factor::RowInsertion,
        name: "row_insertion",
        multiplied: "E_j of each row, the decorated quotient's insertion (Evaluation::ir_forest)",
        certified: "excess per returned row (depth::excess_poles); a dropped row is bounded through hard_numerator and moment_reduction",
        witness: "depth_tests::a_row_without_a_known_bound_is_refused_after_the_first_depth",
    },
    Entry {
        factor: Factor::Compensation,
        name: "compensation",
        multiplied: "Σ_k (-1)^k binom(a-εL, k) of the fixed-mass outer Taylor step (top_mass::compensation_factor)",
        certified: "a polynomial in ε by construction",
        witness: "depth_tests::round5::the_fixed_normalisations_are_regular",
    },
    Entry {
        factor: Factor::MomentReduction,
        name: "moment_reduction",
        multiplied: "the vacuum tensor reduction of the quotient integrals with a monomial insertion (vakint, closed.rs)",
        certified: "the complete coefficient of a polynomial tensor integrand is regular at d = 4 (module doc); its premise is hard_numerator's normalisation",
        witness: "depth_tests::round5::the_fixed_normalisations_are_regular",
    },
];

/// The entry of `factor`.
pub fn entry(factor: Factor) -> &'static Entry {
    FACTORS
        .iter()
        .find(|entry| entry.factor == factor)
        .expect("every factor is listed")
}

/// Why a factor is refused (module doc).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Refusal {
    UnresolvedTensor,
    OpaqueFunction,
    SymbolInDenominator,
    OutsideList,
}

impl Refusal {
    /// The `reason=` label (and the `IR_HARD_DEPTH_SUMMARY` field).
    pub fn label(self) -> &'static str {
        match self {
            Self::UnresolvedTensor => "unresolved_tensor",
            Self::OpaqueFunction => "opaque_function",
            Self::SymbolInDenominator => "symbol_in_denominator",
            Self::OutsideList => "outside_list",
        }
    }
}

/// A refused factor: its class, the factor and the offending piece.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Refused {
    pub class: Refusal,
    pub factor: &'static str,
    pub detail: String,
}

impl std::fmt::Display for Refused {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "IRH-DEPTH regularity: {} refused ({}): {}",
            self.factor,
            self.class.label(),
            self.detail
        )
    }
}

type Normalised = std::result::Result<Atom, (Refusal, String)>;

/// A tensor structure: a spenso tensor, representation or index symbol.
fn structural(symbol: Symbol) -> bool {
    let tags = &SPENSO_TAG;
    symbol.has_tag(&tags.tensor)
        || symbol.has_tag(&tags.representation)
        || symbol.has_tag(&tags.index)
        || symbol.has_tag(&tags.rank1)
        || symbol.has_tag(&tags.upper)
        || symbol.has_tag(&tags.lower)
}

fn holds_structure(view: AtomView<'_>) -> bool {
    let mut found = false;
    view.visitor(&mut |view| {
        if let AtomView::Fun(function) = view
            && structural(function.get_symbol())
        {
            found = true;
        }
        !found
    });
    found
}

/// No symbol, no function: a pure number.
fn numeric(view: AtomView<'_>) -> bool {
    let mut plain = true;
    view.visitor(&mut |view| {
        if matches!(view, AtomView::Var(_) | AtomView::Fun(_)) {
            plain = false;
        }
        plain
    });
    plain
}

/// A known constant function of pure numbers.
fn constant(symbol: Symbol) -> bool {
    [
        Symbol::LOG,
        Symbol::EXP,
        Symbol::SQRT,
        Symbol::SIN,
        Symbol::COS,
        Symbol::ABS,
        Symbol::CONJ,
    ]
    .contains(&symbol)
        || symbol.get_stripped_name() == "zeta"
}

fn short(view: AtomView<'_>) -> String {
    let text = view.to_owned().to_canonical_string();
    if text.chars().count() > 240 {
        format!("{}…", text.chars().take(240).collect::<String>())
    } else {
        text
    }
}

/// Lorentz slots in `d` dimensions (conventional dimensional
/// regularisation): `mink(n, i)` → `mink(d, i)`.
fn lorentz_in_d(expression: &Atom) -> Atom {
    let minkowski = Minkowski {}.to_symbolic([GS.dim]).get_symbol().unwrap();
    expression.replace_map(|view, _, out| {
        if let AtomView::Fun(function) = view
            && function.get_symbol() == minkowski
            && function.get_nargs() >= 1
        {
            let mut arguments = function.iter().map(|argument| argument.to_owned()).collect::<Vec<_>>();
            arguments[0] = Atom::var(GS.dim);
            **out = symbolica::atom::FunctionBuilder::new(minkowski)
                .add_args(&arguments)
                .finish();
        }
    })
}

/// Contract a closed tensor expression (every index summed) to a scalar in
/// `d = 4-2ε`; `None` when an index stays open, the algebra fails, or a
/// tensor structure remains.
pub fn contract_closed(expression: &Atom) -> Option<Atom> {
    let expression = lorentz_in_d(expression);
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let cooked = cooking.try_cook(expression.as_view()).ok()?;
    if !cooked.list_dangling::<AbstractIndex>().ok()?.is_empty() {
        return None;
    }
    let contracted = SymbolicTensor::infer_indexed(cooked)
        .ok()?
        .undo_trace()
        .ok()?
        .undo_chain()
        .ok()?
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::canonical()),
            color: Some(ColorSimplifySettings::default()),
            contract: AlgebraContraction::Dots,
            ..Default::default()
        })
        .ok()?;
    if !contracted
        .expression()
        .list_dangling::<AbstractIndex>()
        .ok()?
        .is_empty()
    {
        return None;
    }
    let scalar = cooking.uncook(contracted.expression().as_view());
    let scalar = Integrated::dimensionally_regularized(&scalar).expand();
    (!holds_structure(scalar.as_view())).then_some(scalar)
}

/// The normalised form of `view` (module doc); `denominator`: the view sits
/// in a denominator; `hard`: a factor of the hard part, integrated and
/// TRUNCATED before `d = 4-2ε` is substituted in it, so a `d` or a closed
/// trace in one of its denominators is a hidden pole of `h_j` that the
/// truncation cannot see (refused, not contracted).
fn walk(view: AtomView<'_>, denominator: bool, hard: bool) -> Normalised {
    let epsilon = GS.dim_epsilon;
    match view {
        AtomView::Num(_) => Ok(view.to_owned()),
        AtomView::Var(variable) => {
            let symbol = variable.get_symbol();
            if symbol == GS.dim {
                Ok(Atom::num(4) - Atom::num(2) * Atom::var(epsilon))
            } else if symbol == epsilon
                || !denominator
                || symbol == Symbol::PI
                || symbol == Symbol::E
            {
                Ok(view.to_owned())
            } else {
                Err((
                    Refusal::SymbolInDenominator,
                    symbol.get_name().to_owned(),
                ))
            }
        }
        AtomView::Fun(function) => {
            let symbol = function.get_symbol();
            if structural(symbol) {
                if denominator {
                    return Err((Refusal::UnresolvedTensor, short(view)));
                }
                // A tensor structure in a numerator: a polynomial factor (its
                // contractions give polynomials in d). Representation
                // dimensions are metadata.
                return Ok(view.to_owned());
            }
            if constant(symbol) && function.iter().all(numeric) {
                return Ok(view.to_owned());
            }
            if function.iter().any(holds_structure) {
                return Err((Refusal::UnresolvedTensor, short(view)));
            }
            Err((Refusal::OpaqueFunction, short(view)))
        }
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            let AtomView::Num(number) = exponent else {
                return Err((Refusal::OpaqueFunction, short(view)));
            };
            let number = number.get_coeff_view().to_owned();
            let integer = number.is_integer();
            if !integer && !numeric(base) {
                let dependent = base.contains_symbol(epsilon) || base.contains_symbol(GS.dim);
                if dependent || holds_structure(base) {
                    return Err((Refusal::OpaqueFunction, short(view)));
                }
            }
            let inverted = match &number {
                symbolica::coefficient::Coefficient::Complex(value) => {
                    value.is_real() && value.re.is_negative()
                }
                _ => false,
            };
            let inside = denominator != inverted;
            if inverted && hard && holds_structure(base) {
                return Err((Refusal::UnresolvedTensor, short(base)));
            }
            if inverted && hard && base.contains_symbol(GS.dim) {
                // A denominator in d of the hard part: admitted only when it
                // does not vanish at d = 4 (then the truncation of h_j in a
                // symbolic d loses nothing), e.g. 1/(d-1), 1/(d+2).
                let at_four = Integrated::dimensionally_regularized(&base.to_owned());
                if super::depth::excess_poles(&(Atom::one() / at_four)) != Some(0) {
                    return Err((
                        Refusal::SymbolInDenominator,
                        format!(
                            "d in a denominator of the hard part that vanishes at d = 4 (h_j is truncated before d = 4-2ε): {}",
                            short(view)
                        ),
                    ));
                }
            }
            let base = if inverted && holds_structure(base) {
                // A tensor structure in a denominator: only a closed trace,
                // contracted, is admitted.
                contract_closed(&base.to_owned())
                    .ok_or_else(|| (Refusal::UnresolvedTensor, short(base)))?
            } else {
                base.to_owned()
            };
            let base = walk(base.as_view(), inside, hard)?;
            Ok(base.pow(exponent.to_owned()))
        }
        AtomView::Mul(product) => {
            let mut result = Atom::one();
            for factor in product.iter() {
                result *= walk(factor, denominator, hard)?;
            }
            Ok(result)
        }
        AtomView::Add(sum) => {
            let terms = sum
                .iter()
                .map(|term| walk(term, denominator, hard))
                .collect::<std::result::Result<Vec<_>, _>>()?;
            Ok(Atom::add_many(terms.iter()))
        }
    }
}

/// The normalised form of one factor of the list (module doc): every closed
/// trace in a denominator contracted, `d = 4-2ε` in every scalar; refused
/// in the classes of [`Refusal`].
pub fn normalise(factor: Factor, expression: &Atom) -> std::result::Result<Atom, Refused> {
    let hard = factor == Factor::HardNumerator;
    walk(expression.as_view(), false, hard).map_err(|(class, detail)| Refused {
        class,
        factor: entry(factor).name,
        detail,
    })
}

/// `normalise` for a factor named by its log name; a name outside
/// [`FACTORS`] is refused (`outside_list`).
pub fn normalise_named(name: &str, expression: &Atom) -> std::result::Result<Atom, Refused> {
    match FACTORS.iter().find(|entry| entry.name == name) {
        Some(entry) => normalise(entry.factor, expression),
        None => Err(Refused {
            class: Refusal::OutsideList,
            factor: "unlisted",
            detail: name.to_owned(),
        }),
    }
}

/// The excess (`depth::excess_poles`) of a listed factor after
/// normalisation; `Ok(None)` when its order cannot be read (the depth
/// module's `unbounded`).
pub fn certified_excess(
    factor: Factor,
    expression: &Atom,
) -> std::result::Result<Option<usize>, Refused> {
    Ok(super::depth::excess_poles(&normalise(factor, expression)?))
}

/// A graph weight is regular at `d = 4` (beta's public boundary): its
/// normalised excess is zero. `Err` names the refusal or the pole order.
pub fn ensure_regular_weight(name: &str, expression: &Atom) -> std::result::Result<(), String> {
    match normalise_named(name, expression) {
        Err(refused) => Err(refused.to_string()),
        Ok(normalised) => match super::depth::excess_poles(&normalised) {
            Some(0) => Ok(()),
            Some(poles) => Err(format!(
                "IRH-DEPTH regularity: {name} has a pole of order {poles} at d = 4 ({}); a graph weight must be regular (the counterterm's finite orders are not kept)",
                short(expression.as_view())
            )),
            None => Err(format!(
                "IRH-DEPTH regularity: the ε order of {name} cannot be read ({})",
                short(expression.as_view())
            )),
        },
    }
}

/// The factors of one IR correction (the run-time entries of [`FACTORS`]).
pub(super) struct CorrectionFactors<'a> {
    /// `Graph::numerator(support, hard)`: the quotient's vertex and edge
    /// numerators.
    pub quotient: &'a Atom,
    /// `Graph::numerator(hard, ∅)`.
    pub hard: &'a Atom,
    /// `Graph::global_prefactor.num`.
    pub prefactor: &'a Atom,
    /// `Graph::global_prefactor.projector`.
    pub projector: &'a Atom,
}

/// Certify one correction (module doc): the normalised product of the
/// factors whose excess enters `B` (quotient numerators, prefactor,
/// projector; the root closure and the rows are read by the depth module),
/// or the first refusal.
pub(super) fn certify_correction(
    factors: &CorrectionFactors<'_>,
) -> std::result::Result<Atom, Refused> {
    let quotient = normalise(Factor::QuotientVertexNumerator, factors.quotient)?;
    // The hard numerator: its grammar only (its scalar poles are h_j's).
    let _ = normalise(Factor::HardNumerator, factors.hard)?;
    let prefactor = normalise(Factor::GraphNumeratorPrefactor, factors.prefactor)?;
    let projector = normalise(Factor::GraphProjector, factors.projector)?;
    Ok(quotient * prefactor * projector)
}

/// [`certify_correction`] for one correction of `graph`: `quotient` is
/// `Graph::numerator(support, hard)` (already built by the caller), the
/// hard numerator, prefactor and projector are read from `graph`.
pub(super) fn certify_graph_correction(
    graph: &crate::graph::Graph,
    hard: &linnet::half_edge::subgraph::SuBitGraph,
    quotient: &Atom,
) -> eyre::Result<std::result::Result<Atom, Refused>> {
    use crate::uv::UltravioletGraph;
    #[cfg(test)]
    if ROUND3.with(|slot| slot.get()) {
        return Ok(Ok(quotient.clone()));
    }
    let empty: linnet::half_edge::subgraph::SuBitGraph = graph.underlying.empty_subgraph();
    let hard_numerator = graph
        .numerator(hard, &empty)
        .get_single_atom()
        .map_err(|error| eyre::eyre!("the hard numerator: {error}"))?;
    Ok(certify_correction(&CorrectionFactors {
        quotient,
        hard: &hard_numerator,
        prefactor: &graph.global_prefactor.num,
        projector: &graph.global_prefactor.projector,
    }))
}

#[cfg(test)]
thread_local! {
    static ROUND3: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Tests: run `f` without this module's certificate on this thread (the
/// round-3 planner: `B` read on the factors as written, the graph
/// prefactor not read), to show each witness failing before round 5.
#[cfg(test)]
pub(super) fn as_round3<T>(f: impl FnOnce() -> T) -> T {
    struct Restore(bool);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0;
            ROUND3.with(|slot| slot.set(previous));
        }
    }
    let _restore = Restore(ROUND3.with(|slot| slot.replace(true)));
    f()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::initialisation::test_initialise;
    use spenso::{
        network::library::symbolic::ETS,
        structure::representation::RepName,
    };
    use symbolica::parse;

    fn eps() -> Atom {
        Atom::var(GS.dim_epsilon)
    }

    fn metric(a: usize, b: usize, dimension: impl Into<spenso::structure::dimension::Dimension>) -> Atom {
        let mink = Minkowski {}.new_rep(dimension);
        ETS.metric(
            mink.to_symbolic([crate::numerator::aind::Aind::Hedge(a as u16, 1).into()]),
            mink.to_symbolic([crate::numerator::aind::Aind::Hedge(b as u16, 1).into()]),
        )
    }

    fn excess(expression: &Atom) -> std::result::Result<Option<usize>, Refusal> {
        certified_excess(Factor::QuotientVertexNumerator, expression).map_err(|refused| refused.class)
    }

    #[test]
    fn the_list_is_closed() {
        let mut names = std::collections::BTreeSet::new();
        for entry in FACTORS {
            assert!(names.insert(entry.name), "{} twice", entry.name);
            assert!(!entry.multiplied.is_empty() && !entry.certified.is_empty());
            assert!(!entry.witness.is_empty());
            assert_eq!(super::entry(entry.factor).name, entry.name);
        }
        assert_eq!(
            normalise_named("vertex_weight", &Atom::one()).map_err(|r| r.class),
            Err(Refusal::OutsideList)
        );
    }

    /// The mechanical matrix on one factor: ε⁻¹ in three spellings (explicit,
    /// a `1/(d-4)` form, a closed metric trace with fresh dummy indices,
    /// in 4 or in d slots) has excess 1; their regular neighbours 0.
    #[test]
    fn every_spelling_of_a_pole_is_read() -> eyre::Result<()> {
        test_initialise()?;
        let d = Atom::var(GS.dim);
        let x = parse!("x");
        for four in [true, false] {
            let trace = if four {
                metric(81, 82, 4)
            } else {
                metric(81, 82, GS.dim)
            }
            .pow(2);
            let hidden = Atom::num(-2) / (&trace - Atom::num(4));
            assert_eq!(excess(&(&x * &hidden)), Ok(Some(1)), "{hidden}");
            let neighbour = Atom::num(-2) / (&trace - Atom::num(3));
            assert_eq!(excess(&(&x * &neighbour)), Ok(Some(0)), "{neighbour}");
            // A closed trace in a numerator is a polynomial in d.
            assert_eq!(excess(&(&x * &trace)), Ok(Some(0)));
        }
        assert_eq!(excess(&(&x / eps())), Ok(Some(1)));
        assert_eq!(excess(&(&x * Atom::num(-2) / (&d - Atom::num(4)))), Ok(Some(1)));
        // An open index in a denominator, a symbol, an opaque function.
        let open = metric(81, 82, GS.dim);
        assert_eq!(excess(&(&x / (&open - Atom::num(4)))), Err(Refusal::UnresolvedTensor));
        assert_eq!(excess(&(&x / (parse!("y") - Atom::num(4)))), Err(Refusal::SymbolInDenominator));
        let opaque = symbolica::function!(symbolica::symbol!("irh_regularity_f"), parse!("y"));
        assert_eq!(excess(&(&x * opaque)), Err(Refusal::OpaqueFunction));
        assert_eq!(excess(&(&x * eps().pow(parse!("y")))), Err(Refusal::OpaqueFunction));
        assert_eq!(
            excess(&(&x * eps().pow(Atom::num(1) / Atom::num(2)))),
            Err(Refusal::OpaqueFunction)
        );
        // Constants and numerator symbols are admitted.
        assert_eq!(excess(&(&x * parse!("zeta(3)") / (&d - Atom::num(1)))), Ok(Some(0)));
        // The hard part is truncated before d = 4-2ε: a d-denominator that
        // vanishes at d = 4, or a trace in a denominator, is refused there; a
        // regular one (1/(d-1)) and a d or a trace in a numerator are not.
        let hard = |expression: &Atom| {
            certified_excess(Factor::HardNumerator, expression).map_err(|refused| refused.class)
        };
        let trace = metric(81, 82, GS.dim).pow(2);
        assert_eq!(hard(&(&x / (&d - Atom::num(4)))), Err(Refusal::SymbolInDenominator));
        assert_eq!(hard(&(&x / (&trace - Atom::num(4)))), Err(Refusal::UnresolvedTensor));
        assert_eq!(hard(&(&x / (&d - Atom::num(1)))), Ok(Some(0)));
        assert_eq!(hard(&(&x * (&d - Atom::num(4)) * &trace)), Ok(Some(0)));
        assert_eq!(hard(&(&x / eps())), Ok(Some(1)));
        Ok(())
    }
}
