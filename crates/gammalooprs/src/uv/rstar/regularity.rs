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
    AdditionalNormalisation,
    VakintNormalisation,
    IrDiracClosure,
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
pub const FACTORS: [Entry; 17] = [
    Entry {
        factor: Factor::QuotientVertexNumerator,
        name: "quotient_vertex_numerator",
        multiplied: "Graph::numerator(support, hard) (numerator/uninit.rs fill_in_reduced: graph[node].get_num()), the integrand of every forest step of Z̃_F(G/γ; E_j) (Evaluation::ir_forest, quotient)",
        certified: "round 8b: its declared order k_v (normalise + excess of THIS vertex numerator alone, d = 4-2ε; 0 by type without ε) at run time, per correction, before H(γ) is integrated (regularity::certify_graph_correction); k_v > 0 declines the child route (declared_order); a scalar d is 4-2ε before every expansion of the forest steps (C1)",
        witness: "depth_tests::round8b::round8b_b_quotient_vertex_spellings_equal_the_literal",
    },
    Entry {
        factor: Factor::QuotientEdgeNumerator,
        name: "quotient_edge_numerator",
        multiplied: "Graph::numerator(support, hard) (fill_in_reduced: e.data.num), as the vertex numerators",
        certified: "as quotient_vertex_numerator, per edge (k_e of this edge numerator alone)",
        witness: "depth_tests::round8b::round8b_c_quotient_edge_spellings_equal_the_literal",
    },
    Entry {
        factor: Factor::HardNumerator,
        name: "hard_numerator",
        multiplied: "Graph::numerator(hard, ∅): the integrand of H(γ) (Evaluation::integrate_ir_hard_term); its tensor leaves are the insertions E_j",
        certified: "normalise at run time before H(γ) is integrated (grammar of every possible E_j); its explicit scalar poles stay in h_j, exact through ε^(B-1) because d = 4-2ε before the single-mass truncation (round 8b, C1: integrate.rs single_mass_truncation); a d-denominator that vanishes at d = 4 (symbol_in_denominator) or any tensor in a denominator (unresolved_tensor, not contracted) is still refused",
        witness: "depth_tests::round8b::round8b_a_hard_vertex_spellings_equal_the_literal",
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
        certified: "as quotient_vertex_numerator / hard_numerator (it IS a vertex numerator); regular by construction: an N-loop insertion over an N-loop scalar chain, truncated at ε^(N-1) after d = 4-2ε (round 8b, C1: insertion_coefficient_truncation), so k = 0; insertion group weights are rejected with a pole (PropagatorInsertionEvaluator::validate_weight, a weight with a tensor structure normalised here first: MATRIX-R45)",
        witness: "depth::tests::a_surrogate_coefficient_has_no_excess",
    },
    Entry {
        factor: Factor::GraphNumeratorPrefactor,
        name: "numerator_prefactor",
        multiplied: "not by R*: the caller's root projector (beta src/integrate.rs project_physical_root_mode: tensor · σ(global_prefactor.num), and closed_root_scalar_normalization), i.e. onto h_j z_j after the outer Taylor step",
        certified: "round 8b + 8a (integration): its declared order k_prefactor = NormalisedWeight::declared_k of regularity::normalise_weight (the one weight normaliser, the prefactor alone, never multiplied into another factor) when a root projector applies it: 0 for a regular weight; a pole, a non-finite, undecidable or zero value is refused (the correction, in its class); not read by the unweighted API, which does not apply it; beta normalises it at entry with the same normaliser and refuses a non-regular one there",
        witness: "depth_tests::round8b::round8b_d_a_regular_prefactor_is_no_strict_alarm",
    },
    Entry {
        factor: Factor::GraphProjector,
        name: "projector",
        multiplied: "not by R* nor by beta's R* path (beta only multiplies it in its direct integration, src/inputs.rs)",
        certified: "round 8b: NOT read (it is applied by no R* route: the review's finding 4, a zero projector erased every quotient excess); beta refuses a non-regular one at its public boundary",
        witness: "depth_tests::round8b::round8b_d_the_stored_projector_moves_nothing",
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
        certified: "round 8b: its declared order k_closure (normalise + excess at run time); k_closure > 0 declines",
        witness: "depth::tests::the_bound_adds_the_quotient_loops_and_every_excess",
    },
    Entry {
        factor: Factor::RootNormalisation,
        name: "root_normalisation",
        multiplied: "the caller's closed-root normalisation c = transverse·σ(prefactor)/(i q²)/(d−1) (beta closed_root_scalar_normalization) after Taylor",
        certified: "σ(prefactor) as numerator_prefactor; i q² is the external invariant (ε-free), d−1 is regular at d = 4; beta builds it from these alone. On the child route the closed-root factor is applied through the root projector callback, so it is measured as root_projector (MATRIX-R45 review)",
        witness: "depth_tests::round5::the_fixed_normalisations_are_regular",
    },
    Entry {
        factor: Factor::RootProjector,
        name: "root_projector",
        multiplied: "the caller's root projector callback on the complete outer polynomial (Evaluation::ir_outer_integrated)",
        certified: "round 8b (release ruling R1): only in its TYPED form X ↦ factor · σ(X) (ClosedRootNormalization given as data by the _scalar entry points: beta's closed metric root, c = transverse·σ(prefactor)/(i q²)/(d−1); gammaloop's own UV child closures 1/(d(d−1)(d+2)), 1/(Tr[1] d), 1/Tr[1], and the surrogate closed root 1): k_projector is the order of factor at d = 4-2ε, and the child route APPLIES this form (Evaluation::child_root_projection), so the certified scalar is the one multiplied; an opaque callback is never probed: the child route declines (opaque_projector)",
        witness: "depth_tests::round8b::round8b_e_an_opaque_callback_declines_without_a_probe",
    },
    Entry {
        factor: Factor::RowInsertion,
        name: "row_insertion",
        multiplied: "E_j of each row, the decorated quotient's insertion (Evaluation::ir_forest)",
        certified: "round 8b (C4): ε-free by construction (a leaf of the hard numerator's tensor structure, its scalars in h_j: moments::hard_moments), checked on every returned row (depth::eps_free); a row that is not declines (insertion_not_eps_free); a dropped row has a vanishing coefficient at depth B",
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
        certified: "the complete coefficient of a polynomial tensor integrand is regular at d = 4 (module doc); its premise is hard_numerator's normalisation. An intermediate coefficient is not (a rank-ten Weingarten coefficient has 1/(d-4)): an R* insertion coefficient's Laurent order is also read after d = 4-2ε before its closure integral (integrate.rs insertion_coefficient_valuation), and the outer coefficient of every row is measured after the tensor finish and the root projector (Evaluation::ir_outer_integrated)",
        witness: "depth_tests::round5::a_rank_ten_reduction_coefficient_pole_is_read",
    },
    Entry {
        factor: Factor::AdditionalNormalisation,
        name: "additional_normalisation",
        multiplied: "uv_settings.vakint.additional_normalization, parsed and multiplied as N^loops into every vacuum integral: the quotient K's and the hard part H(γ) (approx/integrated.rs to_vakint, rstar/expression.rs physical remainder)",
        certified: "normalise (hard semantics: it enters H(γ) before the truncation) at run time per correction; it must be regular (excess 0), else the correction is refused (pole_in_fixed_factor)",
        witness: "depth_tests::round5::a_poleful_additional_normalisation_is_refused",
    },
    Entry {
        factor: Factor::VakintNormalisation,
        name: "vakint_normalisation",
        multiplied: "uv_settings.vakint.normalization, the named loop normalisation vakint multiplies into every vacuum integral (settings.rs true_settings)",
        certified: "the named schemes (MSbar, FMFTandMATAD, GScheme, pySecDec) are fixed code, regular at d = 4 (π, exp and Γ factors of ε per loop); any other (custom) string is refused on the child route (outside_list)",
        witness: "depth_tests::round5::a_poleful_additional_normalisation_is_refused",
    },
    Entry {
        factor: Factor::IrDiracClosure,
        name: "ir_dirac_closure",
        multiplied: "GAMMALOOP_IR_DIRAC_CLOSURE (on by default): the closure γ(column,row,μ) multiplied into every remainder of H(γ) (times_momentum_free), and its normalisation 1/Tr[1] or 1/(Tr[1] d) multiplied into h_j after the truncation (ir_dirac::reassemble)",
        certified: "normalise (hard semantics) at run time per correction; the normalisation must be regular (excess 0), else the correction is refused (pole_in_fixed_factor); the representation dimensions behind Tr[1] are polynomial slots (normalise)",
        witness: "depth_tests::round5::a_poleful_dirac_closure_normalisation_is_refused",
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
    /// A fixed factor that must be regular (a normalisation multiplied into
    /// every vacuum integral or into `h_j` after the truncation) has a pole
    /// or an order that cannot be read.
    PoleInFixedFactor,
    /// A division by zero or a non-finite number (symbolica's infinity or an
    /// indeterminate): only the weight normaliser ([`normalise_weight`])
    /// raises it (review N1, round 8a).
    NonFinite,
}

impl Refusal {
    /// The `reason=` label (and the `IR_HARD_DEPTH_SUMMARY` field).
    pub fn label(self) -> &'static str {
        match self {
            Self::UnresolvedTensor => "unresolved_tensor",
            Self::OpaqueFunction => "opaque_function",
            Self::SymbolInDenominator => "symbol_in_denominator",
            Self::OutsideList => "outside_list",
            Self::PoleInFixedFactor => "pole_in_fixed_factor",
            Self::NonFinite => "non_finite",
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

/// An argument of a tensor structure that cannot hide a pole: a number, a
/// symbol, a polynomial in symbols with numeric coefficients (non-negative
/// integer powers, no function), or a tensor structure whose every argument
/// is a slot.
fn slot(view: AtomView<'_>) -> bool {
    match view {
        AtomView::Fun(function) => {
            structural(function.get_symbol()) && function.iter().all(slot)
        }
        view => polynomial(view),
    }
}

/// A polynomial in symbols with numeric coefficients.
fn polynomial(view: AtomView<'_>) -> bool {
    match view {
        AtomView::Num(_) | AtomView::Var(_) => true,
        AtomView::Fun(_) => false,
        AtomView::Add(sum) => sum.iter().all(polynomial),
        AtomView::Mul(product) => product.iter().all(polynomial),
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            let AtomView::Num(number) = exponent else {
                return false;
            };
            let number = number.get_coeff_view().to_owned();
            let negative = match &number {
                symbolica::coefficient::Coefficient::Complex(value) => {
                    !value.is_real() || value.re.is_negative()
                }
                _ => true,
            };
            number.is_integer() && !negative && polynomial(base)
        }
    }
}

/// Whether `expression` holds a tensor structure (a spenso tensor,
/// representation or index symbol) anywhere.
pub fn has_tensor_structure(expression: &Atom) -> bool {
    holds_structure(expression.as_view())
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
/// truncation cannot see (refused, not contracted). `finite` (the weight
/// normaliser only; the correction factors keep `false`): a non-finite
/// number, a denominator that normalises to zero, or a sub-expression whose
/// normal form is non-finite is refused AT that sub-expression
/// ([`Refusal::NonFinite`]), so a 1/0 cannot vanish under a further
/// inversion (`(1 + 1/0)^-1` is 0 in symbolica).
fn walk(view: AtomView<'_>, denominator: bool, hard: bool, finite: bool) -> Normalised {
    let result = walk_inner(view, denominator, hard, finite)?;
    if finite && !result.is_finite() {
        return Err((
            Refusal::NonFinite,
            format!("{} normalises to {}", short(view), short(result.as_view())),
        ));
    }
    Ok(result)
}

fn walk_inner(view: AtomView<'_>, denominator: bool, hard: bool, finite: bool) -> Normalised {
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
                // contractions give polynomials in d) when every argument is
                // a slot (`slot`): a number, a symbol, a polynomial (a
                // representation dimension d or 4-2ε) or a structure of
                // slots. Any other argument, a scalar argument of a tensor or
                // a dimension such as -2/(d-4) or 1/ε (a trace of the unit or
                // a closure normalisation 1/Tr[1] turns it into a pole), is
                // refused (regularity review of MATRIX-R45, findings 3, 4).
                if let Some(argument) = function.iter().find(|argument| !slot(*argument)) {
                    return Err((
                        Refusal::UnresolvedTensor,
                        format!(
                            "a tensor argument that is not a slot: {} in {}",
                            short(argument),
                            short(view)
                        ),
                    ));
                }
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
            let base = walk(base.as_view(), inside, hard, finite)?;
            if finite && inverted && (base.is_zero() || base.together().is_zero()) {
                return Err((
                    Refusal::NonFinite,
                    format!(
                        "a division by zero: the denominator of {} is 0 at d = 4-2ε",
                        short(view)
                    ),
                ));
            }
            Ok(base.pow(exponent.to_owned()))
        }
        AtomView::Mul(product) => {
            let mut result = Atom::one();
            for factor in product.iter() {
                result *= walk(factor, denominator, hard, finite)?;
            }
            Ok(result)
        }
        AtomView::Add(sum) => {
            let terms = sum
                .iter()
                .map(|term| walk(term, denominator, hard, finite))
                .collect::<std::result::Result<Vec<_>, _>>()?;
            Ok(Atom::add_many(terms.iter()))
        }
    }
}

/// The normalised form of one factor of the list (module doc): every closed
/// trace in a denominator contracted, `d = 4-2ε` in every scalar; refused
/// in the classes of [`Refusal`].
pub fn normalise(factor: Factor, expression: &Atom) -> std::result::Result<Atom, Refused> {
    let hard = matches!(
        factor,
        Factor::HardNumerator | Factor::AdditionalNormalisation | Factor::IrDiracClosure
    );
    walk(expression.as_view(), false, hard, false).map_err(|(class, detail)| Refused {
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

/// A graph weight is regular at `d = 4` (beta's public boundary): it
/// normalises ([`normalise_weight`]). `Err` names the refusal or the pole
/// order. A thin wrapper: the check and the value are [`normalise_weight`]'s.
pub fn ensure_regular_weight(name: &str, expression: &Atom) -> std::result::Result<(), String> {
    normalise_weight(name, expression)
        .map(|_| ())
        .map_err(|refused| refused.to_string())
}

/// IRH-DEPTH round 8 (Astra r47, finding 2): a graph weight as ONE
/// normalised scalar, the object that is both checked and multiplied.
///
/// Only [`normalise_weight`] constructs one (and [`NormalisedWeight::one`],
/// the literal unit), so a value of this type has passed the contract:
/// * every closed metric or Dirac trace is contracted (fresh dummy indices
///   included, Lorentz slots in `d`), wherever it sits (numerator,
///   denominator, nested);
/// * every spelling of the dimension in use is `4 − 2ε`: gammaloop's (and
///   feynkit's) `dim`, a Lorentz slot `mink(n, ·)` of a contracted trace,
///   and vakint's FORM symbols `vakint::d` (= `4 − 2 vakint::ep`) and
///   `vakint::ep` (= `ε`);
/// * the result is simplified (`together`) to one rational function of `ε`
///   with constant coefficients (numbers, `π`, `e`, and the known constant
///   functions of numbers, e.g. `zeta(3)`), whose order at `ε = 0` is then
///   decided exactly on the whole scalar (not a sum of factor bounds:
///   `(d² − 16)/(d − 4) = d + 4` is regular);
/// * that order is `≥ 0` (regular).
///
/// REFUSED ([`WeightRefused`], named class, before any projection): a form
/// the normaliser cannot decide (an unknown function head, a symbolic or
/// non-integer exponent of an `ε`-dependent base, an open index or a
/// tensor that does not contract to a scalar, any symbol other than `ε`,
/// `π`, `e` and the dimension spellings), a pole (order `< 0`), an order
/// that cannot be read, a NON-FINITE value or a division by zero at any
/// stage (round 8a, review N1: `(d − 4 + 2ε)⁻¹` normalised to `⧞`, which
/// has no `ε` and read as order 0), a value whose zero cannot be decided
/// (a function, a non-integer power or an inexact number left after
/// simplification: the value must be a rational function of `ε` with exact
/// coefficients over `Q(i)`, `π` and `e`), and a weight that normalises to
/// exactly ZERO. A
/// zero is refused, not returned: a weight of 0 makes the whole diagram 0,
/// which is exactly how a lost value looks (the r47 finding-2 defect turned
/// the unit weight into a zero result); a graph whose weights cancel is
/// removed by the generator with a recorded reason, never by a zero weight
/// reaching the integration.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NormalisedWeight {
    value: Atom,
    order: i64,
}

impl NormalisedWeight {
    /// The literal unit (normalised by construction).
    pub fn one() -> Self {
        Self {
            value: Atom::one(),
            order: 0,
        }
    }

    /// The normalised scalar: the ONLY value a caller multiplies.
    pub fn value(&self) -> &Atom {
        &self.value
    }

    /// Its order at `ε = 0` (`≥ 0`).
    pub fn order(&self) -> i64 {
        self.order
    }

    /// Whether the value is free of `ε` (a constant).
    pub fn is_epsilon_free(&self) -> bool {
        !self.value.contains_symbol(GS.dim_epsilon)
    }

    /// Whether the value is exactly 1.
    pub fn is_one(&self) -> bool {
        self.value.is_one()
    }

    /// IRH-DEPTH 8a + 8b integration: the declared order `k` of this
    /// normalised external factor in the topological certificate
    /// ([`Certificate`]): its pole order `max(0, −order)`. A normalised
    /// weight is regular by contract (a pole is refused by
    /// [`normalise_weight`]), so this is 0; it is computed from the order,
    /// not assumed, so the certificate and the normaliser read ONE number.
    pub fn declared_k(&self) -> usize {
        usize::try_from((-self.order).max(0)).unwrap_or(usize::MAX)
    }
}

/// Why [`normalise_weight`] refused a weight.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum WeightRefusal {
    /// The normaliser does not understand the form (the [`Refusal`] class:
    /// unresolved tensor, opaque function, symbol in a denominator, outside
    /// the list).
    Unsupported(Refusal),
    /// A tensor structure that is not a closed trace (an open index in what
    /// must be a scalar) or that the algebra does not contract.
    OpenIndex,
    /// A symbol whose value is unknown (it could be `ε`).
    UnknownSymbol,
    /// A pole at `d = 4` (order `< 0`).
    Pole,
    /// The order at `ε = 0` cannot be read.
    Unreadable,
    /// The weight normalises to exactly zero.
    Zero,
    /// A division by zero, or a value that is not a finite number
    /// (symbolica's infinity `⧞` or an indeterminate), at any stage.
    NonFinite,
    /// A value whose zero (or pole) the normaliser cannot decide: a function
    /// (e.g. `zeta(3)`, `log(2)`), a non-integer power (`2^(1/2)`) or an
    /// inexact (floating-point) number left after simplification.
    Undecidable,
}

impl WeightRefusal {
    pub fn label(self) -> &'static str {
        match self {
            Self::Unsupported(class) => class.label(),
            Self::OpenIndex => "open_index",
            Self::UnknownSymbol => "unknown_symbol",
            Self::Pole => "pole",
            Self::Unreadable => "unreadable_order",
            Self::Zero => "zero",
            Self::NonFinite => "non_finite",
            Self::Undecidable => "undecidable",
        }
    }
}

/// A refused graph weight: its name, class and the detail.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WeightRefused {
    pub name: String,
    pub class: WeightRefusal,
    pub detail: String,
}

impl std::fmt::Display for WeightRefused {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.class {
            WeightRefusal::Pole => write!(
                f,
                "IRH-DEPTH regularity: {} has a pole of order {} at d = 4; a graph weight must be regular (the counterterm's finite orders are not kept)",
                self.name, self.detail
            ),
            class => write!(
                f,
                "IRH-DEPTH regularity: {} refused ({}): {}",
                self.name,
                class.label(),
                self.detail
            ),
        }
    }
}

impl std::error::Error for WeightRefused {}

/// The constant symbols a normalised weight may hold besides `ε`.
fn weight_symbol(symbol: Symbol) -> bool {
    symbol == GS.dim_epsilon || symbol == Symbol::PI || symbol == Symbol::E
}

/// Whether `value` is a rational function of `ε` (and `π`, `e`) whose
/// coefficients are exact, finite complex rationals: no function, no
/// non-integer power, no float, finite-field or nested rational-polynomial
/// coefficient, no infinity or indeterminate.
fn decidable(value: &Atom) -> std::result::Result<(), (WeightRefusal, String)> {
    use symbolica::coefficient::CoefficientView;
    let mut refusal = None;
    value.visitor(&mut |view| {
        if refusal.is_some() {
            return false;
        }
        match view {
            AtomView::Num(number) => match number.get_coeff_view() {
                CoefficientView::Natural(..) | CoefficientView::Large(..) => {}
                CoefficientView::Infinity(_) | CoefficientView::Indeterminate => {
                    refusal = Some((WeightRefusal::NonFinite, format!("the number {}", short(view))));
                }
                _ => {
                    refusal = Some((
                        WeightRefusal::Undecidable,
                        format!("the inexact number {}", short(view)),
                    ));
                }
            },
            AtomView::Fun(_) => {
                refusal = Some((WeightRefusal::Undecidable, format!("the function {}", short(view))));
            }
            AtomView::Pow(power) => {
                let (_, exponent) = power.get_base_exp();
                let integer = matches!(exponent, AtomView::Num(number)
                    if number.get_coeff_view().to_owned().is_integer());
                if !integer {
                    refusal = Some((
                        WeightRefusal::Undecidable,
                        format!("the non-integer power {}", short(view)),
                    ));
                }
            }
            AtomView::Var(variable) if !weight_symbol(variable.get_symbol()) => {
                refusal = Some((
                    WeightRefusal::UnknownSymbol,
                    format!("the symbol {}", variable.get_symbol().get_name()),
                ));
            }
            _ => {}
        }
        refusal.is_none()
    });
    match refusal {
        Some(refusal) => Err(refusal),
        None => Ok(()),
    }
}

/// THE weight normaliser ([`NormalisedWeight`]): `name` must be a graph
/// weight of the closed list ([`FACTORS`]: `overall_factor`,
/// `numerator_prefactor`, `projector`), else `outside_list`.
pub fn normalise_weight(
    name: &str,
    expression: &Atom,
) -> std::result::Result<NormalisedWeight, WeightRefused> {
    let refuse = |class: WeightRefusal, detail: String| WeightRefused {
        name: name.to_owned(),
        class,
        detail,
    };
    if !matches!(name, "overall_factor" | "numerator_prefactor" | "projector") {
        return Err(refuse(
            WeightRefusal::Unsupported(Refusal::OutsideList),
            format!("{name} is not a graph weight"),
        ));
    }
    // vakint's FORM spellings of the dimension and the regulator.
    let vakint = &vakint::symbols::S;
    let expression = expression
        .replace(vakint.form_dimension)
        .with(Atom::var(GS.dim))
        .replace(vakint.form_epsilon)
        .with(Atom::var(GS.dim_epsilon));
    let unsupported = |(class, detail): (Refusal, String)| match class {
        Refusal::NonFinite => refuse(WeightRefusal::NonFinite, detail),
        class => refuse(WeightRefusal::Unsupported(class), detail),
    };
    if !expression.is_finite() {
        return Err(refuse(
            WeightRefusal::NonFinite,
            format!("{} is not finite", short(expression.as_view())),
        ));
    }
    // Closed traces in denominators (and nested ones) contracted, d = 4-2ε.
    let mut normalised = walk(expression.as_view(), false, false, true).map_err(unsupported)?;
    if holds_structure(normalised.as_view()) {
        // A closed trace in a numerator: the whole scalar is contracted.
        let contracted = contract_closed(&normalised).ok_or_else(|| {
            refuse(
                WeightRefusal::OpenIndex,
                format!(
                    "a tensor structure that is not a closed trace contracting to a scalar: {}",
                    short(normalised.as_view())
                ),
            )
        })?;
        normalised = walk(contracted.as_view(), false, false, true).map_err(unsupported)?;
    }
    let mut unknown = None;
    normalised.visitor(&mut |view| {
        if let AtomView::Var(variable) = view
            && !weight_symbol(variable.get_symbol())
        {
            unknown = Some(variable.get_symbol());
        }
        unknown.is_none()
    });
    if let Some(symbol) = unknown {
        return Err(refuse(
            WeightRefusal::UnknownSymbol,
            format!(
                "{} in {} (a graph weight holds numbers, ε, π, e and constant functions of numbers only)",
                symbol.get_name(),
                short(expression.as_view())
            ),
        ));
    }
    let value = normalised.together();
    // The value must be a rational function of ε with exact, finite
    // coefficients over Q(i) and π, e: then its zero and its order are
    // decided exactly (review N1).
    decidable(&value).map_err(|(class, detail)| {
        refuse(
            class,
            format!("{detail} in {} (normalised {})", short(expression.as_view()), short(value.as_view())),
        )
    })?;
    if value.is_zero() {
        return Err(refuse(
            WeightRefusal::Zero,
            format!(
                "{} normalises to exactly 0 (a zero weight would hide a lost value; the generator removes cancelled graphs)",
                short(expression.as_view())
            ),
        ));
    }
    let order = if value.contains_symbol(GS.dim_epsilon) {
        use symbolica::poly::series::SeriesDepth;
        let unreadable = || {
            refuse(
                WeightRefusal::Unreadable,
                format!("the ε order of {} cannot be read", short(value.as_view())),
            )
        };
        let series = value
            .series(GS.dim_epsilon, 0, SeriesDepth::relative(1))
            .map_err(|_| unreadable())?;
        if series.is_zero() {
            return Err(unreadable());
        }
        let exponent = series.get_trailing_exponent();
        if !exponent.is_integer() {
            return Err(unreadable());
        }
        exponent.numerator().to_i64().ok_or_else(unreadable)?
    } else {
        0
    };
    if order < 0 {
        return Err(refuse(
            WeightRefusal::Pole,
            format!("{} ({})", -order, short(expression.as_view())),
        ));
    }
    Ok(NormalisedWeight { value, order })
}

/// IRH-DEPTH round 8b: one declared order of the topological certificate
/// (TOPOLOGICAL-BOUND.md section 1): a non-negative integer `k` per factor,
/// kept separately, never summed with signs (a positive order of one factor
/// never offsets another factor's pole: the review's finding 4).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct Declared {
    /// The [`FACTORS`] name of the factor kind.
    pub factor: &'static str,
    /// Which factor (`v3`, `e7`, `root`, …).
    pub label: String,
    pub k: usize,
}

/// How the root projection of an evaluation is given (IRH-DEPTH round 8b,
/// release ruling R1).
#[derive(Clone, Copy, Debug)]
pub(super) enum RootProjection<'a> {
    /// No root projector: the value is the unprojected counterterm.
    Absent,
    /// The typed regular form `X ↦ factor · σ(X)` (a fixed routing `σ` of the
    /// external momenta and a declared scalar normalisation); its `k` is the
    /// pole order of `factor` at `d = 4 − 2ε`.
    Typed(&'a Atom),
    /// An opaque callback: nothing certifies it; the child route declines.
    Opaque,
}

/// The topological certificate of one IR correction, computed from the graph
/// and the declared orders BEFORE any expression of the correction is built:
/// `B = L_q + Σ k`, `L_q = loops(support) − loops(hard)`.
#[derive(Clone, Debug)]
pub(super) struct Certificate {
    pub quotient_loops: usize,
    pub declared: Vec<Declared>,
}

impl Certificate {
    /// `B`.
    pub fn bound(&self) -> usize {
        self.quotient_loops + self.declared.iter().map(|declared| declared.k).sum::<usize>()
    }

    /// The first factor with a declared order above 0 (the child route then
    /// declines: TOPOLOGICAL-BOUND.md C2).
    pub fn positive(&self) -> Option<&Declared> {
        self.declared.iter().find(|declared| declared.k > 0)
    }

    /// The nonzero orders as a log field (`-` when every `k` is 0).
    pub fn label(&self) -> String {
        let positive = self
            .declared
            .iter()
            .filter(|declared| declared.k > 0)
            .map(|declared| format!("{}:{}={}", declared.factor, declared.label, declared.k))
            .collect::<Vec<_>>();
        if positive.is_empty() {
            "-".to_owned()
        } else {
            positive.join(",")
        }
    }
}

/// The declared order `k` of one listed factor (round 8b): normalised as in
/// round 5 (a closed trace in a denominator contracted, `d = 4 − 2ε`; an
/// unresolved tensor, an opaque function or a symbol in a denominator is
/// refused), then `max(0, −ord_ε)` of each of its multiplicative factors
/// read on its own ([`super::depth::excess_poles`]). A factor without `ε`
/// (every model Feynman rule: no ε, no d in 80,681 L5 DOTs) has `k = 0` by
/// type. `Err` with `PoleInFixedFactor` when the order cannot be read.
fn declared_order(factor: Factor, label: String, expression: &Atom) -> std::result::Result<Declared, Refused> {
    let k = certified_excess(factor, expression)?.ok_or_else(|| Refused {
        class: Refusal::PoleInFixedFactor,
        factor: entry(factor).name,
        detail: format!("{label}: the ε order cannot be read ({})", short(expression.as_view())),
    })?;
    Ok(Declared {
        factor: entry(factor).name,
        label,
        k,
    })
}

/// IRH-DEPTH 8a + 8b integration: the declared order of a graph WEIGHT (an
/// external factor of the closed list of [`normalise_weight`]: today the
/// numerator prefactor where a root projector applies it) is read through
/// the weight normaliser and nothing else: `k` is
/// [`NormalisedWeight::declared_k`] of the one [`NormalisedWeight`] (0 for a
/// regular weight), and every refusal of the normaliser (a pole, a
/// non-finite or undecidable value, an open index, an unknown symbol, a
/// zero) refuses the correction in its class. One path: the order the
/// certificate adds is the order of the object the caller multiplies (beta
/// stores the normalised value in `global_prefactor.num` at entry), never a
/// second reading of the spelling by [`super::depth::excess_poles`].
fn declared_weight(
    factor: Factor,
    weight: &'static str,
    label: String,
    expression: &Atom,
) -> std::result::Result<Declared, Refused> {
    let normalised = normalise_weight(weight, expression).map_err(|refused| Refused {
        class: match refused.class {
            WeightRefusal::Unsupported(class) => class,
            WeightRefusal::NonFinite => Refusal::NonFinite,
            WeightRefusal::OpenIndex => Refusal::UnresolvedTensor,
            WeightRefusal::Undecidable => Refusal::OpaqueFunction,
            WeightRefusal::UnknownSymbol => Refusal::OutsideList,
            WeightRefusal::Pole | WeightRefusal::Unreadable | WeightRefusal::Zero => {
                Refusal::PoleInFixedFactor
            }
        },
        factor: entry(factor).name,
        detail: format!("{label}: {refused}"),
    })?;
    Ok(Declared {
        factor: entry(factor).name,
        label,
        k: normalised.declared_k(),
    })
}

/// The vertex and edge numerators of the quotient `G/γ` one by one, as
/// `Graph::numerator(support, hard)` multiplies them
/// (`numerator/uninit.rs fill_in_reduced`): every edge of `support − hard`
/// with its end vertices not in `hard`, then the vertices whose crown lies in
/// `support`. Each is read before any product, so no factor's order can
/// offset another's.
fn quotient_factors(
    graph: &crate::graph::Graph,
    support: &linnet::half_edge::subgraph::SuBitGraph,
    hard: &linnet::half_edge::subgraph::SuBitGraph,
) -> Vec<(Factor, String, Atom)> {
    use linnet::half_edge::{
        NodeIndex,
        involution::HedgePair,
        subgraph::{Inclusion, ModifySubSet, SubSetLike, SubSetOps, subset::SubSet},
    };
    let mut factors = Vec::new();
    let mut seen: SubSet<NodeIndex> = SubSet::empty(graph.n_nodes());
    for (node, _, _) in graph.underlying.iter_nodes_of(hard) {
        seen.add(node);
    }
    let not_ignored = support.subtract(hard);
    for (pair, edge, data) in graph.underlying.iter_edges_of(&not_ignored) {
        if let HedgePair::Paired { source, sink } = pair {
            for end in [graph.node_id(source), graph.node_id(sink)] {
                if !seen[end] {
                    seen.add(end);
                    factors.push((
                        Factor::QuotientVertexNumerator,
                        format!("v{}", usize::from(end)),
                        Atom::one() * graph[end].get_num(),
                    ));
                }
            }
            factors.push((
                Factor::QuotientEdgeNumerator,
                format!("e{}", usize::from(edge)),
                data.data.num.value.clone(),
            ));
        }
    }
    let unseen = !seen;
    for node in unseen.included_iter() {
        if graph
            .iter_crown(node)
            .all(|hedge| support.includes(&hedge) || graph[graph[&hedge]].is_dummy)
        {
            factors.push((
                Factor::QuotientVertexNumerator,
                format!("v{}", usize::from(node)),
                Atom::one() * graph[node].get_num(),
            ));
        }
    }
    factors
}

/// The certificate of one IR correction (round 8b): `L_q` and the declared
/// orders of every factor that multiplies the quotient side of the product
/// `h_j · z_j · closure · compensation · root projection`, each on its own:
/// * every quotient vertex and edge numerator (`k_v`, the decorated
///   quotient's integrand);
/// * the root closure (`k_closure`);
/// * the root projection (`k_projector`): 0 without one; the declared scalar
///   normalisation's order for the typed form; an opaque callback is not
///   certified (the child route declines before this is used);
/// * the graph's numerator prefactor (`k_prefactor`), only where it is
///   applied: by the caller's root projector (beta's multiplies
///   `σ(global_prefactor.num)`); the unweighted API does not apply it.
///
/// The graph's stored projector is NOT read: no R* route and no beta R* path
/// applies it (the review's finding 4: a zero projector erased every
/// quotient excess). The hard numerator is checked for its grammar only (its
/// scalar poles are `h_j`'s, exact through `ε^(B−1)` once d = 4 − 2ε before
/// the truncation, C1). `ω`, the hard loops, the rows `E_j` (ε-free by
/// construction, asserted per row), the compensation and the tensor finish
/// add nothing (TOPOLOGICAL-BOUND.md section 2).
pub(super) fn certify_graph_correction(
    graph: &crate::graph::Graph,
    support: &linnet::half_edge::subgraph::SuBitGraph,
    hard: &linnet::half_edge::subgraph::SuBitGraph,
    root_closure: Option<&Atom>,
    projection: RootProjection<'_>,
) -> eyre::Result<std::result::Result<Certificate, Refused>> {
    use crate::uv::UltravioletGraph;
    let hard_loops = graph.n_loops(hard);
    let quotient_loops = graph.n_loops(support).saturating_sub(hard_loops);
    #[cfg(test)]
    if ROUND3.with(|slot| slot.get()) {
        // The round-3 planner: one order, read on the quotient numerator
        // as written (their product), the prefactor not read.
        let quotient = graph
            .numerator(support, hard)
            .get_single_atom()
            .map_err(|error| eyre::eyre!("the quotient numerator: {error}"))?;
        return Ok(Ok(Certificate {
            quotient_loops,
            declared: vec![Declared {
                factor: "round3",
                label: "quotient".to_owned(),
                k: super::depth::excess_poles(&quotient).unwrap_or(0),
            }],
        }));
    }
    let empty: linnet::half_edge::subgraph::SuBitGraph = graph.underlying.empty_subgraph();
    let hard_numerator = graph
        .numerator(hard, &empty)
        .get_single_atom()
        .map_err(|error| eyre::eyre!("the hard numerator: {error}"))?;
    Ok((|| {
        // The hard numerator: its grammar only.
        let _ = normalise(Factor::HardNumerator, &hard_numerator)?;
        let mut declared = Vec::new();
        for (factor, label, expression) in quotient_factors(graph, support, hard) {
            declared.push(declared_order(factor, label, &expression)?);
        }
        if let Some(closure) = root_closure {
            declared.push(declared_order(Factor::RootClosure, "root".to_owned(), closure)?);
        }
        match projection {
            RootProjection::Absent => {}
            RootProjection::Typed(factor) => {
                // The declared normalisation (beta: transverse · σ(prefactor)
                // / (i q²) / (d − 1); a UV child closure's 1/(d (d−1) (d+2)),
                // 1/(Tr[1] d), 1/Tr[1]): its scalar order at d = 4 − 2ε; an
                // external invariant in a denominator is ε-free.
                let k = super::depth::excess_poles(factor).ok_or_else(|| Refused {
                    class: Refusal::PoleInFixedFactor,
                    factor: entry(Factor::RootProjector).name,
                    detail: format!(
                        "the ε order of the typed root normalisation cannot be read ({})",
                        short(factor.as_view())
                    ),
                })?;
                declared.push(Declared {
                    factor: entry(Factor::RootProjector).name,
                    label: "typed".to_owned(),
                    k,
                });
                declared.push(declared_weight(
                    Factor::GraphNumeratorPrefactor,
                    "numerator_prefactor",
                    "prefactor".to_owned(),
                    &graph.global_prefactor.num,
                )?);
            }
            RootProjection::Opaque => {
                declared.push(declared_weight(
                    Factor::GraphNumeratorPrefactor,
                    "numerator_prefactor",
                    "prefactor".to_owned(),
                    &graph.global_prefactor.num,
                )?);
            }
        }
        Ok(Certificate {
            quotient_loops,
            declared,
        })
    })())
}

/// A fixed factor that must be regular ([`Refusal::PoleInFixedFactor`]):
/// normalised (refused in its class), then its excess must be `Some(0)`.
pub fn certify_regular(factor: Factor, expression: &Atom) -> std::result::Result<(), Refused> {
    let normalised = normalise(factor, expression)?;
    match super::depth::excess_poles(&normalised) {
        Some(0) => Ok(()),
        excess => Err(Refused {
            class: Refusal::PoleInFixedFactor,
            factor: entry(factor).name,
            detail: format!(
                "excess {} ({})",
                excess.map_or("unreadable".to_owned(), |excess| excess.to_string()),
                short(expression.as_view())
            ),
        }),
    }
}

/// The named loop normalisations of vakint (entry `vakint_normalisation`).
pub const NAMED_NORMALISATIONS: [&str; 4] = ["MSbar", "FMFTandMATAD", "GScheme", "pySecDec"];

/// The settings factors of one correction (entries `additional_normalisation`
/// and `vakint_normalisation`): the free-form `N` multiplied as `N^loops`
/// into every vacuum integral must parse and be regular; the named loop
/// normalisation must be one of the fixed schemes.
pub fn certify_vakint_normalisations(
    additional: &str,
    named: &str,
) -> std::result::Result<(), Refused> {
    if !NAMED_NORMALISATIONS.contains(&named) {
        return Err(Refused {
            class: Refusal::OutsideList,
            factor: entry(Factor::VakintNormalisation).name,
            detail: format!("a custom loop normalisation: {named}"),
        });
    }
    let parsed = symbolica::try_parse!(additional).map_err(|error| Refused {
        class: Refusal::OpaqueFunction,
        factor: entry(Factor::AdditionalNormalisation).name,
        detail: format!("{additional}: {error}"),
    })?;
    certify_regular(Factor::AdditionalNormalisation, &parsed)
}

/// The IR Dirac closure of one correction (entry `ir_dirac_closure`): the
/// closure multiplied into `H(γ)` is normalised with the hard part's
/// semantics; its normalisation, multiplied into `h_j` after the truncation,
/// must be regular.
pub fn certify_dirac_closure(
    closure: &Atom,
    normalization: &Atom,
) -> std::result::Result<(), Refused> {
    let _ = normalise(Factor::IrDiracClosure, closure)?;
    certify_regular(Factor::IrDiracClosure, normalization)
}

#[cfg(test)]
thread_local! {
    static ROUND3: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Tests: whether [`as_round3`] is active on this thread.
#[cfg(test)]
pub(super) fn is_round3() -> bool {
    ROUND3.with(|slot| slot.get())
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
