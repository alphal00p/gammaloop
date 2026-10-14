use std::{
    cell::RefCell,
    collections::{BTreeMap, HashMap},
    sync::LazyLock,
};

use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, Orientation};

use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    spenso_print_scripted_indexed,
    structure::{
        abstract_index::AIND_SYMBOLS,
        concrete_index::ExpandedIndex,
        representation::{Minkowski, RepName},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    domains::rational::Rational,
    function, get_symbol,
    id::Replacement,
    printer::{PrintState, PrintUserData},
    symbol,
};
use symbolica_utils::{PrintSettingsExt, TypstMode};

use crate::{cff::orientations::GraphOrientation, graph::LoopMomentumBasis};

/// Persistent ownership carried by an EMR momentum through a local 4D UV
/// Taylor expansion. Wire value two is deliberately reserved for the
/// transient denominator marker used while the Taylor operator is acting.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(i64)]
pub(crate) enum UvMomentumProvenanceRole {
    /// A pre-existing numerator factor fixed to its original owner during raw
    /// Taylor recursion. Completed UV projection may certify a denominator
    /// class and transfer this factor to that class's occurrence pool.
    TaylorFixed = 0,
    /// A hard numerator factor created by differentiating a denominator. Raw
    /// provenance permits copies of the same owner; completed UV projection
    /// may certify equivalent signed channels before class-wide dispatch.
    DenominatorDerived = 1,
    /// A physical source momentum reconstructed from its stored hard lift,
    /// independently of occurrence-local contact samples.
    PhysicalSourceFixed = 3,
    /// A soft momentum created by a denominator Taylor derivative. Its exact
    /// off-shell routing may be chosen when the outer numerator is complete;
    /// it is not an occurrence of the differentiated hard denominator.
    DenominatorDerivedSoft = 4,
}

/// An algebraic denominator channel in one completed UV projection view.
/// This namespace is distinct from physical source edges and CFF occurrences.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct UvDenominatorClassId(pub usize);

impl From<bool> for UvMomentumProvenanceRole {
    fn from(denominator_derived: bool) -> Self {
        if denominator_derived {
            Self::DenominatorDerived
        } else {
            Self::TaylorFixed
        }
    }
}

pub struct WildCards {
    pub edgeid_: Symbol,
    pub mom_: Symbol,
    pub mass_: Symbol,
    pub prop_: Symbol,
    pub a_: Symbol,
    pub b_: Symbol,
    pub c_: Symbol,
    pub d_: Symbol,
    pub e_: Symbol,
    pub f_: Symbol,
    pub g_: Symbol,
    pub h_: Symbol,
    pub i_: Symbol,
    pub j_: Symbol,
    pub k_: Symbol,
    pub l_: Symbol,
    pub m_: Symbol,
    pub n_: Symbol,
    pub o_: Symbol,
    pub p_: Symbol,
    pub q_: Symbol,
    pub r_: Symbol,
    pub s_: Symbol,
    pub t_: Symbol,
    pub u_: Symbol,
    pub v_: Symbol,
    pub w_: Symbol,
    pub x_: Symbol,
    pub y_: Symbol,
    pub z_: Symbol,
    pub a__: Symbol,
    pub b__: Symbol,
    pub c__: Symbol,
    pub d__: Symbol,
    pub e__: Symbol,
    pub f__: Symbol,
    pub g__: Symbol,
    pub h__: Symbol,
    pub i__: Symbol,
    pub j__: Symbol,
    pub k__: Symbol,
    pub l__: Symbol,
    pub m__: Symbol,
    pub n__: Symbol,
    pub o__: Symbol,
    pub p__: Symbol,
    pub q__: Symbol,
    pub r__: Symbol,
    pub s__: Symbol,
    pub t__: Symbol,
    pub u__: Symbol,
    pub v__: Symbol,
    pub w__: Symbol,
    pub x__: Symbol,
    pub y__: Symbol,
    pub z__: Symbol,
    pub a___: Symbol,
    pub b___: Symbol,
    pub c___: Symbol,
    pub d___: Symbol,
    pub e___: Symbol,
    pub f___: Symbol,
    pub g___: Symbol,
    pub h___: Symbol,
    pub i___: Symbol,
    pub j___: Symbol,
    pub k___: Symbol,
    pub l___: Symbol,
    pub m___: Symbol,
    pub n___: Symbol,
    pub o___: Symbol,
    pub p___: Symbol,
    pub q___: Symbol,
    pub r___: Symbol,
    pub s___: Symbol,
    pub t___: Symbol,
    pub u___: Symbol,
    pub v___: Symbol,
    pub w___: Symbol,
    pub x___: Symbol,
    pub y___: Symbol,
    pub z___: Symbol,
}

pub struct GammaloopSymbols {
    pub integrand: Symbol,
    /// Exact generalized-3D-representation residue-map identifier selected by
    /// the runtime evaluator. This is independent of physical edge signs.
    pub residue_map_id: Symbol,
    /// wrapper function for the 4d bridge denominators
    pub tree_denom_wrapper: Symbol,

    pub killing_func: Symbol,
    pub is_function: Symbol,
    pub is_symbol: Symbol,
    pub ufozero: Symbol,
    pub _linear: Symbol,
    pub linearize: Symbol,
    pub loop_mom: Symbol,
    pub edgeid: Symbol,
    pub uvaind: Symbol,
    pub edgeaind: Symbol,
    pub dummyaind: Symbol,
    pub hedgeaind: Symbol,
    pub vertexaind: Symbol,
    pub vertexid: Symbol,
    pub source_id: Symbol,
    pub sink_id: Symbol,
    pub orientation_delta: Symbol,
    pub ubar: Symbol,
    pub hfunction_lu_cut: Symbol,
    pub hfunction_left_th: Symbol,
    pub hfunction_right_th: Symbol,
    pub eta: Symbol,
    pub eta_left: Symbol,
    pub eta_right: Symbol,
    pub deta_lu_cut: Symbol,
    pub deta_left_th: Symbol,
    pub deta_right_th: Symbol,
    pub vbar: Symbol,
    pub ose: Symbol,
    /// Positive on-shell energy with no additive energy shift. Arguments carry
    /// the component owner and the invariant mass squared plus spatial norm squared.
    pub on_shell_energy: Symbol,
    /// A CFF energy sum, including its shift, with every dependency in its arguments.
    pub esurface: Symbol,
    pub energy: Symbol,
    pub v: Symbol,
    pub u: Symbol,
    pub color_wrap: Symbol,
    /// Epsilon for dimensional regularization
    pub dim_epsilon: Symbol,

    pub epsilon: Symbol,
    pub epsilonbar: Symbol,
    pub rescale: Symbol,
    /// Bookkeeping scale that restores consumed loop measures in integrated projections.
    pub integrated_loop_scale: Symbol,
    pub rescale_mass: Symbol,
    pub rescale_star: Symbol,
    pub pi: Symbol,

    //Parameters for UV renormalization and localization
    /// UV mass used only as the auxiliary expansion/deformation scale in local UV series.
    ///
    /// This symbol prints as `mUVexp`. Public runtime output should not depend on this
    /// name after the endpoint mass-role collapse.
    pub m_uv_expansion: Symbol,
    /// UV vacuum mass used by integrated counterterms and the vacuum-integral basis.
    ///
    /// This symbol prints as `mUV` and is the public symbolic representative of the
    /// `general.m_uv` runtime setting.
    pub m_uv_vacuum: Symbol,
    /// UV localization scale factor
    pub renormalization_localization_scale: Symbol,
    pub mu_r_sq: Symbol,
    pub numerator_sampling_scale: Symbol,
    pub sign: Symbol,
    pub theta: Symbol,
    pub broadcasting_sqrt: Symbol,
    pub heaviside: Symbol,
    ///for selecting orientations at generation
    pub selected: Symbol,

    /// For analytic UV profiling
    pub expansion: Symbol,
    ///For selecting a concete index.
    pub delta_vec: Symbol,
    /// Canonical spatial component `Q(edge, index___)`.
    pub emr_mom: Symbol,
    /// Component-local on-shell point: boundary owner, signed physical mass,
    /// four fixed velocity components, and optionally a vector index.
    /// The mass remains visible to outer jets.
    pub(crate) os_mom: Symbol,
    /// UV-local provenance stored in the tag slot of `Q(...)` while applying a
    /// Taylor operator. Its arguments are the immutable source edge, its
    /// `UvMomentumProvenanceRole`, and the literal momentum payload. Hard roles
    /// retain their complete projection in the frozen child LMB; newly derived
    /// soft factors retain their crown carrier until outer-CFF routing.
    pub uv_momentum_provenance: Symbol,
    /// Projection-owned denominator class, never a physical edge identifier.
    pub uv_class: Symbol,
    /// Raw `Q3` symbol for pattern construction; stored atoms always carry a
    /// representation or `cind` through the canonical constructor.
    pub emr_vec: Symbol,
    pub dot: Symbol,
    pub external_mom: Symbol,
    pub dim: Symbol,
    pub coeff: Symbol,

    pub localizing_integrand: Symbol,

    pub uv_subgraph: Symbol,
    pub uv_approx: Symbol,
    pub uv_integrate: Symbol,
    pub uv_series: Symbol,
    pub uv_truncate: Symbol,
    /// Marker for the ct term in the UV integrand
    pub ct_marker: Symbol,

    pub nc2_1: Symbol,
    pub expr: Symbol,
    pub num: Symbol,
    /// Denominator wrapper `den(edge, momentum, mass, full_expr)` (no power and
    /// should not be multiplied in but divided).
    pub den: Symbol,
    pub radius_left: Symbol,
    pub radius_star_left: Symbol,
    pub uv_damp_plus_left: Symbol,
    pub uv_damp_minus_left: Symbol,
    pub radius_right: Symbol,
    pub radius_star_right: Symbol,
    pub uv_damp_plus_right: Symbol,
    pub uv_damp_minus_right: Symbol,
    pub thermal_distribution: Symbol,
    pub thermal_weight_wrapper: Symbol,
    pub inverse_temperature: Symbol,
}

impl GammaloopSymbols {
    pub fn collect_orientation_if<'a>(&self, arg: impl Into<AtomOrView<'a>>) -> Atom {
        let arg = arg.into();
        if !arg.as_view().contains_symbol(self.theta) && !arg.as_view().contains_symbol(Symbol::IF)
        {
            return arg.into_owned();
        }
        let mut expression = arg
            .replace(self.sign_theta(W_.a_))
            .with(Symbol::IF.call(Atom::var(W_.a_) + 1))
            .replace(Symbol::IF.call(W_.a_) * Symbol::IF.call(W_.b_))
            .repeat()
            .with(Symbol::IF.call(W_.a_ * W_.b_))
            .replace(Symbol::IF.call(W_.a_) * W_.b___)
            .with(Symbol::IF.call_args([Atom::var(W_.a_), Atom::var(W_.b___), Atom::Zero]))
            .replace(Symbol::IF.call_args([Atom::var(W_.a_), Atom::Zero]))
            .with(Symbol::IF.call_args([Atom::var(W_.a_), Atom::one(), Atom::Zero]))
            .replace(Symbol::IF.call_args([Atom::var(W_.a_)]))
            .with(Symbol::IF.call_args([Atom::var(W_.a_), Atom::one(), Atom::Zero]));
        loop {
            // A zero-arm selector owns its complete scalar weight, including
            // inverse factors which must remain inside the lazy selected arm.
            // Function arguments can be literal alias tags: absorb weights only
            // in arithmetic and IF bodies, without changing those lookup keys.
            let weighted = expression.replace_map(|node, _, output| {
                if let AtomView::Fun(call) = node
                    && call.get_symbol() != Symbol::IF
                {
                    output.set_from_view(&node);
                    return;
                }
                let AtomView::Mul(product) = node else { return };
                let Some(condition) = product.iter().find_map(|factor| {
                    if let AtomView::Fun(branch) = factor
                        && branch.get_symbol() == Symbol::IF
                        && branch.get_nargs() == 3
                        && branch.get(1).is_zero()
                    {
                        Some(branch.get(0))
                    } else {
                        None
                    }
                }) else {
                    return;
                };
                let body = Atom::mul_many(product.iter().map(|factor| {
                    if let AtomView::Fun(branch) = factor
                        && branch.get_symbol() == Symbol::IF
                        && branch.get_nargs() == 3
                        && branch.get(0) == condition
                        && branch.get(1).is_zero()
                    {
                        branch.get(2)
                    } else {
                        factor
                    }
                }));
                **output = Symbol::IF.call_args([condition, Atom::Zero.as_view(), body.as_view()]);
            });
            // Ordinary same-key sums remain collectable inside functions. Their
            // normalizers may create new sums, so both phases share this fixed
            // point rather than rescanning the expression for each branch pair.
            let collected = weighted.replace_map_bottom_up(|node, _, output| {
                if let Some(collected) = Self::merge_orientation_branches(node) {
                    **output = collected;
                }
            });
            if collected == expression {
                return collected;
            }
            // Rebuilding a parent can invoke a normalizer that introduces new
            // conditional sums after its children have already been visited.
            expression = collected;
        }
    }

    fn merge_orientation_branches(node: AtomView<'_>) -> Option<Atom> {
        let AtomView::Add(sum) = node else {
            return None;
        };
        let mut groups = BTreeMap::<_, Vec<_>>::new();
        let mut terms = Vec::<AtomOrView<'_>>::new();
        for term in sum.iter() {
            if let AtomView::Fun(branch) = term
                && branch.get_symbol() == Symbol::IF
                && branch.get_nargs() == 3
                && branch.get(1).is_zero()
            {
                groups
                    .entry(branch.get(0))
                    .or_default()
                    .push((term, branch.get(2)));
            } else {
                terms.push(term.into());
            }
        }
        if groups.values().all(|group| group.len() == 1) {
            return None;
        }
        // Group direct scalar branches after weight absorption. Powers and
        // unrelated factors stay intact. One bulk sum avoids repeated matching
        // and rebuilding a growing body.
        for (condition, group) in groups {
            if group.len() == 1 {
                terms.push(group[0].0.into());
                continue;
            }
            let body = Atom::add_many(group.into_iter().map(|(_, body)| body));
            // Joining bodies can expose another sum of conditionals. Its
            // children were already visited; only collect this new sum root.
            let body = Self::merge_orientation_branches(body.as_view()).unwrap_or(body);
            terms.push(
                Symbol::IF
                    .call_args([condition, Atom::Zero.as_view(), body.as_view()])
                    .into(),
            );
        }
        Some(Atom::add_many(terms))
    }

    pub fn den<'a>(
        &self,
        eid: impl Into<AtomOrView<'a>>,
        mom: impl Into<AtomOrView<'a>>,
        mass: impl Into<AtomOrView<'a>>,
        full_expr: impl Into<AtomOrView<'a>>,
    ) -> Atom {
        self.den.call_args([
            eid.into().as_view(),
            mom.into().as_view(),
            mass.into().as_view(),
            full_expr.into().as_view(),
        ])
    }

    pub(crate) fn orientation_delta<O: GraphOrientation>(&self, orientation: &O) -> Atom {
        let args: Vec<i32> = orientation
            .orientation()
            .iter()
            .map(|(_, t)| match t {
                Orientation::Default => 1,
                Orientation::Reversed => -1,
                Orientation::Undirected => 0,
            })
            .collect_vec();
        FunctionBuilder::new(self.orientation_delta)
            .add_args(&args)
            .finish()
    }

    pub(crate) fn sign_theta<'a>(&self, arg: impl Into<AtomOrView<'a>>) -> Atom {
        let arg = arg.into();

        function!(self.theta, arg.as_view())
    }

    pub(crate) fn sign(&self, edge: EdgeIndex) -> Atom {
        function!(self.sign, Atom::num(edge.0 as i64))
    }
}

pub static W_: LazyLock<WildCards> = LazyLock::new(|| WildCards {
    edgeid_: symbol!("eid_"),
    mom_: symbol!("mom_"),
    mass_: symbol!("mass_"),
    prop_: symbol!("prop_"),
    a_: symbol!("a_"),
    b_: symbol!("b_"),
    c_: symbol!("c_"),
    d_: symbol!("d_"),
    e_: symbol!("e_"),
    f_: symbol!("f_"),
    g_: symbol!("g_"),
    h_: symbol!("h_"),
    i_: symbol!("i_"),
    j_: symbol!("j_"),
    k_: symbol!("k_"),
    l_: symbol!("l_"),
    m_: symbol!("m_"),
    n_: symbol!("n_"),
    o_: symbol!("o_"),
    p_: symbol!("p_"),
    q_: symbol!("q_"),
    r_: symbol!("r_"),
    s_: symbol!("s_"),
    t_: symbol!("t_"),
    u_: symbol!("u_"),
    v_: symbol!("v_"),
    w_: symbol!("w_"),
    x_: symbol!("x_"),
    y_: symbol!("y_"),
    z_: symbol!("z_"),
    a__: symbol!("a__"),
    b__: symbol!("b__"),
    c__: symbol!("c__"),
    d__: symbol!("d__"),
    e__: symbol!("e__"),
    f__: symbol!("f__"),
    g__: symbol!("g__"),
    h__: symbol!("h__"),
    i__: symbol!("i__"),
    j__: symbol!("j__"),
    k__: symbol!("k__"),
    l__: symbol!("l__"),
    m__: symbol!("m__"),
    n__: symbol!("n__"),
    o__: symbol!("o__"),
    p__: symbol!("p__"),
    q__: symbol!("q__"),
    r__: symbol!("r__"),
    s__: symbol!("s__"),
    t__: symbol!("t__"),
    u__: symbol!("u__"),
    v__: symbol!("v__"),
    w__: symbol!("w__"),
    x__: symbol!("x__"),
    y__: symbol!("y__"),
    z__: symbol!("z__"),
    a___: symbol!("a___"),
    b___: symbol!("b___"),
    c___: symbol!("c___"),
    d___: symbol!("d___"),
    e___: symbol!("e___"),
    f___: symbol!("f___"),
    g___: symbol!("g___"),
    h___: symbol!("h___"),
    i___: symbol!("i___"),
    j___: symbol!("j___"),
    k___: symbol!("k___"),
    l___: symbol!("l___"),
    m___: symbol!("m___"),
    n___: symbol!("n___"),
    o___: symbol!("o___"),
    p___: symbol!("p___"),
    q___: symbol!("q___"),
    r___: symbol!("r___"),
    s___: symbol!("s___"),
    t___: symbol!("t___"),
    u___: symbol!("u___"),
    v___: symbol!("v___"),
    w___: symbol!("w___"),
    x___: symbol!("x___"),
    y___: symbol!("y___"),
    z___: symbol!("z___"),
});

macro_rules! spenso_print_uv_unary {
    ($a:ident, $opt:ident, $prefix:expr, $suffix:expr) => {{
        match $opt.custom_print_mode.get("spenso") {
            Some(PrintUserData::Integer(_)) => {
                let AtomView::Fun(f) = $a else {
                    return None;
                };
                if f.get_nargs() != 1 {
                    return None;
                }

                let mut out = $prefix.to_string();
                f.iter()
                    .next()
                    .unwrap()
                    .format(&mut out, $opt, PrintState::new())
                    .unwrap();
                out.push_str($suffix);
                Some(out)
            }
            _ => None,
        }
    }};
}

spenso::symbolica_init_lazy_static! {
pub static GS, GS_INNER: GammaloopSymbols = || GammaloopSymbols {
    renormalization_localization_scale: symbol!("rls"),
    integrand: symbol!("integrand"),
    residue_map_id: symbol!("residue_map_id"),
    tree_denom_wrapper: symbol!("tree_denoms"),
    dim_epsilon: symbol!("ε"),
    killing_func: symbol!(
        "killing_func",
        norm = |f, out| {
            if let AtomView::Fun(_f) = f {
                **out = Atom::one()
            }
        }
    ),
    ufozero: symbol!(
        "UFO::ZERO",
        norm = |_, out| {
            **out = Atom::Zero;
        }
    ),

    localizing_integrand: symbol!("int_loc"),
    uvaind: symbol!(
        "uvind",
        print = spenso::network::tags::tensor_print,
        tags = [SPENSO_TAG.index.clone(), "spenso::index-label:u".to_owned()]
    ),
    edgeaind: feynkit_graph::symbols::edge_index(),
    vertexaind: feynkit_graph::symbols::vertex_index(),
    dummyaind: symbol!(
        "dummy",
        print = spenso::network::tags::tensor_print,
        tags = [SPENSO_TAG.index.clone(), "spenso::index-label:d".to_owned()]
    ),
    hedgeaind: feynkit_graph::symbols::hedge_index(),
    uv_subgraph: symbol!(
        "gammalooprs::uv::subgraph",
        print = |a, opt, _state| {
            match opt.custom_print_mode.get("spenso") {
                Some(PrintUserData::Integer(_)) => {
                    let AtomView::Fun(f) = a else {
                        return None;
                    };
                    if f.get_nargs() != 2 {
                        return None;
                    }

                    let mut args = f.iter();
                    let current = args.next().unwrap();
                    let mut out = String::new();
                    if let AtomView::Var(current) = current
                        && current.get_symbol().get_stripped_name().starts_with("S_")
                    {
                        out.push_str(current.get_symbol().get_stripped_name());
                    } else {
                        out.push_str("S_");
                        current
                            .format(&mut out, opt, PrintState::new())
                            .unwrap();
                    }
                    out.push('⊛');
                    let given = args.next().unwrap();
                    if let AtomView::Var(given) = given
                        && let Some(label) =
                            given.get_symbol().get_stripped_name().strip_prefix("S_")
                    {
                        out.push_str(label);
                    } else {
                        given.format(&mut out, opt, PrintState::new()).unwrap();
                    }
                    Some(out)
                }
                _ => None,
            }
        }
    ),
    uv_approx: symbol!(
        "Approx",
        print = |a, opt, _state| spenso_print_uv_unary!(a, opt, "K[", "]")
    ),
    uv_integrate: symbol!(
        "Integrate",
        print = |a, opt, _state| spenso_print_uv_unary!(a, opt, "⟨", "⟩")
    ),
    uv_series: symbol!(
        "Series",
        print = |a, opt, _state| spenso_print_uv_unary!(a, opt, "Σ(", ")")
    ),
    uv_truncate: symbol!(
        "Truncate",
        print = |a, opt, _state| {
            let prefix = if opt.typst_mode().is_some() {
                r#"op("Tr")("#
            } else {
                "Tr("
            };
            spenso_print_uv_unary!(a, opt, prefix, ")")
        }
    ),
    ct_marker: symbol!(
        "CT",
        print = |a, opt, _state| {
            spenso_print_uv_unary!(a, opt, "", "").map(|out| {
                if opt.color_builtin_symbols {
                    nu_ansi_term::Color::Magenta.paint(out).to_string()
                } else {
                    out
                }
            })
        }
    ),
    is_function: symbol!("is_function"),
    is_symbol: symbol!("is_symbol"),
    nc2_1: symbol!("NC2_1"),
    rescale: symbol!("t";Scalar,Real,Positive),
    integrated_loop_scale: symbol!("uvIntegratedLoopScale";Scalar,Real,Positive),
    rescale_mass: symbol!("t_m";Scalar,Real,Positive),
    _linear: symbol!("_linear";Linear),
    linearize: symbol!(
        "linearize",
        norm = |f, out| {
            if let AtomView::Fun(ff) = f
                && ff.get_nargs() == 1
                && let AtomView::Fun(arg) = ff.iter().next().unwrap()
            {
                let mut args = vec![];
                for a in arg.iter() {
                    if let AtomView::Fun(a) = a
                        && a.get_symbol().get_wildcard_level() > 0
                    {
                        return;
                    }
                    if let AtomView::Var(a) = a
                        && a.get_symbol().get_wildcard_level() > 0
                    {
                        return;
                    }

                    args.push(a);
                }
                let expr = GS
                    ._linear
                    .call_args(args.as_slice())
                    .replace(GS._linear.call_args([W_.a___]))
                    .with(arg.get_symbol().call_args([W_.a___]));
                println!(":{}->{}", f, expr);
                **out = expr;
            }
        }
    ),
    broadcasting_sqrt: symbol!(
        "broadcasting_sqrt",
        tag = SPENSO_TAG.broadcast,
        der = |a, _, out| {
            **out = Atom::num(1) / (Atom::num(2) * a);
        }
    ),
    heaviside: symbol!("heaviside"),
    expansion: symbol!("expansion"),
    rescale_star: symbol!("t⃰"),
    hfunction_lu_cut: symbol!("h_lu_cut"),
    hfunction_left_th: symbol!("h_left_th"),
    hfunction_right_th: symbol!("h_right_th"),
    eta: symbol!("η"),
    eta_left: symbol!("η_left"),
    eta_right: symbol!("η_right"),
    deta_lu_cut: symbol!("∇η_lu_cut"),
    deta_left_th: symbol!("∇η_left_th"),
    deta_right_th: symbol!("∇η_right_th"),
    edgeid: symbol!("eid"),
    vertexid: symbol!("vid"),
    source_id: symbol!("source"),
    sink_id: symbol!("sink"),
    sign: symbol!("σ"; Scalar),
    selected: symbol!("selected"),
    theta: symbol!("θ"),
    m_uv_expansion: symbol!(
        "mUVexp",
        print = |a, opt, _state| {
            let AtomView::Var(_a) = a else {
                return None;
            };
            match opt.custom_print_mode.get("spenso") {
                Some(PrintUserData::Integer(_)) if opt.typst_mode().is_some() => {
                    Some("m_\"UVexp\"".to_string())
                }
                _ => None,
            }
        }
    ),
    m_uv_vacuum: symbol!(
        "mUV",
        print = |a, opt, _state| {
            let AtomView::Var(_a) = a else {
                return None;
            };
            match opt.custom_print_mode.get("spenso") {
                Some(PrintUserData::Integer(_)) if opt.typst_mode().is_some() => {
                    Some("m_\"UV\"".to_string())
                }
                _ => None,
            }
        }
    ),
    mu_r_sq: symbol!(
        "μᵣ²",
        print = |a, opt, _state| {
            let AtomView::Var(_a) = a else {
                return None;
            };
            match opt.custom_print_mode.get("spenso") {
                Some(PrintUserData::Integer(_)) if opt.typst_mode().is_some() => {
                    Some("mu_r^2".to_string())
                }
                _ => None,
            }
        }
    ),
    numerator_sampling_scale: symbol!("M"),
    delta_vec: ETS.delta,
    expr: symbol!(
        "expr",
        print = |_, opt, _state| {
            Some(match opt.typst_mode()? {
                TypstMode::Math | TypstMode::Markup => "#expr".to_string(),
                TypstMode::Code => "expr".to_string(),
            })
        }
    ),
    num: symbol!("num"),
    den: feynkit_graph::symbols::denominator(),
    ubar: feynkit_graph::symbols::ubar(),
    vbar: feynkit_graph::symbols::vbar(),
    dot: symbol!("dot"),
    dim: feynkit_graph::symbols::dimension(),
    v: feynkit_graph::symbols::v(),
    u: feynkit_graph::symbols::u(),
    emr_mom: feynkit_graph::symbols::momentum(),
    // The original index-free momentum is provenance metadata, not a tensor child.
    uv_momentum_provenance: symbol!("gammalooprs::uv::momentum_provenance"; Scalar),
    os_mom: symbol!(
        "gammalooprs::uv::on_shell_point",
        tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
    ),
    uv_class: symbol!("gammalooprs::uv::class"),
    orientation_delta: symbol!("orientation_delta"),
    emr_vec: symbol!(
        "Q3",
        norm = |f, out| {
            if let AtomView::Fun(ff) = f
                && ff.get_nargs() == 2
            {
                let mut iter = ff.iter();
                let eid = iter.next().unwrap();
                let structure = iter.next().unwrap();
                let component = if let AtomView::Fun(cind) = structure
                    && cind.get_symbol() == AIND_SYMBOLS.cind
                    && cind.get_nargs() == 1
                {
                    cind.iter()
                        .next()
                        .and_then(|i| i64::try_from(i).ok())
                        .and_then(|i| usize::try_from(i).ok())
                } else {
                    None
                };
                if i64::try_from(eid).is_ok_and(|eid| eid >= 0)
                    && let Some(i @ 0..=3) = component
                {
                    if i == 0 {
                        **out = Atom::Zero;
                    } else {
                        let cind = AIND_SYMBOLS.cind.call_args([i as i64]);
                        **out = get_symbol!("Q")
                            .unwrap()
                            .call_args([eid, cind.as_view()])
                    }
                }
            }
        },
        print = |a, opt, _state| {
            spenso_print_scripted_indexed!(a, opt, "q⃗", "arrow(q)")
        },
        tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
    ),
    ose: feynkit_cff::symbols::on_shell(),
    on_shell_energy: symbol!(
        "OnShellEnergy"; Scalar;
        norm = |view, out| {
            if let AtomView::Fun(call) = view && call.get_nargs() == 2 {
                // Series reconstruction repeatedly normalizes identical energy
                // invariants. Keep a bounded cache per worker, with no borrow
                // held during symbolic operations or recursive normalization.
                thread_local! {
                    static CACHE: RefCell<(HashMap<Atom, Atom>, usize)> = RefCell::default();
                }
                let input = call.get(1);
                let cached = CACHE.with_borrow(|(cache, _)| cache.get(input.get_data()).cloned());
                // Canonicalize opposite routed momenta inside squared factors.
                // Distribute numeric coefficients only; keep the invariant's
                // momentum sums squared and graph numerators outside this call.
                let invariant = cached.unwrap_or_else(|| {
                    let invariant = input.replace_map(|part, _, output| {
                        if let AtomView::Pow(power) = part {
                            let (base, exponent) = power.get_base_exp();
                            if i64::try_from(exponent) == Ok(2) {
                                let direct = base.expand_num();
                                let reversed = (-&direct).expand_num();
                                **output = direct.min(reversed).pow(2);
                            }
                        }
                    }).expand_num();
                    let bytes = input.get_byte_size() + invariant.as_view().get_byte_size();
                    const MAX_BYTES: usize = 1024 * 1024;
                    if bytes <= MAX_BYTES {
                        CACHE.with_borrow_mut(|(cache, retained_bytes)| {
                            if cache.len() >= 1024 || *retained_bytes + bytes > MAX_BYTES {
                                cache.clear();
                                *retained_bytes = 0;
                            }
                            cache.insert(input.to_owned(), invariant.clone());
                            *retained_bytes += bytes;
                        });
                    }
                    invariant
                });
                if invariant.as_view() != input {
                    **out = call.get_symbol().call_args([call.get(0).to_owned(), invariant]);
                }
            }
        },
        der = |view, argument, out| {
            if let AtomView::Fun(call) = view && call.get_nargs() == 2 {
                **out = if argument == 1 {
                    Atom::num((1, 2)) / view
                } else {
                    Atom::Zero
                };
            }
        }
    ),
    esurface: symbol!(
        "Esurface"; Scalar;
        norm = |view, out| {
            if let AtomView::Fun(call) = view {
                let sum = Atom::add_many(call.iter());
                if let AtomView::Add(terms) = sum.as_view() {
                    // Pull explicit Laurent monomials outside the function
                    // before series expansion. Otherwise a function callback
                    // receives already truncated pole-valued arguments and can
                    // lose absolute precision when their leading terms cancel.
                    let mut poles = std::collections::BTreeMap::<Symbol, i64>::new();
                    let mut record = |factor: AtomView<'_>| {
                        if let AtomView::Pow(power) = factor
                            && let AtomView::Var(variable) = power.get_base()
                            && let Ok(exponent) = i64::try_from(power.get_exp())
                            && exponent < 0
                        {
                            poles.entry(variable.get_symbol())
                                .and_modify(|power| *power = (*power).min(exponent))
                                .or_insert(exponent);
                        }
                    };
                    for term in terms.iter() {
                        if let AtomView::Mul(product) = term {
                            product.iter().for_each(&mut record);
                        } else {
                            record(term);
                        }
                    }
                    if let Some((variable, exponent)) = poles.first_key_value() {
                        let factor = Atom::var(*variable).pow(*exponent);
                        let regular = call.get_symbol().call_args(
                            terms.iter().map(|term| term / factor.as_view())
                        );
                        **out = factor * regular;
                        return;
                    }
                    if !call.iter().eq(terms.iter()) {
                        **out = call.get_symbol().call_args(terms.iter());
                    }
                } else {
                    **out = sum;
                }
            }
        },
        der = |_, _, out| { **out = Atom::one(); },
        eval = symbolica::atom::EvaluationInfo::new()
            .register(|arguments: &[f64]| arguments.iter().sum::<f64>())
            .register(|arguments: &[symbolica::domains::float::Complex<symbolica::domains::float::Float>]| {
                arguments.iter().cloned().reduce(|left, right| left + right)
                    .expect("a retained CFF surface has at least two summands")
            })
    ),
    energy: feynkit_cff::symbols::energy(),
    external_mom: feynkit_graph::symbols::external_momentum(),
    loop_mom: feynkit_graph::symbols::loop_momentum(),
    epsilon: feynkit_graph::symbols::epsilon(),
    pi: Symbol::PI,
    color_wrap: symbol!("color"),
    epsilonbar: feynkit_graph::symbols::epsilonbar(),
    coeff: symbol!("coef"),
    radius_left: symbol!("r_left"),
    radius_star_left: symbol!("r⃰_left"),
    uv_damp_plus_left: symbol!("damp_plus_left"),
    uv_damp_minus_left: symbol!("damp_minus_left"),
    radius_right: symbol!("r_right"),
    radius_star_right: symbol!("r⃰_right"),
    uv_damp_plus_right: symbol!("damp_plus_right"),
    uv_damp_minus_right: symbol!("damp_minus_right"),
    thermal_distribution: three_dimensional_reps::symbols::thermal_distribution_symbol(),
    thermal_weight_wrapper: *three_dimensional_reps::symbols::THERMAL_WEIGHT_WRAPPER,
    inverse_temperature: symbol!("β"),
};
}

impl GammaloopSymbols {
    pub fn wrap_tree_denoms<'a>(&self, arg: impl Into<AtomOrView<'a>>) -> Atom {
        self.tree_denom_wrapper.call_args([arg.into().as_view()])
    }

    pub fn do_dot_product_in_sqrt<'a>(&self, arg: impl Into<AtomOrView<'a>>) -> Atom {
        let a = arg.into();

        let mink = Minkowski {}.new_rep(4).to_symbolic([Atom::var(W_.i_)]);

        let pat = function!(self.emr_vec, W_.a_, mink) * function!(self.emr_vec, W_.b_, mink);
        let rhs = function!(
            self.emr_mom,
            W_.a_,
            Atom::from(ExpandedIndex::from_iter([1]))
        ) * function!(
            self.emr_mom,
            W_.b_,
            Atom::from(ExpandedIndex::from_iter([1]))
        ) + function!(
            self.emr_mom,
            W_.a_,
            Atom::from(ExpandedIndex::from_iter([2]))
        ) * function!(
            self.emr_mom,
            W_.b_,
            Atom::from(ExpandedIndex::from_iter([2]))
        ) + function!(
            self.emr_mom,
            W_.a_,
            Atom::from(ExpandedIndex::from_iter([3]))
        ) * function!(
            self.emr_mom,
            W_.b_,
            Atom::from(ExpandedIndex::from_iter([3]))
        );

        let a = a.replace(pat).with(rhs);

        let pat = function!(self.emr_vec, W_.a_, mink).pow(2);
        let rhs = function!(
            self.emr_mom,
            W_.a_,
            Atom::from(ExpandedIndex::from_iter([1]))
        )
        .pow(2)
            + function!(
                self.emr_mom,
                W_.a_,
                Atom::from(ExpandedIndex::from_iter([2]))
            )
            .pow(2)
            + function!(
                self.emr_mom,
                W_.a_,
                Atom::from(ExpandedIndex::from_iter([3]))
            )
            .pow(2);
        // * GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([1])))
        // + GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([2])))
        //     * GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([2])))
        // + GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([3])))
        //     * GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([3])));

        // let dot = self.emr_vec(e, mink.to_atom()) * self.emr_vec(e, mink.to_atom())
        a.replace(pat).with(rhs)

        // a.replace_map(|a, ctx, out| {
        //     if let AtomView::Pow(p) = a {
        //         let (a, b) = p.get_base_exp();
        //         let Ok(exp) = Rational::try_from(b) else {
        //             return;
        //         };

        //         if exp.denominator() == Integer::from(2) {
        //             **out = self
        //                 .broadcasting_sqrt
        //                 .call_args([a])
        //                 .pow(Atom::num(exp.numerator()));
        //         }
        //     }
        // });
    }

    pub fn to_broadcasting_sqrt<'a>(&self, arg: impl Into<AtomOrView<'a>>) -> Atom {
        let a = arg.into();
        a.replace_map(|a, _ctx, out| {
            if let AtomView::Pow(p) = a {
                let (a, b) = p.get_base_exp();
                let Ok(exp) = Rational::try_from(b) else {
                    return;
                };

                if exp.denominator() == 2 {
                    **out = self
                        .broadcasting_sqrt
                        .call_args([a])
                        .pow(Atom::num(exp.numerator()));
                }
            }
        })
    }

    pub(crate) fn heaviside<'a>(&self, arg: impl Into<AtomOrView<'a>>) -> Atom {
        let a = arg.into();
        function!(self.heaviside, a.as_view())
    }

    pub(crate) fn emr_mom<'a>(&self, e: EdgeIndex, arg: impl Into<AtomOrView<'a>>) -> Atom {
        let a = arg.into();
        function!(self.emr_mom, usize::from(e) as i64, a.as_view())
    }

    pub(crate) fn uv_momentum_provenance_tag<'a>(
        &self,
        edge: impl Into<AtomOrView<'a>>,
        role: impl Into<UvMomentumProvenanceRole>,
        momentum: impl Into<AtomOrView<'a>>,
    ) -> Atom {
        let role = role.into();
        self.uv_momentum_provenance.call_args([
            edge.into().as_view(),
            Atom::num(role as i64).as_view(),
            momentum.into().as_view(),
        ])
    }

    pub(crate) fn uv_momentum_provenance_data(
        &self,
        argument: AtomView<'_>,
    ) -> Option<(EdgeIndex, UvMomentumProvenanceRole, Atom)> {
        let AtomView::Fun(provenance) = argument else {
            return None;
        };
        if provenance.get_symbol() != self.uv_momentum_provenance || provenance.get_nargs() != 3 {
            return None;
        }
        let edge = EdgeIndex(usize::try_from(provenance.get(0)).ok()?);
        let role = match i64::try_from(provenance.get(1)).ok()? {
            0 => UvMomentumProvenanceRole::TaylorFixed,
            1 => UvMomentumProvenanceRole::DenominatorDerived,
            3 => UvMomentumProvenanceRole::PhysicalSourceFixed,
            4 => UvMomentumProvenanceRole::DenominatorDerivedSoft,
            _ => return None,
        };
        Some((edge, role, provenance.get(2).to_owned()))
    }

    pub(crate) fn uv_class_ref(&self, class: UvDenominatorClassId) -> Atom {
        self.uv_class.call(class.0)
    }

    pub(crate) fn uv_class_data(&self, argument: AtomView<'_>) -> Option<UvDenominatorClassId> {
        let AtomView::Fun(class) = argument else {
            return None;
        };
        (class.get_symbol() == self.uv_class && class.get_nargs() == 1)
            .then(|| usize::try_from(class.get(0)).ok().map(UvDenominatorClassId))
            .flatten()
    }

    /// Attach the same component indices to every vector in a linear momentum.
    pub(crate) fn indexed_momentum(&self, momentum: &Atom, indices: &[Atom]) -> Atom {
        momentum.replace_map(|view, _, output| {
            let AtomView::Fun(vector) = view else {
                return;
            };
            let arity = if vector.get_symbol() == self.emr_mom {
                1
            } else if vector.get_symbol() == self.os_mom {
                6
            } else {
                return;
            };
            if vector.get_nargs() == arity {
                let mut component = FunctionBuilder::new(vector.get_symbol());
                for argument in vector.iter() {
                    component = component.add_arg(argument);
                }
                for index in indices {
                    component = component.add_arg(index);
                }
                **output = component.finish();
            }
        })
    }

    /// A fixed component-owned point. The signed physical mass is kept as one
    /// scalar argument so enclosing Taylor operations also scale its four-vector.
    pub(crate) fn on_shell_momentum(
        &self,
        owner: EdgeIndex,
        signed_mass: &Atom,
        velocity: &[Atom; 4],
    ) -> Atom {
        FunctionBuilder::new(self.os_mom)
            .add_arg(usize::from(owner))
            .add_arg(signed_mass)
            .add_args(velocity)
            .finish()
    }

    /// Lower indexed on-shell points, retaining bare vector metadata for
    /// exact denominator reconstruction.
    pub(crate) fn resolve_os_momenta(&self, atom: &Atom) -> Atom {
        atom.replace_map(|view, _, output| {
            if let AtomView::Fun(point) = view
                && point.get_symbol() == self.os_mom
                && point.get_nargs() == 7
            {
                let vector = Atom::add_many((0..4).map(|component| {
                    point.get(component + 2).to_owned()
                        * self
                            .delta_vec
                            .call_args([self.cind(component), point.get(6).to_owned()])
                }));
                **output = point.get(1) * vector.as_view();
            }
        })
    }

    /// Remove only the UV Taylor provenance from EMR tags. This is used at
    /// algebraic boundaries such as Vakint and by exact test oracles; the local
    /// 4D-to-CFF route deliberately retains the tag until energy ownership has
    /// been planned.
    pub(crate) fn erase_uv_momentum_provenance(&self, atom: &Atom) -> Atom {
        atom.replace_map(|view, _, output| {
            let AtomView::Fun(momentum) = view else {
                return;
            };
            if momentum.get_symbol() != self.emr_mom || momentum.get_nargs() < 1 {
                return;
            }
            let AtomView::Fun(provenance) = momentum.get(0) else {
                return;
            };
            if provenance.get_symbol() != self.uv_momentum_provenance || provenance.get_nargs() != 3
            {
                return;
            }
            let indices = momentum
                .iter()
                .skip(1)
                .map(|index| index.to_owned())
                .collect::<Vec<_>>();
            **output = self.indexed_momentum(&provenance.get(2).to_owned(), &indices);
        })
    }

    pub(crate) fn localizing_integrand(&self, lmb: &LoopMomentumBasis) -> Atom {
        // Normalize each factor with int d^3k / (|k|^2 + rls^2)^2 = pi^2 / rls,
        // so the localizing integrand itself integrates to one.
        let pi_atom = (Symbol::PI).to_atom();
        let normalization_term_integral = (pi_atom.pow(2)) / GS.renormalization_localization_scale;

        let mut res = Atom::one();

        for l in lmb.loop_edges.iter() {
            // Orientation ownership is applied by `Localizer`: orientation-local 3D uses one
            // deterministic valid host, while selector-free 3D and local 4D explicitly sum their
            // source orientations. It is therefore not an extra factor in this normalized radial
            // kernel.
            res /= normalization_term_integral.as_view();

            let spatial_norm_sq = function!(self.emr_mom, l.0, GS.cind(1)).pow(2)
                + function!(self.emr_mom, l.0, GS.cind(2)).pow(2)
                + function!(self.emr_mom, l.0, GS.cind(3)).pow(2);

            // Per-orientation CFF localizer of the normalized cubic tadpole.
            let denominator = spatial_norm_sq
                + GS.renormalization_localization_scale * GS.renormalization_localization_scale;
            res /= denominator.as_view() * denominator.as_view();
        }

        if res.is_one() {
            res
        } else {
            function!(self.localizing_integrand, res)
        }
    }

    /// Construct the canonical represented or indexed spatial momentum `Q3(e, structure)`.
    ///
    /// `e` is a non-negative graph edge ID. Stored abstract vectors use a
    /// four-dimensional Minkowski representation. Tensor execution transiently
    /// emits `Q3(e,cind(i))`; the symbol normalizer immediately maps that image
    /// to `Q(e,cind(i))` (or zero for the temporal component), so it is never a
    /// persistent second convention. An indexed Minkowski argument still carries
    /// an abstract tensor index and must not be interpreted as a component. The
    /// trailing argument may instead be one wildcard when constructing a
    /// replacement pattern. This API deliberately has no bare `Q3(e)` form.
    pub(crate) fn emr_vec<'a>(&self, e: EdgeIndex, structure: impl Into<AtomOrView<'a>>) -> Atom {
        function!(
            GS.emr_vec,
            usize::from(e) as i64,
            structure.into().as_view()
        )
    }

    pub(crate) fn cind(&self, e: usize) -> Atom {
        AIND_SYMBOLS.cind.call_args([e as i64])
    }

    pub(crate) fn energy_delta<'a>(&self, index: impl Into<AtomOrView<'a>>) -> Atom {
        self.delta_vec
            .call_args([self.cind(0), index.into().into_owned()])
    }

    /// Keep a surface or one of its Taylor coefficients opaque to scalar
    /// factor algebra, without hiding any physical dependency behind an ID.
    pub(crate) fn wrap_esurface(&self, sum: &Atom) -> Atom {
        if let AtomView::Add(terms) = sum.as_view() {
            self.esurface.call_args(terms.iter())
        } else {
            sum.clone()
        }
    }

    /// Apply only to the denominator-only CFF, before attaching a numerator.
    pub(crate) fn retain_esurfaces(&self, cff: &Atom) -> Atom {
        cff.replace_map(|view, _, out| {
            if let AtomView::Pow(power) = view {
                let (base, exponent) = power.get_base_exp();
                if matches!(base, AtomView::Add(_))
                    && i64::try_from(exponent).is_ok_and(|power| power < 0)
                    && (base.contains_symbol(self.ose) || base.contains_symbol(self.energy))
                {
                    **out = self.wrap_esurface(&base.to_owned()).pow(exponent);
                }
            }
        })
    }

    /// Lower on an evaluator or diagnostic copy, after all forest operations.
    /// This unfolds sums only; it never distributes a numerator product.
    pub(crate) fn lower_esurfaces(&self, atom: &Atom) -> Atom {
        if !atom.contains_symbol(self.esurface) {
            return atom.clone();
        }
        atom.replace_map_bottom_up(|view, _, out| {
            if let AtomView::Fun(call) = view
                && call.get_symbol() == self.esurface
            {
                **out = Atom::add_many(call.iter());
            }
        })
    }

    pub(crate) fn ose_full(
        &self,
        e: EdgeIndex,
        lmb_id: EdgeIndex,
        e_mass: Atom,
        index: Option<Atom>,
    ) -> Atom {
        let ose = function!(
            self.on_shell_energy,
            lmb_id.0,
            e_mass.pow(2)
                + self.emr_vec(e, self.cind(1)).pow(2)
                + self.emr_vec(e, self.cind(2)).pow(2)
                + self.emr_vec(e, self.cind(3)).pow(2)
        );

        if let Some(index) = index {
            ose * self.energy_delta(index)
        } else {
            ose
        }
    }

    pub(crate) fn split_mom_pattern(
        &self,
        e: EdgeIndex,
        lmb_id: EdgeIndex,
        e_mass: Atom,
    ) -> Replacement {
        let id = Minkowski {}.to_symbolic([W_.a__]);
        Replacement::new(
            self.emr_mom(e, &id).to_pattern(),
            self.emr_vec(e, &id) + self.ose_full(e, lmb_id, e_mass, Some(id)) * sign_atom(e),
        )
    }

    /// Split Q(i,mu) into its spatial part and the already-generated on-shell
    /// energy, with the production sign attached to the time component. The
    /// typed 4D CFF source owns mass-dependent on-shell reconstruction before
    /// this late numerator/export boundary. The split is
    /// `Q(i,mink(4,mu)) -> Q3(i,mink(4,mu)) + OSE(i)σ(i)δ(cind(0),mink(4,mu))`.
    pub fn split_mom_pattern_simple(&self, e: EdgeIndex) -> Replacement {
        let eidc = usize::from(e) as i64;
        let index = Minkowski {}.to_symbolic([W_.a__]);
        Replacement::new(
            self.emr_mom(e, &index).to_pattern(),
            function!(GS.ose, eidc) * sign_atom(e) * self.energy_delta(&index)
                + function!(GS.emr_vec, eidc, &index),
        )
    }

    pub(crate) fn ose(&self, e: EdgeIndex) -> Atom {
        function!(GS.ose, usize::from(e) as i64)
    }
}

pub(crate) fn sign_atom(eid: EdgeIndex) -> Atom {
    FunctionBuilder::new(symbol!("σ"))
        .add_arg(usize::from(eid) as i64)
        .finish()
}

#[cfg(test)]
mod tests {
    use insta::assert_snapshot;
    use spenso::shadowing::symbolica_utils::LogPrint;
    use spenso::shadowing::symbolica_utils::SpensoPrintSettings;

    use super::*;

    #[test]
    fn fixed_on_shell_point_scales_with_its_physical_mass_in_every_component() {
        crate::initialisation::test_initialise().unwrap();
        let mass = symbol!("on_shell_point_scaling::mass"; Scalar);
        let parameter = symbol!("on_shell_point_scaling::parameter"; Scalar);
        let mass_atom = Atom::var(mass);
        let scale = Atom::var(parameter);
        // A unit timelike velocity with every spatial component nonzero.
        let velocity = [
            Atom::num((13, 8)),
            Atom::num(1),
            Atom::num((1, 2)),
            Atom::num((5, 8)),
        ];
        for sign in [1, -1] {
            let point =
                GS.on_shell_momentum(EdgeIndex(7), &(Atom::num(sign) * &mass_atom), &velocity);
            let scaled_point = point.replace(mass).with(&scale * &mass_atom);
            assert_eq!(
                scaled_point,
                GS.on_shell_momentum(
                    EdgeIndex(7),
                    &(Atom::num(sign) * &scale * &mass_atom),
                    &velocity
                ),
            );
            assert_eq!(
                GS.resolve_os_momenta(&point),
                point,
                "bare metadata stays available for denominator certificates"
            );
            for (component, value) in velocity.iter().enumerate() {
                let indexed = GS.indexed_momentum(&point, &[GS.cind(component)]);
                let resolved = GS.resolve_os_momenta(&indexed);
                let expected = Atom::num(sign) * &mass_atom * value;
                assert_eq!(resolved, expected);
                assert_eq!(
                    GS.resolve_os_momenta(
                        &GS.indexed_momentum(&scaled_point, &[GS.cind(component)])
                    ),
                    &scale * &expected,
                );
                assert_eq!(
                    resolved.replace(mass).with(&scale * &mass_atom),
                    &scale * expected
                );
            }
            assert_ne!(
                point,
                GS.on_shell_momentum(EdgeIndex(8), &(Atom::num(sign) * &mass_atom), &velocity)
            );
        }
    }

    #[test]
    fn esurface_preserves_momentum_and_mass_substitutions() {
        let t = symbol!("esurface_substitution_test::t"; Scalar);
        let m = symbol!("esurface_substitution_test::m"; Scalar);
        let q = GS.emr_vec(EdgeIndex(3), GS.cind(1));
        let shift = GS.emr_mom(EdgeIndex(0), GS.cind(0));
        let energy = function!(GS.on_shell_energy, 3, q.pow(2) + Atom::var(m).pow(2));
        let sum = energy + &shift;
        let surface = GS.wrap_esurface(&sum);
        assert!(surface.contains_symbol(GS.esurface));
        for scale in [Atom::var(t), Atom::var(t).pow(-1)] {
            let replacements = [
                Replacement::new(q.clone(), &q * &scale),
                Replacement::new(Atom::var(m), Atom::var(m) * &scale),
                Replacement::new(shift.clone(), &shift * &scale),
            ];
            let transformed = surface.replace_multiple(&replacements);
            assert!(transformed.contains_symbol(t));
            // Laurent normalization can factor a power of t out of the sum.
            // Check scalar algebraic equality, including the chain rule.
            assert!(
                (GS.lower_esurfaces(&transformed) - sum.replace_multiple(&replacements))
                    .expand()
                    .is_zero()
            );
            assert!(
                (GS.lower_esurfaces(&transformed.derivative(t))
                    - sum.replace_multiple(&replacements).derivative(t))
                .expand()
                .is_zero()
            );
        }
    }

    #[test]
    fn esurface_retention_is_idempotent_and_preserves_energy_shifts() {
        let first = GS.ose(EdgeIndex(1));
        let second = GS.ose(EdgeIndex(2));
        let shift = GS.emr_mom(EdgeIndex(0), GS.cind(0));
        let cff = (&first + &second - &shift).pow(-2) / &first;
        let retained = GS.retain_esurfaces(&cff);
        assert!(retained.contains_symbol(GS.esurface));
        assert_eq!(GS.retain_esurfaces(&retained), retained);
        assert_eq!(GS.lower_esurfaces(&retained), cff);
        assert_eq!(
            GS.lower_esurfaces(&retained.replace(shift.to_pattern()).with(0)),
            cff.replace(shift.to_pattern()).with(0),
        );
    }

    #[test]
    fn esurface_series_preserves_cancelling_and_singular_leading_terms() {
        use symbolica::poly::series::SeriesDepth;

        let t = symbol!("esurface_laurent_test::t"; Scalar);
        let x = Atom::var(t);
        for sum in [
            (Atom::one() + &x).pow((1, 2)) - 1,
            (Atom::one() + &x).pow((1, 2)) - 1 - &x / 2,
            &x + x.pow(-1),
            (Atom::one() + &x).pow((1, 2)) / &x - x.pow(-1) - Atom::num((1, 2)),
            (Atom::one() + &x).pow((1, 2)) / x.pow(2) - x.pow(-2) - Atom::num((1, 2)) / &x
                + Atom::num((1, 8)),
            (Atom::one() + &x).pow((1, 2)) / &x - x.pow(-1),
        ] {
            let wrapped = GS.wrap_esurface(&sum);
            assert!((GS.lower_esurfaces(&wrapped) - &sum).expand().is_zero());
            for power in [-2, -1, 1, 2] {
                let actual = wrapped
                    .pow(power)
                    .series(t, 0, SeriesDepth::absolute(2))
                    .unwrap();
                let expected = sum
                    .pow(power)
                    .series(t, 0, SeriesDepth::absolute(2))
                    .unwrap();
                assert_eq!(actual.absolute_order(), expected.absolute_order());
                // Denominator-only identity: no graph numerator is distributed.
                assert!(
                    (GS.lower_esurfaces(&actual.to_atom()) - expected.to_atom())
                        .expand()
                        .is_zero()
                );
            }
        }
    }

    #[test]
    fn esurface_coefficients_accept_an_enclosing_momentum_deformation() {
        use symbolica::poly::series::SeriesDepth;

        let t = symbol!("esurface_nested_test::t"; Scalar);
        let u = symbol!("esurface_nested_test::u"; Scalar);
        let q = symbol!("esurface_nested_test::q"; Scalar);
        let invariant = Atom::var(q).pow(2) + 1;
        let sum = (&invariant + Atom::var(t)).pow((1, 2)) + &invariant;
        let first = GS
            .wrap_esurface(&sum)
            .pow(-1)
            .series(t, 0, SeriesDepth::absolute(2))
            .unwrap()
            .to_atom();
        let expected_first = sum
            .pow(-1)
            .series(t, 0, SeriesDepth::absolute(2))
            .unwrap()
            .to_atom();
        assert!(first.contains_symbol(GS.esurface));
        let deform = |atom: Atom| atom.replace(q).with(Atom::var(u) + 2);
        let actual = deform(first)
            .series(u, 0, SeriesDepth::absolute(1))
            .unwrap()
            .to_atom();
        let expected = deform(expected_first)
            .series(u, 0, SeriesDepth::absolute(1))
            .unwrap()
            .to_atom();
        assert!((GS.lower_esurfaces(&actual) - expected).expand().is_zero());
    }

    #[test]
    fn compact_energy_derivatives_keep_shift_and_mass_dependence() {
        crate::initialisation::test_initialise().unwrap();
        let (t, u, m, k, p, q) = symbol!(
            "compact_energy_test::t",
            "compact_energy_test::u",
            "compact_energy_test::m",
            "compact_energy_test::k",
            "compact_energy_test::p",
            "compact_energy_test::q"
        );
        let momentum = Atom::var(k) + Atom::var(t) * p + Atom::var(u) * q;
        let transverse = Atom::var(t) * q;
        let longitudinal = Atom::var(u) * p;
        let energy = function!(
            GS.on_shell_energy,
            3,
            Atom::var(m).pow(2) + momentum.pow(2) + transverse.pow(2) + longitudinal.pow(2)
        );
        assert_eq!(
            energy,
            function!(
                GS.on_shell_energy,
                3,
                Atom::var(m).pow(2)
                    + (-&momentum).expand_num().pow(2)
                    + (-&transverse).expand_num().pow(2)
                    + (-&longitudinal).expand_num().pow(2)
            )
        );
        assert_eq!(energy.derivative(m), Atom::var(m) / &energy);
        assert!(
            (energy.derivative(t) - (&momentum * p + &transverse * q) / &energy)
                .expand()
                .is_zero()
        );
        assert!(
            (energy.derivative(u) - (momentum * q + longitudinal * p) / &energy)
                .expand()
                .is_zero()
        );
        assert!(energy.contains_symbol(t) && energy.contains_symbol(u));

        let lower = |atom: Atom| {
            atom.replace(function!(GS.on_shell_energy, W_.a_, W_.prop_))
                .with(Atom::var(W_.prop_).pow((1, 2)))
        };
        let specialize = |atom: Atom| {
            atom.replace(m)
                .with(3)
                .replace(k)
                .with(4)
                .replace(p)
                .with(2)
                .replace(q)
                .with(1)
        };
        // Positive mass and a nonzero routed spatial momentum fix the same
        // energy branch in both independent constructions.
        for exponent in [-3, -1, 1, 2, 3] {
            let compact = specialize(energy.pow(exponent));
            let explicit = lower(compact.clone());
            let actual = compact
                .series(t, Atom::Zero, 3)
                .unwrap()
                .to_atom()
                .series(u, Atom::Zero, 2)
                .unwrap()
                .to_atom();
            let expected = explicit
                .series(t, Atom::Zero, 3)
                .unwrap()
                .to_atom()
                .series(u, Atom::Zero, 2)
                .unwrap()
                .to_atom();
            // Only this small scalar energy oracle is distributed; graph
            // numerators never enter this comparison.
            assert!(
                (lower(actual) - expected).expand().is_zero(),
                "power={exponent}"
            );
        }
    }

    #[test]
    fn compact_energy_invariant_matches_high_order_laurent_coefficients() {
        crate::initialisation::test_initialise().unwrap();
        let t = symbol!("compact_energy_laurent_test::t");
        let parameter = Atom::var(t);
        // Three independently shifted components feed one differentiated
        // invariant. This reproduces the order-six probe needed by a t^-6 body.
        let squared = (&parameter + 1).pow(2)
            + (&parameter * 2 + 2).pow(2)
            + (Atom::num(2) - &parameter).pow(2);
        let energy = function!(GS.on_shell_energy, 3, &squared);
        let actual = (energy / parameter.pow(6))
            .series(t, Atom::Zero, 0)
            .unwrap();
        let expected = (squared.pow((1, 2)) / parameter.pow(6))
            .series(t, Atom::Zero, 0)
            .unwrap();
        assert_eq!(actual.terms().count(), 7);
        // Native square-root series can include zero half-integer slots.
        // Compare the exact polynomials rather than zipping their storage grids.
        assert_eq!(
            actual
                .to_atom()
                .replace(function!(GS.on_shell_energy, 3, 9))
                .with(3),
            expected.to_atom()
        );
    }

    #[test]
    fn compact_energy_normalization_preserves_owners_after_cache_eviction() {
        crate::initialisation::test_initialise().unwrap();
        let t = symbol!("compact_energy_cache_test::t");
        let parameter = Atom::var(t);
        let invariant = ((&parameter + 1).pow(2) + 3).pow(2) + 4;
        let cold = function!(GS.on_shell_energy, 3, &invariant);
        let other_owner = function!(GS.on_shell_energy, 7, &invariant);
        let AtomView::Fun(cold_call) = cold.as_view() else {
            panic!("Expected an owned energy function");
        };
        let AtomView::Fun(other_call) = other_owner.as_view() else {
            panic!("Expected an owned energy function");
        };
        assert_eq!(other_call.get(0), Atom::num(7).as_view());
        assert_eq!(cold_call.get(1), other_call.get(1));
        assert_ne!(cold, other_owner);

        // More than the entry budget evicts earlier invariants. Eviction must
        // only repeat normalization, including for nested squared expressions.
        for offset in 0..1050 {
            let _ = function!(GS.on_shell_energy, 3, &invariant + offset);
        }
        assert_eq!(function!(GS.on_shell_energy, 3, &invariant), cold);
        assert_eq!(function!(GS.on_shell_energy, 7, &invariant), other_owner);

        let lower = |atom: Atom| {
            atom.replace(function!(GS.on_shell_energy, W_.a_, W_.prop_))
                .with(Atom::var(W_.prop_).pow((1, 2)))
        };
        let expected = invariant
            .pow((1, 2))
            .series(t, Atom::Zero, 3)
            .unwrap()
            .to_atom();
        for energy in [cold, other_owner] {
            let actual = energy.series(t, Atom::Zero, 3).unwrap().to_atom();
            // The coefficients contain radicals: expansion alone need not
            // identify equivalent roots such as sqrt(1/5) and sqrt(5)/5.
            let (_, difference) = (lower(actual) - &expected)
                .to_polynomial_in_algebraic_extension::<u16>(t, &[])
                .unwrap();
            assert!(difference.is_zero());
        }
    }

    #[test]
    fn orientation_collection_preserves_selector_free_factorization() {
        let (a, b, c, d, x) = symbol!(
            "selector_free_a",
            "selector_free_b",
            "selector_free_c",
            "selector_free_d",
            "selector_free_x"
        );
        let expression = (Atom::var(a) + b).pow(7) * (Atom::var(c) + d).pow(5)
            / (Atom::var(x) + Atom::one()).pow(2);
        assert_eq!(GS.collect_orientation_if(expression.clone()), expression);
        assert_eq!(GS.collect_orientation_if(expression.as_view()), expression);
    }

    #[test]
    fn orientation_collection_merges_same_key_scalar_contributions() {
        let (key, a, b, c, d, x) = symbol!(
            "selector_merge_key",
            "selector_merge_a",
            "selector_merge_b",
            "selector_merge_c",
            "selector_merge_d",
            "selector_merge_x"
        );
        let condition = Atom::var(key) - 4;
        let bodies = [
            (Atom::var(a) + b).pow(3) * (Atom::var(c) + d).pow(2),
            Atom::var(x).pow(-1),
            Atom::var(c) + 7,
        ];
        // Unrelated siblings and other branch forms keep their own conditions.
        let unrelated = Atom::var(d)
            + Symbol::IF.call_args([Atom::var(key) - 9, Atom::Zero, Atom::var(x).pow(-2)])
            + Symbol::IF.call_args([condition.clone(), Atom::var(a), Atom::Zero])
            + Symbol::IF.call_args([condition.clone(), Atom::var(b), Atom::var(c)]);
        let expression = Atom::add_many(
            bodies
                .iter()
                .map(|body| Symbol::IF.call_args([condition.clone(), Atom::Zero, body.clone()])),
        ) + &unrelated;
        let expected =
            Symbol::IF.call_args([condition, Atom::Zero, Atom::add_many(bodies)]) + unrelated;
        let collected = GS.collect_orientation_if(expression);
        assert_eq!(collected, expected);
        assert_eq!(GS.collect_orientation_if(collected.clone()), collected);
    }

    fn pairwise_orientation_branches(expression: &Atom) -> Atom {
        expression
            .replace(
                Symbol::IF.call_args([Atom::var(W_.a_), Atom::Zero, Atom::var(W_.b_)])
                    + Symbol::IF.call_args([Atom::var(W_.a_), Atom::Zero, Atom::var(W_.c_)]),
            )
            .repeat()
            .with(Symbol::IF.call_args([Atom::var(W_.a_), Atom::Zero, Atom::var(W_.b_) + W_.c_]))
    }

    #[test]
    fn orientation_collection_matches_pairwise_rule_for_nested_branches() {
        let (key, nested, a, b, c, x, opaque) = symbol!(
            "selector_nested_key",
            "selector_nested_inner",
            "selector_nested_a",
            "selector_nested_b",
            "selector_nested_c",
            "selector_nested_x",
            "selector_nested_opaque"
        );
        let branch =
            |condition, body: Atom| Symbol::IF.call_args([Atom::var(condition), Atom::Zero, body]);
        let inner = branch(nested, Atom::var(a)) + branch(nested, Atom::var(b));
        let outer =
            branch(key, branch(nested, Atom::var(a))) + branch(key, branch(nested, Atom::var(b)));
        let spectators = Atom::var(x) * branch(key, Atom::var(a))
            + 2 * branch(key, Atom::var(b))
            + branch(key, Atom::var(c)).pow(2)
            + Symbol::IF.call_args([Atom::var(key), Atom::var(a), Atom::Zero])
            + Symbol::IF.call_args([Atom::var(key), Atom::var(b), Atom::var(c)]);
        for expression in [
            outer.clone(),
            function!(opaque, &outer, &inner),
            outer.pow(3),
        ] {
            let collected = GS.collect_orientation_if(&expression);
            assert_eq!(collected, pairwise_orientation_branches(&expression));
            assert_eq!(GS.collect_orientation_if(&collected), collected);
        }
        let spectator = (Atom::var(x) + c).pow(7);
        let powered_branch = branch(key, Atom::var(c)).pow(2);
        for expression in [&spectator * &outer, &inner + &spectators, spectators] {
            let collected = GS.collect_orientation_if(&expression);
            let previous = pairwise_orientation_branches(&expression);
            // Weight absorption intentionally changes the old normal form.
            // Exhaust every selector truth value, retaining symbolic payloads;
            // only these small synthetic scalar oracle copies are distributed.
            for key_value in [0, 1] {
                for nested_value in [0, 1] {
                    let select = |atom: &Atom| {
                        atom.replace(key)
                            .with(key_value)
                            .replace(nested)
                            .with(nested_value)
                            .expand()
                    };
                    assert_eq!(select(&collected), select(&previous));
                }
            }
            for retained in [&spectator, &powered_branch] {
                if expression.contains(retained) {
                    assert!(collected.contains(retained));
                }
            }
            assert_eq!(GS.collect_orientation_if(&collected), collected);
        }
    }

    #[test]
    fn orientation_collection_revisits_callback_created_branches() {
        let (key, nested, a, b, c, d, wrapper) = symbol!(
            "selector_callback_key",
            "selector_callback_nested",
            "selector_callback_a",
            "selector_callback_b",
            "selector_callback_c",
            "selector_callback_d",
            "selector_callback_wrapper"
        );
        let branch =
            |condition, body: Atom| Symbol::IF.call_args([Atom::var(condition), Atom::Zero, body]);
        let introduced = function!(
            wrapper,
            branch(nested, Atom::var(c)) + branch(nested, Atom::var(d))
        );
        let callback = symbol!(
            "selector_callback_parent",
            norm = move |node, output| {
                if let AtomView::Fun(parent) = node
                    && parent.get_nargs() == 1
                    && let AtomView::Fun(argument) = parent.get(0)
                    && argument.get_symbol() == Symbol::IF
                {
                    **output = introduced.clone();
                }
            }
        );
        // The callback leaves the original sum alone, then creates a new sum
        // below a fresh wrapper when its argument becomes one collected IF.
        let expression = function!(
            callback,
            branch(key, Atom::var(a)) + branch(key, Atom::var(b))
        );
        let expected = function!(wrapper, branch(nested, Atom::var(c) + d));
        let collected = GS.collect_orientation_if(&expression);
        assert_eq!(collected, expected);
        assert_eq!(collected, pairwise_orientation_branches(&expression));
        assert_eq!(GS.collect_orientation_if(&collected), collected);
    }

    #[test]
    fn orientation_collection_preserves_repeated_and_cancelling_payloads() {
        let (key, a, b, x) = symbol!(
            "selector_cancel_key",
            "selector_cancel_a",
            "selector_cancel_b",
            "selector_cancel_x"
        );
        let branch = |body: Atom| Symbol::IF.call_args([Atom::var(key), Atom::Zero, body]);
        let first = (Atom::var(a) + b).pow(5);
        let second = Atom::var(x).pow(-1);
        for expression in [
            branch(first.clone()) + branch(-&first),
            branch(first.clone()) + branch(second.clone()) + branch(&first + &second),
            2 * branch(first.clone()) + branch(second.clone()) + branch(-second),
        ] {
            let collected = GS.collect_orientation_if(&expression);
            let previous = pairwise_orientation_branches(&expression);
            // Pairwise normalization may turn two equal calls into 2*IF and
            // stop matching them. Bulk collection can keep that 2 inside the
            // selected body instead. Expand only these small synthetic scalar
            // branches for an exact oracle, never a graph numerator.
            let pattern = Symbol::IF
                .call_args([Atom::var(key), Atom::var(W_.a_), Atom::var(W_.b_)])
                .to_pattern();
            for selected in [W_.a_, W_.b_] {
                assert_eq!(
                    collected
                        .replace(&pattern)
                        .with(Atom::var(selected))
                        .expand(),
                    previous
                        .replace(&pattern)
                        .with(Atom::var(selected))
                        .expand(),
                );
            }
            assert_eq!(GS.collect_orientation_if(&collected), collected);
        }
    }

    #[test]
    fn orientation_collection_bulk_groups_many_branches_without_expansion() {
        let (key, body, a, b, spectator) = symbol!(
            "selector_bulk_key",
            "selector_bulk_body",
            "selector_bulk_a",
            "selector_bulk_b",
            "selector_bulk_spectator"
        );
        let factor = (Atom::var(a) + b).pow(7);
        let branches = (0..3004)
            .map(|index| {
                let condition = function!(key, index % 4);
                let payload = function!(body, index) * &factor;
                Symbol::IF.call_args([condition, Atom::Zero, payload])
            })
            .collect::<Vec<_>>();
        let unrelated = function!(spectator, &factor).pow(-1);
        let expected = Atom::add_many((0..4).map(|group| {
            Symbol::IF.call_args([
                function!(key, group),
                Atom::Zero,
                Atom::add_many(
                    (group..3004)
                        .step_by(4)
                        .map(|index| function!(body, index) * &factor),
                ),
            ])
        })) + &unrelated;
        for expression in [
            Atom::add_many(branches.iter()) + &unrelated,
            Atom::add_many(branches.iter().rev()) + &unrelated,
        ] {
            let collected = GS.collect_orientation_if(expression);
            assert_eq!(collected, expected);
            assert_eq!(GS.collect_orientation_if(&collected), collected);
        }
    }

    #[test]
    fn uv_provenance_is_scalar_metadata_with_exact_momentum_roundtrip() {
        use idenso::tensor::SymbolicTensor;
        use spenso::structure::representation::{Minkowski, RepName};

        crate::initialisation::test_initialise().unwrap();
        let momentum = GS.emr_mom.call(2) - GS.emr_mom.call(3);
        let provenance =
            GS.uv_momentum_provenance_tag(7, UvMomentumProvenanceRole::TaylorFixed, &momentum);
        assert_eq!(
            GS.uv_momentum_provenance_data(provenance.as_view()),
            Some((
                EdgeIndex(7),
                UvMomentumProvenanceRole::TaylorFixed,
                momentum.clone()
            ))
        );
        let mink = Minkowski {}.new_rep(GS.dim);
        for slot in [mink.to_symbolic([]), mink.to_symbolic([Atom::num(1)])] {
            let tagged = function!(GS.emr_mom, &provenance, &slot);
            let typed = SymbolicTensor::infer(tagged.clone()).unwrap();
            assert_eq!(typed.expression(), &tagged);
            assert_eq!(
                GS.erase_uv_momentum_provenance(typed.expression()),
                GS.indexed_momentum(&momentum, &[slot])
            );
        }
    }

    #[test]
    fn orientation_collection_merges_nested_same_key_contributions() {
        let (key, a, b, c, d, w) = symbol!(
            "nested_selector_key",
            "nested_selector_a",
            "nested_selector_b",
            "nested_selector_c",
            "nested_selector_d",
            "nested_selector_w"
        );
        let select = |offset: i64, body: Atom| {
            Symbol::IF.call_args([Atom::var(key) - offset, Atom::Zero, body])
        };
        // Same-key siblings inside a factorized sum merge where they are, and
        // a merged body merges the same-key siblings it now contains.
        let inner = select(1, select(2, Atom::var(a)))
            + select(1, select(2, Atom::var(b)))
            + select(3, Atom::var(c));
        let expression = Atom::var(w) * (Atom::var(d) + Atom::var(w) * inner)
            + select(1, Atom::var(c))
            + select(1, Atom::var(d));
        let expected = Atom::var(w)
            * (Atom::var(d)
                + Atom::var(w)
                    * (select(1, select(2, Atom::var(a) + b)) + select(3, Atom::var(c))))
            + select(1, Atom::var(c) + d);
        let collected = GS.collect_orientation_if(expression);
        assert_eq!(collected, expected);
        assert_eq!(GS.collect_orientation_if(collected.clone()), collected);
    }

    #[test]
    fn orientation_collection_merges_weighted_factorized_contributions() {
        let (key, a, b, c, d, x, alias) = symbol!(
            "weighted_selector_key",
            "weighted_selector_a",
            "weighted_selector_b",
            "weighted_selector_c",
            "weighted_selector_d",
            "weighted_selector_x",
            "weighted_selector_alias"
        );
        let condition = Atom::var(key) - 4;
        let select = |body: Atom| Symbol::IF.call_args([condition.clone(), Atom::Zero, body]);
        let numerator = (Atom::var(a) + b).pow(3);
        let weight = (Atom::var(c) + d).pow(2);
        let handle = alias.call(7);
        let expression = select(numerator.clone()) * &weight
            + Atom::i() * &handle * select(Atom::var(x).pow(-1));
        let expected = select(numerator * weight + Atom::i() * &handle * Atom::var(x).pow(-1));
        let collected = GS.collect_orientation_if(expression);
        assert_eq!(collected, expected);
        assert!(collected.contains(&handle));
        assert_eq!(GS.collect_orientation_if(collected.clone()), collected);
    }

    #[test]
    fn orientation_collection_reaches_nested_weighted_guards() {
        let (key, a, b, c, x, y, z, w) = symbol!(
            "nested_weighted_key",
            "nested_weighted_a",
            "nested_weighted_b",
            "nested_weighted_c",
            "nested_weighted_x",
            "nested_weighted_y",
            "nested_weighted_z",
            "nested_weighted_w"
        );
        let select = |body: Atom| Symbol::IF.call_args([Atom::var(key), Atom::Zero, body]);
        let expression = Atom::var(w) * (select(Atom::var(a)) * x + select(Atom::var(b)) * y)
            + select(Atom::var(c)) * z;
        let expected =
            select(Atom::var(w) * (Atom::var(a) * x + Atom::var(b) * y) + Atom::var(c) * z);
        let collected = GS.collect_orientation_if(expression);
        assert_eq!(collected, expected);
        assert_eq!(GS.collect_orientation_if(collected.clone()), collected);
    }

    #[test]
    fn orientation_collection_preserves_distinct_weighted_conditions() {
        let (key, a, b, c, x, y, z) = symbol!(
            "distinct_weighted_key",
            "distinct_weighted_a",
            "distinct_weighted_b",
            "distinct_weighted_c",
            "distinct_weighted_x",
            "distinct_weighted_y",
            "distinct_weighted_z"
        );
        let select = |offset: i64, body: Atom| {
            Symbol::IF.call_args([Atom::var(key) - offset, Atom::Zero, body])
        };
        let ordinary = Symbol::IF.call_args([Atom::var(key), Atom::var(b), Atom::var(c)]);
        let expression = select(1, Atom::var(a)) * x + select(2, Atom::var(b)) * y + &ordinary * z;
        let expected = select(1, Atom::var(a) * x) + select(2, Atom::var(b) * y) + ordinary * z;
        let collected = GS.collect_orientation_if(expression);
        assert_eq!(collected, expected);
        assert_eq!(GS.collect_orientation_if(collected.clone()), collected);
    }

    #[test]
    fn orientation_collection_preserves_literal_function_tags() {
        let (key, a, w, alias) = symbol!(
            "opaque_weighted_key",
            "opaque_weighted_a",
            "opaque_weighted_w",
            "opaque_weighted_alias"
        );
        let tag = Symbol::IF.call_args([Atom::var(key), Atom::Zero, Atom::var(a)]) * Atom::var(w);
        let handle = alias.call(tag);
        assert_eq!(GS.collect_orientation_if(handle.clone()), handle);
    }

    #[test]
    fn orientation_collection_combines_same_key_guard_factors() {
        let (key, a, b, w) = symbol!(
            "guard_product_key",
            "guard_product_a",
            "guard_product_b",
            "guard_product_w"
        );
        let select = |body: Atom| Symbol::IF.call_args([Atom::var(key), Atom::Zero, body]);
        let expression = select(Atom::var(a)) * select(Atom::var(b)) * w;
        let expected = select(Atom::var(a) * b * w);
        let collected = GS.collect_orientation_if(expression);
        assert_eq!(collected, expected);
        assert_eq!(GS.collect_orientation_if(collected.clone()), collected);
    }

    #[test]
    fn canonical_uv_class_references_are_distinct_from_physical_owners() {
        let class = UvDenominatorClassId(7);
        assert_eq!(
            GS.uv_class_data(GS.uv_class_ref(class).as_view()),
            Some(class)
        );
        assert_eq!(GS.uv_class_data(Atom::num(7).as_view()), None);
        assert_eq!(GS.uv_class_data(GS.uv_class.call(-1).as_view()), None);
        assert_eq!(
            GS.uv_class_data(
                GS.uv_class
                    .call_args([Atom::num(7), Atom::num(8)])
                    .as_view()
            ),
            None
        );
        let provenance = GS.uv_momentum_provenance_tag(
            7,
            UvMomentumProvenanceRole::TaylorFixed,
            GS.emr_mom.call(7),
        );
        assert_eq!(GS.uv_class_data(provenance.as_view()), None);
    }

    #[test]
    fn test_print() {
        let p = GS.emr_mom(EdgeIndex(1), Atom::Zero);

        assert_snapshot!(p.log_print(None),@"[35mq₁[0m(0)")
    }

    #[test]
    fn typst_uses_structural_accents_and_scripts() {
        let expressions = [
            (
                GS.ubar.call_args([Atom::num(1), Atom::num(11)]),
                "overline(u)_1^(11)",
            ),
            (
                GS.vbar.call_args([Atom::num(2), Atom::num(12)]),
                "overline(v)_2^(12)",
            ),
            (GS.ose.call_args([Atom::num(3)]), r#"E^("os")_(3)"#),
            (
                GS.emr_vec(EdgeIndex(4), Minkowski {}.new_rep(4).to_symbolic([])),
                "arrow(q)_4",
            ),
            (
                GS.epsilonbar.call_args([Atom::num(5), Atom::num(15)]),
                "overline(epsilon.alt)_5^(15)",
            ),
            (
                GS.emr_mom.call_args([Atom::num(6), Atom::num(16)]),
                "q_6^(16)",
            ),
        ];

        for (expression, expected) in expressions {
            assert_eq!(
                expression
                    .printer(SpensoPrintSettings::typst_options())
                    .to_string(),
                expected
            );
        }
    }

    #[test]
    fn typst_syntax_requires_real_typst_mode() {
        let mut symbolica = SpensoPrintSettings::typst().nice_symbolica();
        symbolica.color_builtin_symbols = false;

        assert_eq!(
            GS.ubar
                .call_args([Atom::num(1), Atom::num(11)])
                .printer(symbolica.clone())
                .to_string(),
            "u̅₁^(11)"
        );
        assert_eq!(
            Atom::var(GS.m_uv_vacuum).printer(symbolica).to_string(),
            "mUV"
        );

        assert_eq!(
            Atom::var(GS.m_uv_expansion)
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            r#"m_"UVexp""#
        );
        assert_eq!(
            Atom::var(GS.m_uv_vacuum)
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            r#"m_"UV""#
        );
        assert_eq!(
            Atom::var(GS.mu_r_sq)
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            "mu_r^2"
        );
    }
}
