use std::sync::LazyLock;

use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, Orientation};

use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    shadowing::symbolica_utils::SpensoPrintSettings,
    structure::{
        abstract_index::AIND_SYMBOLS,
        concrete_index::ExpandedIndex,
        representation::{Minkowski, RepName, Representation},
        slot::{DummyAind, IsAbstractSlot},
    },
    utils::{to_subscript, to_superscript},
};
use symbolica::prelude::*;
use symbolica::printer::{PrintState, PrintUserData};

use crate::{cff::orientations::GraphOrientation, graph::LoopMomentumBasis, numerator::aind::Aind};

use super::symbolica_ext::CallSymbol;

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
    pub override_if: Symbol,
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
    pub deta_lu_cut: Symbol,
    pub deta_left_th: Symbol,
    pub deta_right_th: Symbol,
    pub vbar: Symbol,
    pub ose: Symbol,
    pub energy: Symbol,
    pub v: Symbol,
    pub u: Symbol,
    pub color_wrap: Symbol,
    /// Epsilon for dimensional regularization
    pub dim_epsilon: Symbol,

    pub epsilon: Symbol,
    pub epsilonbar: Symbol,
    pub rescale: Symbol,
    /// Bookkeeping scale for integrated vacuum masses and consumed loop measures.
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
    pub sign: Symbol,
    pub theta: Symbol,
    pub broadcasting_sqrt: Symbol,
    ///for selecting orientations at generation
    pub selected: Symbol,

    /// For analytic UV profiling
    pub expansion: Symbol,
    ///For selecting a concete index.
    pub delta_vec: Symbol,
    ///Q(<edgeid>,index___)
    pub emr_mom: Symbol,
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
    /// Structured forest-DOT approximation record.
    pub t_op: Symbol,

    pub nc2_1: Symbol,
    pub num: Symbol,
    ///denominator wrapper, den(<edge_id>,<momentum>,<mass>,<full_expr>) (no power and should not be multiplied in but divided!)
    pub den: Symbol,
    pub radius_left: Symbol,
    pub radius_star_left: Symbol,
    pub uv_damp_plus_left: Symbol,
    pub uv_damp_minus_left: Symbol,
    pub radius_right: Symbol,
    pub radius_star_right: Symbol,
    pub uv_damp_plus_right: Symbol,
    pub uv_damp_minus_right: Symbol,
}

impl GammaloopSymbols {
    pub fn collect_orientation_if<'a>(
        &self,
        arg: impl Into<AtomOrView<'a>>,
        with_override: bool,
    ) -> Atom {
        arg.into()
            .replace(self.sign_theta(W_.a_))
            .with(Symbol::IF.f(Atom::var(W_.a_) + 1))
            .replace(Symbol::IF.f(W_.a_) * Symbol::IF.f(W_.b_))
            .repeat()
            .with(Symbol::IF.f(W_.a_ * W_.b_))
            .replace(Symbol::IF.f(W_.a_) * W_.b___)
            .with(Symbol::IF.f([Atom::var(W_.a_), Atom::var(W_.b___), Atom::Zero]))
            .replace(Symbol::IF.f([Atom::var(W_.a_), Atom::Zero]))
            .with(Symbol::IF.f([Atom::var(W_.a_), Atom::one(), Atom::Zero]))
            .replace(Symbol::IF.f([Atom::var(W_.a_)]))
            .with(Symbol::IF.f([Atom::var(W_.a_), Atom::one(), Atom::Zero]))
            .replace(Symbol::IF.f([W_.a_, W_.b_, W_.c_]))
            .with({
                if with_override {
                    Symbol::IF.f([
                        Atom::var(self.override_if),
                        Atom::var(W_.b_),
                        Symbol::IF.f([W_.a_, W_.b_, W_.c_]),
                    ])
                } else {
                    Symbol::IF.f([W_.a_, W_.b_, W_.c_])
                }
            })
    }

    pub fn den<'a>(
        &self,
        eid: impl Into<AtomOrView<'a>>,
        mom: impl Into<AtomOrView<'a>>,
        mass: impl Into<AtomOrView<'a>>,
        full_expr: impl Into<AtomOrView<'a>>,
    ) -> Atom {
        self.den.f(&[
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

macro_rules! spenso_print_scripted_indexed {
    ($a:ident, $opt:ident, $symbol:expr) => {{
        match $opt.custom_print_mode.get("spenso") {
            Some(PrintUserData::Integer(i)) => {
                let SpensoPrintSettings {
                    parens,
                    symbol_scripts,
                    commas,
                    with_dim,
                    ..
                } = SpensoPrintSettings::from(*i as usize);

                let AtomView::Fun(f) = $a else {
                    return None;
                };

                let mut argiter = f.iter();
                let id = argiter.next().unwrap();
                let Ok(i) = usize::try_from(id) else {
                    return None;
                };

                let mut out = $symbol.to_string();
                out.push_str(&to_subscript(i as isize));
                if $opt.color_builtin_symbols {
                    out = nu_ansi_term::Color::Magenta.paint(out).to_string();
                }

                let mut printed_args = false;
                for arg in argiter {
                    let hidden_representation = matches!(
                        arg,
                        AtomView::Fun(a)
                            if a.get_symbol().has_tag(&SPENSO_TAG.representation)
                                && a.get_nargs() == 1
                                && !with_dim
                    );
                    if hidden_representation {
                        continue;
                    }

                    if printed_args {
                        out.push(if commas { ',' } else { ' ' });
                    } else {
                        if symbol_scripts {
                            out.push('^');
                        }
                        if parens {
                            out.push('(');
                        }
                        printed_args = true;
                    }

                    arg.format(&mut out, $opt, PrintState::new()).unwrap();
                }
                if printed_args && parens {
                    out.push(')');
                }
                Some(out)
            }
            _ => None,
        }
    }};
}

macro_rules! spenso_print_simple_indexed {
    ($a:ident, $opt:ident, $symbol:expr) => {{
        match $opt.custom_print_mode.get("spenso") {
            Some(PrintUserData::Integer(_)) => {
                let AtomView::Fun(f) = $a else {
                    return None;
                };

                let mut out = $symbol.to_string();
                let mut args = f.iter();

                let id = args.next().unwrap();
                let Ok(i) = usize::try_from(id) else {
                    return None;
                };

                out.push_str(&to_subscript(i as isize));
                let mut first = true;
                for arg in args {
                    if first {
                        first = false;
                        out.push('(');
                    } else {
                        out.push(',');
                    }
                    arg.format(&mut out, $opt, PrintState::new()).unwrap();
                }
                if !first {
                    out.push(')');
                }
                Some(out)
            }
            _ => None,
        }
    }};
}

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
        print = |a, opt, _state| {
            match opt.custom_print_mode.get("spenso") {
                Some(PrintUserData::Integer(_i)) => {
                    let AtomView::Fun(f) = a else {
                        return None;
                    };

                    let mut out = "ᵘ".to_string();
                    let mut first = true;
                    for arg in f.iter() {
                        let Ok(i) = isize::try_from(arg) else {
                            return None;
                        };

                        if !first {
                            out.push('.');
                        } else {
                            first = false;
                        }
                        out.push_str(&to_superscript(i));
                    }
                    Some(out)
                }
                _ => None,
            }
        },
        tags = [SPENSO_TAG.index.clone()]
    ),
    edgeaind: symbol!(
        "edge",
        print = |a, opt, _state| {
            match opt.custom_print_mode.get("spenso") {
                Some(PrintUserData::Integer(_i)) => {
                    let AtomView::Fun(f) = a else {
                        return None;
                    };

                    let mut out = "ᵉ".to_string();
                    let mut first = true;
                    for arg in f.iter() {
                        let Ok(i) = isize::try_from(arg) else {
                            return None;
                        };

                        if !first {
                            out.push('.');
                        }
                        first = false;

                        out.push_str(&to_superscript(i));
                    }
                    Some(out)
                }
                _ => None,
            }
        },
        tags = [SPENSO_TAG.index.clone()]
    ),
    vertexaind: symbol!(
        "vertex",
        print = |a, opt, _state| {
            match opt.custom_print_mode.get("spenso") {
                Some(PrintUserData::Integer(_i)) => {
                    let AtomView::Fun(f) = a else {
                        return None;
                    };

                    let mut out = "ᵛ".to_string();

                    let mut first = true;
                    for arg in f.iter() {
                        let Ok(i) = isize::try_from(arg) else {
                            return None;
                        };

                        if !first {
                            out.push('.');
                        }
                        first = false;

                        out.push_str(&to_superscript(i));
                    }
                    Some(out)
                }
                _ => None,
            }
        },
        tags = [SPENSO_TAG.index.clone()]
    ),
    dummyaind: symbol!(
        "dummy",
        print = |a, opt, _state| {
            match opt.custom_print_mode.get("spenso") {
                Some(PrintUserData::Integer(_i)) => {
                    let AtomView::Fun(f) = a else {
                        return None;
                    };

                    let mut out = "ᵈ".to_string();
                    let mut first = true;
                    for arg in f.iter() {
                        let Ok(i) = isize::try_from(arg) else {
                            return None;
                        };

                        if !first {
                            out.push('.');
                        }
                        first = false;

                        out.push_str(&to_superscript(i));
                    }
                    Some(out)
                }
                _ => None,
            }
        },
        tags = [SPENSO_TAG.index.clone()]
    ),
    hedgeaind: symbol!(
        "hedge",
        print = |a, opt, _state| {
            match opt.custom_print_mode.get("spenso") {
                Some(PrintUserData::Integer(i)) => {
                    let AtomView::Fun(f) = a else {
                        return None;
                    };
                    let SpensoPrintSettings {
                        index_subscripts, ..
                    } = SpensoPrintSettings::from(*i as usize);

                    let mut out = "".to_string();
                    let mut first = true;
                    for arg in f.iter() {
                        let Ok(i) = isize::try_from(arg) else {
                            return None;
                        };

                        if !first {
                            out.push('.');
                        }
                        first = false;
                        if index_subscripts {
                            out.push_str(&to_superscript(i));
                        } else {
                            out.push_str(&to_subscript(i));
                        }
                    }
                    Some(out)
                }
                _ => None,
            }
        },
        tags = [SPENSO_TAG.index.clone()]
    ),
    override_if: symbol!("override_if"),
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
        "gammalooprs::uv::Approx",
        print = |a, opt, _state| spenso_print_uv_unary!(a, opt, "K[", "]")
    ),
    uv_integrate: symbol!(
        "gammalooprs::uv::Integrate",
        print = |a, opt, _state| spenso_print_uv_unary!(a, opt, "⟨", "⟩")
    ),
    uv_series: symbol!(
        "gammalooprs::uv::Series",
        print = |a, opt, _state| spenso_print_uv_unary!(a, opt, "Σ(", ")")
    ),
    uv_truncate: symbol!(
        "gammalooprs::uv::Truncate",
        print = |a, opt, _state| spenso_print_uv_unary!(a, opt, "Tr(", ")")
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
    t_op: symbol!(
        "T",
        print = |a, opt, _state| {
            if !opt.mode.is_typst() {
                return None;
            }

            let mut out = "#T".to_string();
            if let AtomView::Fun(f) = a {
                out.push('(');
                let mut first = true;
                for arg in f.iter() {
                    if first {
                        first = false;
                    } else {
                        out.push(',');
                    }
                    arg.format(&mut out, opt, PrintState::new()).ok()?;
                }
                out.push(')');
            }

            Some(out)
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
                    .f(args.as_slice())
                    .replace(GS._linear.f(&[W_.a___]))
                    .with(arg.get_symbol().f(&[W_.a___]));
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
    expansion: symbol!("expansion"),
    rescale_star: symbol!("t⃰"),
    hfunction_lu_cut: symbol!("h_lu_cut"),
    hfunction_left_th: symbol!("h_left_th"),
    hfunction_right_th: symbol!("h_right_th"),
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
                Some(PrintUserData::Integer(i)) => {
                    let SpensoPrintSettings { .. } = SpensoPrintSettings::from(*i as usize);
                    if SpensoPrintSettings::from(*i as usize).is_typst() {
                        Some("m_\"UVexp\"".to_string())
                    } else {
                        None
                    }
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
                Some(PrintUserData::Integer(i)) => {
                    let SpensoPrintSettings { .. } = SpensoPrintSettings::from(*i as usize);
                    if SpensoPrintSettings::from(*i as usize).is_typst() {
                        Some("m_\"UV\"".to_string())
                    } else {
                        None
                    }
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
                Some(PrintUserData::Integer(i)) => {
                    let SpensoPrintSettings { .. } = SpensoPrintSettings::from(*i as usize);
                    if SpensoPrintSettings::from(*i as usize).is_typst() {
                        Some("mu_r^2".to_string())
                    } else {
                        None
                    }
                }
                _ => None,
            }
        }
    ),
    delta_vec: ETS.delta,
    num: symbol!("num"),
    den: symbol!(
        "denom",
        der = |_, arg, out| {
            if arg != 3 {
                **out = Atom::Zero;
            } else {
                **out = Atom::num(1);
            }
        }
    ),
    ubar: symbol!(
        "ubar",
        print = |a, opt, _state| { spenso_print_scripted_indexed!(a, opt, "u̅") },
        tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
    ),
    vbar: symbol!(
        "vbar",
        print = |a, opt, _state| { spenso_print_scripted_indexed!(a, opt, "v̅") },
        tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
    ),
    dot: symbol!("dot"),
    dim: symbol!("dim"),
    v: symbol!(
        "v",
        print = |a, opt, _state| { spenso_print_scripted_indexed!(a, opt, "v") },
        tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
    ),
    u: symbol!(
        "u",
        print = |a, opt, _state| { spenso_print_scripted_indexed!(a, opt, "u") },
        tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
    ),
    emr_mom: symbol!(
        "Q",
        print = |a, opt, _state| { spenso_print_scripted_indexed!(a, opt, "q") },
        tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
    ),
    orientation_delta: symbol!("orientation_delta"),
    emr_vec: symbol!(
        "Q3",
        norm = |f, out| {
            if let AtomView::Fun(ff) = f
                && ff.get_nargs() == 2
            {
                let mut iter = ff.iter();
                let eid = iter.next().unwrap();
                if let AtomView::Fun(cind) = iter.next().unwrap()
                    && cind.get_symbol() == AIND_SYMBOLS.cind
                    && let Some(i) = cind.iter().next()
                    && let Ok(i) = i64::try_from(i)
                {
                    if i == 0 {
                        **out = Atom::Zero;
                    } else {
                        **out = get_symbol!("Q").unwrap().f(&[eid, cind.as_view()])
                    }
                }
            }
        },
        print = |a, opt, _state| { spenso_print_scripted_indexed!(a, opt, "q⃗") },
        tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
    ),
    ose: symbol!(
        "OSE"; Scalar;
        print = |a, opt, _state| { spenso_print_simple_indexed!(a, opt, "Eᵒˢ") },
            der = |_, arg, out| {
                if arg == 1 {
                    **out = Atom::num(1);
                }
            }
    ),
    energy: symbol!(
        "E",
        print = |a, opt, _state| { spenso_print_simple_indexed!(a, opt, "E") }
    ),
    external_mom: symbol!(
        "P",
        print = |a, opt, _state| { spenso_print_scripted_indexed!(a, opt, "p") },
        tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
    ),
    loop_mom: symbol!(
        "K",
        print = |a, opt, _state| { spenso_print_scripted_indexed!(a, opt, "k") },
        tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
    ),
    epsilon: symbol!(
        "ϵ",
        print = |a, opt, _state| { spenso_print_scripted_indexed!(a, opt, "ϵ") },
        tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
    ),
    pi: Symbol::PI,
    color_wrap: symbol!("color"),
    epsilonbar: symbol!(
        "ϵbar",
        print = |a, opt, _state| { spenso_print_scripted_indexed!(a, opt, "ϵ̅") },
        tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
    ),
    coeff: symbol!("coef"),
    radius_left: symbol!("r_left"),
    radius_star_left: symbol!("r⃰_left"),
    uv_damp_plus_left: symbol!("damp_plus_left"),
    uv_damp_minus_left: symbol!("damp_minus_left"),
    radius_right: symbol!("r_right"),
    radius_star_right: symbol!("r⃰_right"),
    uv_damp_plus_right: symbol!("damp_plus_right"),
    uv_damp_minus_right: symbol!("damp_minus_right"),
};
}

impl GammaloopSymbols {
    pub fn integrand<O: GraphOrientation>(&self, i: usize, orientation: &O) -> Atom {
        let args = orientation
            .orientation()
            .iter()
            .map(|(_, ori)| match ori {
                Orientation::Default => Atom::num(1),
                Orientation::Reversed => Atom::num(-1),
                Orientation::Undirected => Atom::num(0),
            })
            .collect_vec();

        FunctionBuilder::new(self.integrand)
            .add_arg(i)
            .add_args(&args)
            .finish()
    }

    pub fn wrap_tree_denoms<'a>(&self, arg: impl Into<AtomOrView<'a>>) -> Atom {
        self.tree_denom_wrapper.f(&[arg.into().as_view()])
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

        // let dot = self.emr_vec_index(e, mink.to_atom()) * self.emr_vec_index(e, mink.to_atom())
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
        //                 .f(&[a])
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
                        .f(&[a])
                        .pow(Atom::num(exp.numerator()));
                }
            }
        })
    }

    pub(crate) fn emr_mom<'a>(&self, e: EdgeIndex, arg: impl Into<AtomOrView<'a>>) -> Atom {
        let a = arg.into();
        function!(self.emr_mom, usize::from(e) as i64, a.as_view())
    }

    pub(crate) fn localizing_integrand(&self, lmb: &LoopMomentumBasis) -> Atom {
        // Normalize each factor with int d^3k / (|k|^2 + rls^2)^2 = pi^2 / rls,
        // so the localizing integrand itself integrates to one.
        let pi_atom = (Symbol::PI).to_atom();
        let normalization_term_integral = (pi_atom.pow(2)) / GS.renormalization_localization_scale;

        let mut res = Atom::one();

        for l in lmb.loop_edges.iter() {
            //TODO: Add orientation localisation prefactor (Sum of valid orientation thetas)/(number of valid orientations)
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

    pub(crate) fn emr_vec(&self, e: EdgeIndex) -> Atom {
        function!(GS.emr_vec, usize::from(e) as i64)
    }

    pub(crate) fn emr_vec_index<'a>(&self, e: EdgeIndex, index: impl Into<AtomOrView<'a>>) -> Atom {
        function!(GS.emr_vec, usize::from(e) as i64, index.into().as_view())
    }

    pub(crate) fn cind(&self, e: usize) -> Atom {
        AIND_SYMBOLS.cind.f([e as i64])
    }

    pub(crate) fn energy_delta<'a>(&self, index: impl Into<AtomOrView<'a>>) -> Atom {
        self.delta_vec.f(&[self.cind(0), index.into().into_owned()])
    }

    pub(crate) fn ose_full(
        &self,
        e: EdgeIndex,
        lmb_id: EdgeIndex,
        e_mass: Atom,
        index: Option<Atom>,
        inner_product: bool,
    ) -> Atom {
        let m2 = &e_mass * &e_mass;

        let mink: Representation<Minkowski> = Minkowski {}.new_rep(4); //.slot(Aind::new_dummy());
        let q3q3 = if inner_product {
            mink.inner_product(self.emr_vec(e), self.emr_vec(e))
        } else {
            let mink = mink.slot::<Aind, Aind>(Aind::new_dummy()).to_atom();

            self.emr_vec_index(e, mink.as_view()) * self.emr_vec_index(e, mink.as_view())
        };

        let ose = function!(self.ose, lmb_id.0, (m2 - q3q3)).pow((1, 2));

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
        inner_product: bool,
    ) -> Replacement {
        let id = Minkowski {}.to_symbolic([W_.a__]);
        Replacement::new(
            self.emr_mom(e, &id).to_pattern(),
            self.emr_vec_index(e, &id)
                + self.ose_full(e, lmb_id, e_mass, Some(id), inner_product) * sign_atom(e),
        )
    }

    /// Add the sign by splitting Q(i,mu)-> Q3(i,mu)+OSE(i)*σ(i)*δ(cind(0),mu)
    pub fn split_mom_pattern_simple(&self, e: EdgeIndex) -> Replacement {
        let eidc = usize::from(e) as i64;
        let index = Minkowski {}.to_symbolic([W_.a__]);
        Replacement::new(
            self.emr_mom(e, &index).to_pattern(),
            function!(GS.ose, eidc) * sign_atom(e) * self.energy_delta(&index)
                + function!(GS.emr_vec, eidc, &index),
        )
    }

    pub(crate) fn add_parametric_sign(&self, e: EdgeIndex) -> Replacement {
        Replacement::new(
            self.emr_mom(e, AIND_SYMBOLS.cind.f(&[Atom::Zero]))
                .to_pattern(),
            sign_atom(e) * self.ose(e),
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

    use crate::utils::symbolica_ext::LogPrint;

    use super::*;

    #[test]
    fn test_print() {
        let p = GS.emr_mom(EdgeIndex(1), Atom::Zero);

        assert_snapshot!(p.log_print(None),@"[35mq₁[0m(0)")
    }
}
