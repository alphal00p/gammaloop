use std::sync::LazyLock;

use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, Orientation};

use spenso::{
    network::parsing::SPENSO_TAG,
    structure::{
        abstract_index::AIND_SYMBOLS,
        concrete_index::ExpandedIndex,
        representation::{Minkowski, RepName},
        slot::{DummyAind, IsAbstractSlot, Slot},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    domains::rational::Rational,
    function,
    id::Replacement,
    symbol,
};

use crate::{cff::expression::GraphOrientation, numerator::aind::Aind};

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
    pub is_function: Symbol,
    pub is_symbol: Symbol,
    pub ufozero: Symbol,
    pub _linear: Symbol,
    pub linearize: Symbol,
    pub spensocind: Symbol,
    pub loop_mom: Symbol,
    pub edgeid: Symbol,
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
    pub epsilon: Symbol,
    pub epsilonbar: Symbol,
    pub rescale: Symbol,
    pub rescale_star: Symbol,
    pub pi: Symbol,
    pub m_uv: Symbol,
    pub m_uv_int: Symbol,
    pub mu_r_sq: Symbol,
    pub sign: Symbol,
    pub theta: Symbol,
    pub broadcasting_sqrt: Symbol,
    ///for selecting orientations at generation
    pub selected: Symbol,

    ///For selecting a concete index.
    pub delta_vec: Symbol,
    pub emr_mom: Symbol,
    pub emr_vec: Symbol,
    pub dot: Symbol,
    pub external_mom: Symbol,
    pub dim: Symbol,
    pub coeff: Symbol,

    pub if_sigma: Symbol,

    pub nc2_1: Symbol,
    pub top: Symbol,
    pub num: Symbol,
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
    pub(crate) fn sign_symbol(&self, edge: EdgeIndex) -> Symbol {
        symbol!(format!("{}{}", self.sign, edge))
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

pub static GS: LazyLock<GammaloopSymbols> = LazyLock::new(|| GammaloopSymbols {
    ufozero: symbol!(
        "UFO::ZERO",
        norm = |_, out| {
            **out = Atom::Zero;
        }
    ),
    if_sigma: symbol!("if_sigma"),
    is_function: symbol!("is_function"),
    is_symbol: symbol!("is_symbol"),
    nc2_1: symbol!("NC2_1"),
    rescale: symbol!("t"),
    _linear: symbol!("_linear";Linear),
    linearize: symbol!(
        "linearize",
        norm = |f, out| {
            if let AtomView::Fun(ff) = f
                && ff.get_nargs() == 1
                && let AtomView::Fun(arg) = ff.iter().next().unwrap()
            {
                let args = arg.iter().collect_vec();
                **out = GS
                    ._linear
                    .f(args.as_slice())
                    .replace(GS._linear.f(&[W_.a___]))
                    .with(arg.get_symbol().f(&[W_.a___]));
            }
        }
    ),
    broadcasting_sqrt: symbol!(
        "broadcasting_sqrt",
        tag = SPENSO_TAG.tag,
        der = |a, _, out| {
            **out = Atom::num(1) / (Atom::num(2) * a);
        }
    ),
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
    sign: symbol!("σ"),
    selected: symbol!("selected"),
    theta: symbol!("θ"),
    spensocind: symbol!("spenso::cind"),
    m_uv: symbol!("mUV"),
    m_uv_int: symbol!("mUVI"),
    mu_r_sq: symbol!(format!("{}::μᵣ²", vakint::NAMESPACE)),
    delta_vec: symbol!(
        "δ",
        norm = |f, out| {
            if let AtomView::Fun(ff) = f
                && ff.get_nargs() == 2
            {
                let mut iter = ff.iter();
                let matchid = iter.next().unwrap();
                if let AtomView::Fun(cind) = iter.next().unwrap()
                    && cind.get_symbol() == AIND_SYMBOLS.cind
                {
                    if cind.as_view() != matchid {
                        **out = Atom::Zero;
                    } else {
                        **out = Atom::num(1);
                    }
                }
            }
        }
    ),
    top: symbol!("Top"),
    num: symbol!("num"),
    den: symbol!("den"),
    ubar: symbol!("ubar"),
    vbar: symbol!("vbar"),
    dot: symbol!("dot"),
    dim: symbol!("dim"),
    v: symbol!("v"),
    u: symbol!("u"),
    emr_mom: symbol!("Q"),
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
                        **out = symbol!("Q").f(&[eid, cind.as_view()])
                    }
                }
            }
        }
    ),
    ose: symbol!("OSE"),
    energy: symbol!("E"),
    external_mom: symbol!("P"),
    loop_mom: symbol!("K"),
    epsilon: symbol!("ϵ"),
    pi: Symbol::PI,
    color_wrap: symbol!("color"),
    epsilonbar: symbol!("ϵbar"),
    coeff: symbol!("coef"),
    radius_left: symbol!("r_left"),
    radius_star_left: symbol!("r⃰_left"),
    uv_damp_plus_left: symbol!("damp_plus_left"),
    uv_damp_minus_left: symbol!("damp_minus_left"),
    radius_right: symbol!("r_right"),
    radius_star_right: symbol!("r⃰_right"),
    uv_damp_plus_right: symbol!("damp_plus_right"),
    uv_damp_minus_right: symbol!("damp_minus_right"),
});

impl GammaloopSymbols {
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

        let pat = function!(self.emr_vec, W_.a_, mink).npow(2);
        let rhs = function!(
            self.emr_mom,
            W_.a_,
            Atom::from(ExpandedIndex::from_iter([1]))
        )
        .npow(2)
            + function!(
                self.emr_mom,
                W_.a_,
                Atom::from(ExpandedIndex::from_iter([2]))
            )
            .npow(2)
            + function!(
                self.emr_mom,
                W_.a_,
                Atom::from(ExpandedIndex::from_iter([3]))
            )
            .npow(2);
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

    pub(crate) fn emr_vec(&self, e: EdgeIndex) -> Atom {
        function!(GS.emr_vec, usize::from(e) as i64)
    }

    pub(crate) fn emr_vec_index<'a>(&self, e: EdgeIndex, index: impl Into<AtomOrView<'a>>) -> Atom {
        function!(GS.emr_vec, usize::from(e) as i64, index.into().as_view())
    }

    pub(crate) fn cind(&self, e: usize) -> Atom {
        AIND_SYMBOLS.cind.f([e as i64])
    }
    pub(crate) fn delta_vec<'a>(&self, e: usize, index: impl Into<AtomOrView<'a>>) -> Atom {
        function!(GS.delta_vec, self.cind(e), index.into().as_view())
    }

    pub(crate) fn energy_delta<'a>(&self, index: impl Into<AtomOrView<'a>>) -> Atom {
        self.delta_vec
            .f(&[self.spensocind.f(&[0]), index.into().into_owned()])
    }

    pub(crate) fn ose_full(&self, e: EdgeIndex, e_mass: Atom, index: Option<Atom>) -> Atom {
        let eidc = usize::from(e) as i64;
        let m2 = &e_mass * &e_mass;

        let mink: Slot<Minkowski, Aind> = Minkowski {}.new_rep(4).slot(Aind::new_dummy());
        let ose = function!(
            self.ose,
            eidc,
            GS.emr_vec(e),
            m2,
            (m2 - self.emr_vec_index(e, mink.to_atom()) * self.emr_vec_index(e, mink.to_atom()))
        )
        .npow((1, 2));

        if let Some(index) = index {
            ose * self.energy_delta(index)
        } else {
            ose
        }
    }

    pub(crate) fn split_mom_pattern(&self, e: EdgeIndex, e_mass: Atom) -> Replacement {
        let id = Minkowski {}.to_symbolic([W_.a__]);
        Replacement::new(
            self.emr_mom(e, &id).to_pattern(),
            self.emr_vec_index(e, &id) + self.ose_full(e, e_mass, Some(id)) * sign_atom(e),
        )
    }

    pub(crate) fn split_mom_pattern_simple(&self, e: EdgeIndex) -> Replacement {
        let eidc = usize::from(e) as i64;
        Replacement::new(
            self.emr_mom(e, Minkowski {}.to_symbolic([W_.a__]))
                .to_pattern(),
            function!(GS.ose, eidc, Minkowski {}.to_symbolic([W_.a__])) * sign_atom(e)
                + function!(GS.emr_vec, eidc, Minkowski {}.to_symbolic([W_.a__])),
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
