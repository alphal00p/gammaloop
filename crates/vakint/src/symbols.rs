use std::{collections::HashSet, sync::LazyLock};

use crate::utils::vakint_macros::vk_symbol;
use regex::Regex;
use symbolica::{
    atom::{Atom, AtomOrView, Symbol},
    function, symbol,
};

pub static METRIC_SYMBOL: &str = "g";
pub static LOOP_MOMENTUM_SYMBOL: &str = "k";
pub static EXTERNAL_MOMENTUM_SYMBOL: &str = "p";
pub static DOT_SYMBOL: &str = "dot";
pub static MOMENTUM_WITH_INDEX_RE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(&format!(
        r"^(?:{}::)?(?:{}|{})\d+$",
        regex::escape(crate::NAMESPACE),
        regex::escape(LOOP_MOMENTUM_SYMBOL),
        regex::escape(EXTERNAL_MOMENTUM_SYMBOL),
    ))
    .expect("invalid momentum-with-index regex")
});

#[allow(dead_code)]
pub struct VakintSymbols {
    pub uedge: Symbol,
    pub dot: Symbol,
    pub dot_pow: Symbol,
    pub vkdot: Symbol,
    pub g: Symbol,
    pub g_form: Symbol,
    pub x: Symbol,
    pub y: Symbol,
    pub xa: Atom,
    pub ya: Atom,
    pub x_: Symbol,
    pub y_: Symbol,
    pub x_a: Atom,
    pub y_a: Atom,
    pub one: Atom,
    pub zero: Atom,
    pub error_flag_symbol: Symbol,
    pub error_flag: Atom,
    pub n_loops: Atom,
    pub p: Symbol,
    pub k: Symbol,
    pub id1_: Symbol,
    pub id2_: Symbol,
    pub id1_a: Atom,
    pub id2_a: Atom,
    pub pow: Symbol,
    pub pow_: Symbol,
    pub fun_: Symbol,
    pub fun_a: Atom,
    pub any___: Symbol,
    pub any_a___: Atom,
    pub n_: Symbol,
    pub a_: Symbol,
    pub b_: Symbol,
    pub c_: Symbol,
    pub lambda: Symbol,
    pub cmplx_i: Symbol,
    pub lambda_a: Atom,
    pub prop: Symbol,
    pub edge: Symbol,
    pub mom: Symbol,
    pub dot_dummy_ind: Symbol,
    pub topo: Symbol,
    pub metric: Symbol,
    pub uvprop: Symbol,
    pub vxs: Symbol,
    pub vec: Symbol,
    pub vec1: Symbol,
    pub form_epsilon: Symbol,
    pub form_dimension: Symbol,
    pub float_marker: Symbol,
}

pub static S: LazyLock<VakintSymbols> = LazyLock::new(|| VakintSymbols {
    uedge: symbol!(format!("{}::uedge",crate::NAMESPACE); Symmetric),
    dot: symbol!(
        format!("{}::dot",crate::NAMESPACE); Symmetric,Linear
    ),
    dot_dummy_ind: vk_symbol!("dot_dummy_ind"),
    dot_pow: vk_symbol!("dot_pow"),
    vkdot: symbol!(
        format!("{}::vkdot",crate::NAMESPACE);  Symmetric, Linear
    ),
    g: symbol!(format!("{}::{}",crate::NAMESPACE,METRIC_SYMBOL); Symmetric),
    g_form: symbol!(format!("{}::g",crate::NAMESPACE); Symmetric),
    x: vk_symbol!("x"),
    y: vk_symbol!("y"),
    xa: Atom::var(vk_symbol!("xa")),
    ya: Atom::var(vk_symbol!("ya")),
    x_: vk_symbol!("x_"),
    y_: vk_symbol!("y_"),
    x_a: Atom::var(vk_symbol!("x_")),
    y_a: Atom::var(vk_symbol!("y_")),
    one: Atom::num(1),
    zero: Atom::Zero,
    error_flag_symbol: vk_symbol!("ERROR"),
    error_flag: Atom::var(vk_symbol!("ERROR")),
    n_loops: Atom::var(vk_symbol!("n_loops")),
    p: vk_symbol!(EXTERNAL_MOMENTUM_SYMBOL),
    k: vk_symbol!(LOOP_MOMENTUM_SYMBOL),
    id1_: vk_symbol!("id1_"),
    id2_: vk_symbol!("id2_"),
    id1_a: Atom::var(vk_symbol!("id1_")),
    id2_a: Atom::var(vk_symbol!("id2_")),
    pow: vk_symbol!("pow"),
    pow_: vk_symbol!("pow_"),
    fun_: vk_symbol!("fun_"),
    fun_a: Atom::var(vk_symbol!("fun_")),
    any___: vk_symbol!("any___"),
    any_a___: Atom::var(vk_symbol!("any___")),
    n_: vk_symbol!("n_"),
    a_: vk_symbol!("a_"),
    b_: vk_symbol!("b_"),
    c_: vk_symbol!("c_"),
    cmplx_i: vk_symbol!("𝑖"),
    lambda: vk_symbol!("VakintLambdaScalingAnalysis"),
    lambda_a: Atom::var(vk_symbol!("VakintLambdaScalingAnalysis")),
    prop: vk_symbol!("prop"),
    edge: vk_symbol!("edge"),
    mom: vk_symbol!("k"),
    topo: vk_symbol!("topo"),
    metric: vk_symbol!("g"),
    uvprop: vk_symbol!("uvprop"),
    vxs: vk_symbol!("vxs"),
    vec: vk_symbol!("vec"),
    vec1: vk_symbol!("vec1"),
    form_epsilon: vk_symbol!("ep"),
    form_dimension: vk_symbol!("d"),
    float_marker: vk_symbol!("float"),
});

pub static SYMBOL_REGISTRY: LazyLock<HashSet<Symbol>> = LazyLock::new(|| {
    let s = &*S;

    [
        s.uedge,
        s.dot,
        s.dot_pow,
        s.vkdot,
        s.g,
        s.g_form,
        s.x,
        s.y,
        s.x_,
        s.y_,
        s.error_flag_symbol,
        s.p,
        s.k,
        s.id1_,
        s.id2_,
        s.pow,
        s.pow_,
        s.fun_,
        s.any___,
        s.n_,
        s.a_,
        s.b_,
        s.c_,
        s.lambda,
        s.cmplx_i,
        s.prop,
        s.edge,
        s.mom,
        s.dot_dummy_ind,
        s.topo,
        s.metric,
        s.uvprop,
        s.vxs,
        s.vec,
        s.vec1,
        s.form_epsilon,
        s.form_dimension,
        s.float_marker,
    ]
    .into_iter()
    .collect()
});
impl VakintSymbols {
    pub fn should_symbol_be_escaped_in_form(&self, symbol: &Symbol) -> bool {
        symbol.get_namespace() != crate::NAMESPACE
            || (!SYMBOL_REGISTRY.contains(symbol)
                && !MOMENTUM_WITH_INDEX_RE.is_match(symbol.get_name()))
    }

    // pub fn should_symbol_be_escaped_in_form(&self, symbol: &Symbol) -> bool {
    //     return symbol.get_namespace() != crate::NAMESPACE;
    // }

    pub fn dot<'a, A: Into<AtomOrView<'a>>>(&self, a: A, b: A) -> Atom {
        let a = a.into();
        let b = b.into();
        function!(self.dot, a.as_view(), b.as_view())
    }

    pub fn dot_pow<
        'a,
        A: Into<AtomOrView<'a>>,
        B: Into<AtomOrView<'a>>,
        C: Into<AtomOrView<'a>>,
    >(
        &self,
        a: A,
        b: B,
        pow: C,
    ) -> Atom {
        let a = a.into();
        let b = b.into();
        let pow = pow.into();
        function!(self.dot_pow, a.as_view(), b.as_view(), pow.as_view())
    }

    pub fn dot_dummy_ind<'a, A: Into<AtomOrView<'a>>>(&self, a: A) -> Atom {
        let a = a.into();
        function!(self.dot_dummy_ind, a.as_view())
    }
}
