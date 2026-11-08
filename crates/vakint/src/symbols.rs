use std::{collections::HashSet, sync::LazyLock};

use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use regex::Regex;
use spenso::network::tags::SPENSO_TAG;
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, Symbol},
    function,
    id::Pattern,
    symbol,
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

// These names identify recorded coefficients, not independent decimal tables.
pub(crate) const PERIOD_BINDINGS: [(&str, &str, i64); 5] = [
    ("PR4fin", "PR4", 0),
    ("PR4dfin", "PR4d", 0),
    ("PR11dep0", "PR11d", 0),
    ("PR4ep1", "PR4", 1),
    ("PR4dep1", "PR4d", 1),
];

#[allow(dead_code)]
pub struct VakintSymbols {
    pub uedge: Symbol,
    pub dot: Symbol,
    pub dot_pow: Symbol,
    pub vkdot: Symbol,
    /// Euclidean scalar product inside constant definitions; distinct from `dot`.
    pub euclidean_dot: Symbol,
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
    // Fixed parser-resolved vocabulary. These retain the aliases, attributes,
    // tags and user data of the original vk_symbol! lookups. Initialize them
    // once with the rest of Vakint's vocabulary, never inside an atom visitor.
    pub pg: Symbol,
    pub p_master: Symbol,
    pub p_gamma: Symbol,
    pub numerical_source: Symbol,
    pub rustred_master: Symbol,
    pub rustred_coefficient: Symbol,
    pub rustred_finite_literal: Symbol,
    pub rustred_constant: Symbol,
    pub rustred_input: Symbol,
    pub input_coefficient: Symbol,
    pub hr: Symbol,
    pub hi: Symbol,
    pub gpl_real: Symbol,
    pub gpl_imag: Symbol,
    pub polylog_real: Symbol,
    pub polylog_imag: Symbol,
    pub sqrt3: Symbol,
    pub cl2: Symbol,
    pub s2: Symbol,
    pub zeta: Symbol,
    pub sqrt: Symbol,
    pub om1: Symbol,
    pub euler_gamma: Symbol,
    pub d6: Symbol,
    pub t1ep: Symbol,
    pub t1ep2: Symbol,
    pub cl4: Symbol,
    pub s6: Symbol,
    pub mzv: Symbol,
    pub q_: Symbol,
    pub muvsq_: Symbol,
    pub nl_: Symbol,
    pub nr_: Symbol,
    pub props_: Symbol,
    pub idx___: Symbol,
    pub zeta_shorthands: [Symbol; 5],
    pub named_constants: [Symbol; 12],
    pub matad_bars: [(Symbol, Symbol, i64); 6],
    pub period_bindings: [(Symbol, Symbol, i64); 5],
    pub sixth_roots: [Atom; 4],
    pub inverse_sqrt3: Atom,
    pub rationalized_inverse_sqrt3: Atom,
    pub imaginary_inverse_sqrt3: Atom,
    pub rationalized_imaginary_inverse_sqrt3: Atom,
    pub epsilon_dimension: Atom,
    pub rustred_dimension: Atom,
    pub real_simplification_rules: [(Pattern, Pattern); 7],
    pub partial_numerical_input: Symbol,
    pub materialize_spectator: Symbol,
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
    euclidean_dot: symbol!(
        format!("{}::EuclideanDot",crate::NAMESPACE); Symmetric, Linear
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
    p: SPENSO_TAG.rank_one_tensor_symbol(&format!(
        "{}::{}",
        crate::NAMESPACE,
        EXTERNAL_MOMENTUM_SYMBOL
    )),
    k: SPENSO_TAG.rank_one_tensor_symbol(&format!(
        "{}::{}",
        crate::NAMESPACE,
        LOOP_MOMENTUM_SYMBOL
    )),
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
    // Unqualified 𝑖 is a numeric literal in Symbolica. Vakint keeps its
    // namespaced placeholder until the numerical/evaluator boundary.
    cmplx_i: symbol!(format!("{}::𝑖", crate::NAMESPACE)),
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
    float_marker: vk_symbol!("vakintFloat"),
    pg: vk_symbol!("PG"),
    p_master: vk_symbol!("PMaster"),
    p_gamma: vk_symbol!("PGamma"),
    numerical_source: vk_symbol!("NumericalSource"),
    rustred_master: vk_symbol!("RustRedMaster"),
    rustred_coefficient: vk_symbol!("RustRedCoefficient"),
    rustred_finite_literal: vk_symbol!("RustRedFiniteLiteral"),
    rustred_constant: vk_symbol!("RustRedConstant"),
    rustred_input: vk_symbol!("rustred_input"),
    input_coefficient: vk_symbol!("input_coefficient"),
    hr: vk_symbol!("HR"),
    hi: vk_symbol!("HI"),
    gpl_real: vk_symbol!("GPLReal"),
    gpl_imag: vk_symbol!("GPLImag"),
    polylog_real: vk_symbol!("PolylogReal"),
    polylog_imag: vk_symbol!("PolylogImag"),
    sqrt3: vk_symbol!("sqrt3"),
    cl2: vk_symbol!("cl2"),
    s2: vk_symbol!("S2"),
    zeta: vk_symbol!("Zeta"),
    sqrt: vk_symbol!("Sqrt"),
    om1: vk_symbol!("om1"),
    euler_gamma: vk_symbol!("EulerGamma"),
    d6: vk_symbol!("D6"),
    t1ep: vk_symbol!("T1ep"),
    t1ep2: vk_symbol!("T1ep2"),
    cl4: vk_symbol!("cl4"),
    s6: vk_symbol!("s6"),
    mzv: vk_symbol!("MZV"),
    q_: vk_symbol!("q_"),
    muvsq_: vk_symbol!("mUVsq_"),
    nl_: vk_symbol!("nl_"),
    nr_: vk_symbol!("nr_"),
    props_: vk_symbol!("props_"),
    idx___: vk_symbol!("idx___"),
    zeta_shorthands: ["z2", "z3", "z4", "z5", "z6"].map(|name| vk_symbol!(name)),
    named_constants: crate::inspection::constants::NAMED_CONSTANTS.map(|name| vk_symbol!(name)),
    matad_bars: crate::inspection::constants::MATAD_BARS
        .map(|(name, family, power)| (vk_symbol!(name), vk_symbol!(family), power)),
    period_bindings: PERIOD_BINDINGS
        .map(|(name, basis, power)| (vk_symbol!(name), vk_symbol!(basis), power)),
    sixth_roots: [
        "(1+3^(1/2)*𝑖)/2",
        "(-1+3^(1/2)*𝑖)/2",
        "(-1-3^(1/2)*𝑖)/2",
        "(1-3^(1/2)*𝑖)/2",
    ]
    .map(|value| vk_parse!(value).unwrap()),
    inverse_sqrt3: vk_parse!("(1/3)^(1/2)").unwrap(),
    rationalized_inverse_sqrt3: vk_parse!("3^(1/2)/3").unwrap(),
    imaginary_inverse_sqrt3: vk_parse!("(1/3)^(1/2)*𝑖").unwrap(),
    rationalized_imaginary_inverse_sqrt3: vk_parse!("3^(1/2)*𝑖/3").unwrap(),
    epsilon_dimension: vk_parse!("4-2*ep").unwrap(),
    rustred_dimension: vk_parse!("rustred::d").unwrap(),
    real_simplification_rules: [
        ("x_^y_", "1"),
        ("log(x__*exp(EulerGamma))", "EulerGamma+log(x__)"),
        ("exp(x_)^p_", "exp(p_*x_)"),
        ("log(exp(x_))", "x_"),
        ("exp(log(x_))", "x_"),
        ("log(x_^p_)", "p_*log(x_)"),
        ("exp(x_)*exp(y_)", "exp(x_+y_)"),
    ]
    .map(|(left, right)| {
        (
            vk_parse!(left).unwrap().to_pattern(),
            vk_parse!(right).unwrap().to_pattern(),
        )
    }),
    partial_numerical_input: vk_symbol!("partial_numerical_input"),
    materialize_spectator: vk_symbol!("RustRedMaterializeSpectator"),
});

pub static SYMBOL_REGISTRY: LazyLock<HashSet<Symbol>> = LazyLock::new(|| {
    let s = &*S;

    [
        s.uedge,
        s.dot,
        s.dot_pow,
        s.vkdot,
        s.euclidean_dot,
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

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::atom::{NamespacedSymbol, SymbolBuilder, UserData};

    #[test]
    fn global_parser_symbols_preserve_aliases_and_warmed_visitors() {
        const CHILD: &str = "VAKINT_GLOBAL_PARSER_SYMBOLS_TEST";
        if std::env::var_os(CHILD).is_none() {
            let output = std::process::Command::new(std::env::current_exe().unwrap())
                .args([
                    "--exact",
                    "symbols::tests::global_parser_symbols_preserve_aliases_and_warmed_visitors",
                    "--nocapture",
                ])
                .env(CHILD, "1")
                .output()
                .unwrap();
            assert!(
                output.status.success()
                    && String::from_utf8_lossy(&output.stdout).contains("1 passed"),
                "{}\n{}",
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr)
            );
            return;
        }
        // Initialize the real parser aliases before S, as a caller importing
        // symbols with metadata may do. A fresh default SymbolBuilder for the
        // short names would incorrectly lose these identities or conflict.
        let names = [
            "PG",
            "PMaster",
            "PGamma",
            "NumericalSource",
            "RustRedConstant",
            "PR4fin",
            "PR4",
        ];
        let mut aliases = Vec::new();
        for (index, name) in names.iter().enumerate() {
            let data = UserData::Atom(Atom::num(index + 11));
            let symbol = SymbolBuilder::new(NamespacedSymbol::parse(&format!(
                "vakint::ParserIdentity{index}"
            )))
            .with_aliases([*name])
            .with_tags(["vakint_test::parser_identity"])
            .with_user_data(data.clone())
            .build()
            .unwrap();
            aliases.push((symbol, data));
        }
        let globals = &*S;
        let resolved = [
            globals.pg,
            globals.p_master,
            globals.p_gamma,
            globals.numerical_source,
            globals.rustred_constant,
            globals.period_bindings[0].0,
            globals.period_bindings[0].1,
        ];
        for ((actual, name), (expected, data)) in resolved.into_iter().zip(names).zip(&aliases) {
            assert_eq!(actual, *expected);
            assert_eq!(actual, vk_symbol!(name));
            assert_eq!(actual, vk_symbol!(format!("vakint::{name}")));
            assert_eq!(actual.get_data(), data);
            assert!(actual.has_tag("vakint_test::parser_identity"));
        }
        let source = globals.numerical_source.call((1, 2, 3));
        let foreign = vk_symbol!("foreign_sources::NumericalSource").call((1, 2, 3));
        let reference = globals.rustred_constant.call((1, 2));
        let period = Atom::var(globals.period_bindings[0].0);
        let check = || {
            assert_eq!(crate::p_integral::pg_symbol(), globals.pg);
            assert_eq!(crate::p_integral::p_master_symbol(), globals.p_master);
            assert_eq!(crate::p_integral::p_gamma_symbol(), globals.p_gamma);
            assert!(crate::precision::legacy::is_source(source.as_view()));
            assert!(!crate::precision::legacy::is_source(foreign.as_view()));
            assert!(!crate::precision::legacy::is_source(
                Atom::var(globals.numerical_source).as_view()
            ));
            assert!(crate::precision::snapshot_reference(reference.as_view()));
            assert!(
                crate::rustred_evaluation::masters::MasterCatalog::is_recorded_period(
                    period.as_view()
                )
            );
            assert!(
                !crate::rustred_evaluation::masters::MasterCatalog::is_recorded_period(
                    source.as_view()
                )
            );
        };
        check();
        let before = symbolica::state::State::symbol_iter().count();
        std::thread::scope(|scope| {
            for _ in 0..8 {
                scope.spawn(|| {
                    for _ in 0..256 {
                        check();
                    }
                });
            }
        });
        assert_eq!(
            symbolica::state::State::symbol_iter().count(),
            before,
            "warmed fixed-vocabulary traversals must not create symbols"
        );
        // Registry writes themselves are private in Symbolica; these checks
        // cover observable identity/metadata and concurrent lookup behavior.
    }

    #[test]
    fn global_fixed_atoms_and_rewrite_rules_match_the_original_parser() {
        let globals = &*S;
        assert_eq!(Atom::i(), vk_parse!("𝑖").unwrap());
        for (actual, source) in globals.sixth_roots.iter().zip([
            "(1+3^(1/2)*𝑖)/2",
            "(-1+3^(1/2)*𝑖)/2",
            "(-1-3^(1/2)*𝑖)/2",
            "(1-3^(1/2)*𝑖)/2",
        ]) {
            assert_eq!(actual, &vk_parse!(source).unwrap());
        }
        for (actual, source) in [
            (&globals.inverse_sqrt3, "(1/3)^(1/2)"),
            (&globals.rationalized_inverse_sqrt3, "3^(1/2)/3"),
            (&globals.imaginary_inverse_sqrt3, "(1/3)^(1/2)*𝑖"),
            (&globals.rationalized_imaginary_inverse_sqrt3, "3^(1/2)*𝑖/3"),
            (&globals.epsilon_dimension, "4-2*ep"),
            (&globals.rustred_dimension, "rustred::d"),
        ] {
            assert_eq!(actual, &vk_parse!(source).unwrap());
        }
        // Exercise every cached rewrite rule against its previous text-parsed
        // implementation, preserving order and the conditional zero exponent.
        for input in [
            "x^0+log(x*exp(EulerGamma))",
            "exp(x)^p+log(exp(y))",
            "exp(log(x))+log(x^p)+exp(x)*exp(y)",
        ] {
            let input = vk_parse!(input).unwrap();
            let mut expected = crate::utils::replace_until_stable(
                input.as_view(),
                &vk_parse!("x_^y_").unwrap().to_pattern(),
                &vk_parse!("1").unwrap().to_pattern(),
                Some(&symbolica::id::Condition::from((
                    vk_symbol!("y_"),
                    crate::eq_condition(0),
                ))),
                None,
            );
            for (left, right) in [
                ("log(x__*exp(EulerGamma))", "EulerGamma+log(x__)"),
                ("exp(x_)^p_", "exp(p_*x_)"),
                ("log(exp(x_))", "x_"),
                ("exp(log(x_))", "x_"),
                ("log(x_^p_)", "p_*log(x_)"),
                ("exp(x_)*exp(y_)", "exp(x_+y_)"),
            ] {
                expected = crate::utils::replace_until_stable(
                    expected.as_view(),
                    &vk_parse!(left).unwrap().to_pattern(),
                    &vk_parse!(right).unwrap().to_pattern(),
                    None,
                    None,
                );
            }
            expected = expected.collect::<i8>(vk_parse!("EulerGamma").unwrap());
            assert_eq!(crate::utils::simplify_real(input.as_view()), expected);
        }
    }
}
