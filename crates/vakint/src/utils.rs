use ahash::HashMap;
use regex::Regex;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::CoefficientView,
    domains::float::Complex,
    id::{Condition, MatchSettings, Pattern, PatternRestriction},
};
use vakint_macros::vk_parse;

use crate::{
    PYSECDEC_ATTRIBUTE_END, PYSECDEC_ATTRIBUTE_SEPARATOR, PYSECDEC_ATTRIBUTE_START,
    PYSECDEC_NAMESPACE_SEPARATOR, PYSECDEC_UNDERSCORE, VAKINT_SYMBOL_UNDRESSING_RE, VakintError,
    VakintSettings, eq_condition,
};

pub(crate) mod vakint_macros {
    /// Parse `$s` in the vakint namespace. A string literal is parsed once per
    /// call site and cloned afterwards; symbols never change once created.
    macro_rules! vk_parse {
        ($s: literal) => {{
            static PARSED: std::sync::LazyLock<Result<symbolica::atom::Atom, String>> =
                std::sync::LazyLock::new(|| symbolica::try_parse!($s, default_namespace = crate::NAMESPACE));
            PARSED.clone()
        }};
        ($s: expr) => {{ symbolica::try_parse!($s, default_namespace = crate::NAMESPACE) }};
        ($s: expr, $ns: expr) => {{ symbolica::try_parse!($s, default_namespace = $ns) }};
    }
    // macro_rules! vk_symbol {
    //     ($s: expr) => {{
    //         if format!("{}", $s).starts_with(format!("{}::", crate::NAMESPACE).as_str()) {
    //             symbolica::symbol!($s)
    //         } else {
    //             symbolica::symbol!(format!("{}::{}", crate::NAMESPACE, $s))
    //         }
    //     }};
    // }

    /// The symbol `$s` in the vakint namespace. A string literal is parsed once
    /// per call site: hot visitors and predicates name fixed symbols, and the
    /// parser (with two formatted strings) dominated them.
    macro_rules! vk_symbol {
        ($s: literal) => {{
            static SYMBOL: std::sync::LazyLock<symbolica::atom::Symbol> =
                std::sync::LazyLock::new(|| $crate::utils::vakint_macros::vk_symbol_uncached!($s));
            *SYMBOL
        }};
        ($s: expr) => {{ $crate::utils::vakint_macros::vk_symbol_uncached!($s) }};
    }

    macro_rules! vk_symbol_uncached {
        ($s: expr) => {{
            if format!("{}", $s).starts_with(format!("{}::", crate::NAMESPACE).as_str()) {
                match symbolica::try_parse!($s) {
                    Ok::<Atom, String>(a) => match a {
                        Atom::Var(sym) => sym.get_symbol(),
                        _ => panic!("Parsed atom is not a symbol: {}.", $s),
                    },
                    Err(e) => panic!("Could not parse symbol from string: {}, Error: {}", $s, e),
                }
            } else {
                match symbolica::try_parse!(format!("{}", $s), default_namespace = crate::NAMESPACE)
                {
                    Ok::<Atom, String>(a) => match a {
                        Atom::Var(sym) => sym.get_symbol(),
                        _ => panic!("Parsed atom is not a symbol: {}.", $s),
                    },
                    Err(e) => panic!("Could not parse symbol from string: {}, Error: {}", $s, e),
                }
            }
        }};
    }

    pub(crate) use {vk_parse, vk_symbol, vk_symbol_uncached};
}

#[macro_export]
macro_rules! vakint_parse {
    ($s: expr) => {{
        $crate::Vakint::initialize_vakint_symbols();
        symbolica::try_parse!($s, default_namespace = $crate::NAMESPACE)
    }};
    ($s: expr, $ns: expr) => {{
        $crate::Vakint::initialize_vakint_symbols();
        symbolica::try_parse!($s, default_namespace = $ns)
    }};
}

// #[macro_export]
// macro_rules! vakint_symbol {
//     ($s: expr) => {{
//         if format!("{}", $s).starts_with(format!("{}::", $crate::NAMESPACE).as_str()) {
//             symbolica::symbol!($s)
//         } else {
//             symbolica::symbol!(format!("{}::{}", $crate::NAMESPACE, $s))
//         }
//     }};
// }

#[macro_export]
macro_rules! vakint_symbol {
    ($s: expr) => {{
        $crate::Vakint::initialize_vakint_symbols();
        if format!("{}", $s).starts_with(format!("{}::", $crate::NAMESPACE).as_str()) {
            match symbolica::try_parse!($s) {
                Ok::<_, String>(a) => match a {
                    symbolica::atom::Atom::Var(sym) => sym.get_symbol(),
                    _ => panic!("Parsed atom is not a symbol: {}.", $s),
                },
                Err(e) => panic!("Could not parse symbol from string: {}, Error: {}", $s, e),
            }
        } else {
            match symbolica::try_parse!(format!("{}", $s), default_namespace = $crate::NAMESPACE) {
                Ok::<_, String>(a) => match a {
                    symbolica::atom::Atom::Var(sym) => sym.get_symbol(),
                    _ => panic!("Parsed atom is not a symbol: {}.", $s),
                },
                Err(e) => panic!("Could not parse symbol from string: {}, Error: {}", $s, e),
            }
        }
    }};
}

pub fn replace_until_stable(
    target: AtomView<'_>,
    lhs: &Pattern,
    rhs: &Pattern,
    conditions: Option<&Condition<PatternRestriction>>,
    settings: Option<&MatchSettings>,
) -> Atom {
    // Specifiying settings is not possible anymore given new Symbolica v0.16 API
    assert!(settings.is_none());

    // We do not want to use `repeat()` below because it does not use the same termination condition
    // which we want in our case to be that the atom is stable.

    // let mut res = target.replace(lhs);
    // if let Some(c) = conditions {
    //     res = res.when(c);
    // }
    // let a = res.repeat().with(rhs);
    // a
    let mut replacer = target.replace(lhs);
    if let Some(c) = conditions {
        replacer = replacer.when(c);
    }
    let mut old_res = replacer.with(rhs);
    loop {
        // println!("lhs={}, rhs={}", lhs, rhs);
        let mut replacer = old_res.replace(lhs);
        if let Some(c) = conditions {
            replacer = replacer.when(c);
        }
        let res = replacer.with(rhs);
        if res == old_res {
            break;
        }
        old_res = res;
    }
    old_res
}

pub fn simplify_real(input: AtomView) -> Atom {
    let symbols = &*crate::symbols::S;
    let (left, right) = &symbols.real_simplification_rules[0];
    let mut res = replace_until_stable(
        input,
        left,
        right,
        Some(&Condition::from((symbols.y_, eq_condition(0)))),
        None,
    );
    for (left, right) in &symbols.real_simplification_rules[1..] {
        res = replace_until_stable(res.as_view(), left, right, None, None);
    }
    // With no EulerGamma, exact collection only rebuilds the same coefficient.
    // Powers of ordinary variables recombine after splitting exponent sums.
    // Keep collection for other bases (including exp/log normalization on E)
    // and floating/embedded coefficients, where the result can change.
    if !res.contains_symbol(symbols.euler_gamma) {
        let mut unchanged_by_collection = true;
        res.visitor(&mut |part| {
            unchanged_by_collection &= match part {
                AtomView::Num(number) => matches!(
                    number.get_coeff_view(),
                    CoefficientView::Natural(..) | CoefficientView::Large(..)
                ),
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    !matches!(exponent, AtomView::Add(_))
                        || matches!(base, AtomView::Var(variable) if variable.get_symbol() != Symbol::E)
                }
                _ => true,
            };
            unchanged_by_collection
        });
        if unchanged_by_collection {
            return res;
        }
    }
    // Collect in EulerGamma to explicitly realise the MSbar cancellation.
    res.collect::<i8>(Atom::var(symbols.euler_gamma))
}

pub fn multi_string_replace(text: String, repls: &HashMap<String, String>) -> String {
    // 1. Collect and sort patterns by decreasing length
    let mut patterns: Vec<String> = repls.keys().cloned().collect();
    patterns.sort_unstable_by_key(|s| std::cmp::Reverse(s.len()));

    // 2. Build an alternation like "longer_pattern|shorter"
    let alternation = patterns
        .iter()
        .map(|pat| regex::escape(pat))
        .collect::<Vec<_>>()
        .join("|");
    let re = Regex::new(&alternation).unwrap();

    // 3. Replace all in one go, using the closure to look up the right replacement
    re.replace_all(&text, |caps: &regex::Captures| {
        let m = caps.get(0).unwrap().as_str();
        repls[m].clone()
    })
    .into_owned()
}

pub fn could_match(pattern: &Pattern, target: AtomView) -> bool {
    target
        .pattern_match(
            pattern,
            Some(&Condition::default()),
            Some(&MatchSettings::default()),
        )
        .next()
        .is_some()
}

pub fn get_full_name(symbol: &Symbol) -> String {
    Atom::var(*symbol).to_canonical_string()
}

pub fn undress_vakint_symbols(expression: &str) -> String {
    (*VAKINT_SYMBOL_UNDRESSING_RE)
        .replace_all(expression, "")
        .into_owned()
}

pub fn set_precision_in_float_atom(input: AtomView, settings: &VakintSettings) -> Atom {
    let binary_prec = settings.get_binary_precision();
    if let AtomView::Num(fl_view) = input {
        match fl_view.get_coeff_view() {
            CoefficientView::Float(fl_re, fl_im) => {
                let mut truncated_float_re = fl_re.to_float();
                truncated_float_re.set_prec(binary_prec.min(truncated_float_re.prec()));
                let mut truncated_float_im = fl_im.to_float();
                truncated_float_im.set_prec(binary_prec.min(truncated_float_im.prec()));
                Atom::num(Complex::new(truncated_float_re, truncated_float_im))
            }
            _ => input.to_owned(),
        }
    } else {
        input.to_owned()
    }
}

pub fn set_precision_in_polynomial_atom(
    input: AtomView,
    variable: Symbol,
    settings: &VakintSettings,
) -> Atom {
    let mut res = Atom::Zero;
    let coeffs = input.coefficient_list::<i8>(&[Atom::var(variable)]);
    for (var_with_power, coeff) in coeffs {
        let new_coeff = set_precision_in_float_atom(coeff.as_view(), settings);
        res += new_coeff * var_with_power;
    }

    res
}

pub fn to_vk_symbol(s: &str) -> Result<Symbol, VakintError> {
    match vk_parse!(s, crate::NAMESPACE) {
        Ok::<Atom, String>(a) => match a {
            Atom::Var(sym) => Ok(sym.get_symbol()),
            _ => Err(VakintError::InvalidGenericExpression(format!(
                "Parsed atom is not a symbol: {}",
                a
            ))),
        },
        Err(e) => Err(VakintError::InvalidGenericExpression(format!(
            "Could not parse symbol from string: {}",
            e
        ))),
    }
}

pub fn to_symbol(s: &str) -> Result<Symbol, VakintError> {
    match vk_parse!(s) {
        Ok::<Atom, String>(a) => match a {
            Atom::Var(sym) => Ok(sym.get_symbol()),
            _ => Err(VakintError::InvalidGenericExpression(format!(
                "Parsed atom is not a symbol: {}",
                a
            ))),
        },
        Err(e) => Err(VakintError::InvalidGenericExpression(format!(
            "Could not parse symbol from string: {}",
            e
        ))),
    }
}

pub fn pysecdec_encode(expression: &str) -> String {
    expression
        .replace("::", PYSECDEC_NAMESPACE_SEPARATOR)
        .replace("_", PYSECDEC_UNDERSCORE)
        .replace("{", PYSECDEC_ATTRIBUTE_START)
        .replace(",", PYSECDEC_ATTRIBUTE_SEPARATOR)
        .replace("}", PYSECDEC_ATTRIBUTE_END)
        .replace("^", "**")
}

pub fn pysecdec_decode(expression: &str) -> String {
    expression
        .replace(PYSECDEC_NAMESPACE_SEPARATOR, "::")
        .replace(PYSECDEC_UNDERSCORE, "_")
        .replace(PYSECDEC_ATTRIBUTE_START, "{")
        .replace(PYSECDEC_ATTRIBUTE_SEPARATOR, ",")
        .replace(PYSECDEC_ATTRIBUTE_END, "}")
        .replace("**", "^")
}

pub fn split_linear_atom(a: AtomView, variable: AtomView) -> (Atom, Atom) {
    let split_atom: HashMap<_, _> = a.coefficient_list::<i8>(&[variable]).into_iter().collect();
    if split_atom.len() > 2 {
        panic!(
            "Could not split linear atom '{}' into at most two parts using variable '{}'. Split atom:\n{}",
            a,
            variable,
            split_atom
                .iter()
                .map(|(k, v)| format!("{}: {}", k, v))
                .collect::<Vec<_>>()
                .join("\n")
        );
    } else {
        (
            split_atom
                .get(&variable.to_owned())
                .unwrap_or(&Atom::Zero)
                .to_owned(),
            split_atom
                .get(&Atom::num(1))
                .unwrap_or(&Atom::Zero)
                .to_owned(),
        )
    }
}

#[test]
fn test_set_precision_in_float_atom() {
    println!(
        "Res={}",
        set_precision_in_float_atom(
            vk_parse!("1.23456789113245423`100").unwrap().as_view(),
            &VakintSettings {
                run_time_decimal_precision: 30,
                ..VakintSettings::default()
            }
        )
    );
}

#[cfg(test)]
mod simplify_real_tests {
    use super::*;

    /// The previous complete simplifier, including unconditional collection.
    fn previous(input: AtomView) -> Atom {
        let symbols = &*crate::symbols::S;
        let (left, right) = &symbols.real_simplification_rules[0];
        let mut result = replace_until_stable(
            input,
            left,
            right,
            Some(&Condition::from((symbols.y_, eq_condition(0)))),
            None,
        );
        for (left, right) in &symbols.real_simplification_rules[1..] {
            result = replace_until_stable(result.as_view(), left, right, None, None);
        }
        result.collect::<i8>(Atom::var(symbols.euler_gamma))
    }

    #[test]
    fn exact_gamma_free_simplification_matches_previous_collector() {
        crate::Vakint::initialize_vakint_symbols();
        for input in [
            "0",
            "1",
            "x+y+z",
            "(x+y)*(a+b)",
            "(x+y)^2",
            "(x+y)^(4/3)",
            "x^(a+b)",
            "x^(2*a+4*b)",
            "(m^-2)^(3*ep)*mu^(3*ep)*(x/(1-ep)+y/(1+ep))",
            "(2+3*𝑖)*(x/(1-ep)+y/(1+ep))",
            "123456789012345678901234567890123456789*x/7+y",
            "opaque(x+y,(a+b)^2)*x^(1/2)",
            "exp(log(x))+log(exp(y))+log(x^2)",
        ] {
            let input = vk_parse!(input).unwrap();
            assert!(!input.contains_symbol(crate::symbols::S.euler_gamma));
            assert_eq!(
                simplify_real(input.as_view()),
                previous(input.as_view()),
                "{input}"
            );
        }
        let symbols = &*crate::symbols::S;
        let epsilon = Atom::var(symbols.form_epsilon);
        let factor = vk_parse!("opaque_exact_sum").unwrap();
        let input = Atom::add_many((1..=2048).map(|index| {
            &factor * symbols.rustred_master.call((index,)) / (Atom::num(index) - &epsilon)
        }));
        assert!(!input.contains_symbol(symbols.euler_gamma));
        assert_eq!(simplify_real(input.as_view()), input);
        assert_eq!(simplify_real(input.as_view()), previous(input.as_view()));
    }

    #[test]
    fn additive_mass_exponents_recombine_without_changing_the_result() {
        crate::Vakint::initialize_vakint_symbols();
        for input in [
            "mUV^(4-6*ep)",
            "mUV^(-2-3*ep)*mu^(6*ep)",
            "mUV^(4/3-2*ep+3*ep^2)",
            "mUV^(a/2+b/3+(2+3*𝑖)*ep)",
            "x^(2*a+4*b)",
            "𝜋^(1+a)",
        ] {
            let input = vk_parse!(input).unwrap();
            assert_eq!(previous(input.as_view()), input, "{input}");
            assert_eq!(simplify_real(input.as_view()), input, "{input}");
        }

        let mass_factor = vk_parse!("mUV^(4-6*ep)*mu^(6*ep-2)").unwrap();
        let epsilon = Atom::var(crate::symbols::S.form_epsilon);
        let input = Atom::add_many((1..=128).map(|index| {
            &mass_factor * crate::symbols::S.rustred_master.call((index,))
                / (Atom::num(index) - &epsilon)
        }));
        assert_eq!(previous(input.as_view()), input);
        assert_eq!(simplify_real(input.as_view()), input);
    }

    #[test]
    fn nonvariable_additive_powers_keep_the_collectors_structural_changes() {
        crate::Vakint::initialize_vakint_symbols();
        // Unlike an ordinary variable, a numeric or product base loses its
        // identity when the exponent-one factor is normalized. Keep that path.
        for (input, expected) in [("2^(1+a)", "2*2^a"), ("(x*y)^(1+a)", "x*y*(x*y)^a")] {
            let input = vk_parse!(input).unwrap();
            let expected = vk_parse!(expected).unwrap();
            assert_ne!(input, expected);
            assert_eq!(previous(input.as_view()), expected, "{input}");
            assert_eq!(simplify_real(input.as_view()), expected, "{input}");
        }

        // Built-in E is a variable with its own exp/log rules; i is a numeric
        // coefficient. Both retain collection when their exponent is a sum.
        let exponential = Atom::var(Symbol::E).pow(vk_parse!("x+y").unwrap());
        assert_eq!(
            simplify_real(exponential.as_view()),
            previous(exponential.as_view())
        );
        assert!(matches!(
            vk_parse!("𝑖").unwrap().as_view(),
            AtomView::Num(_)
        ));
        for input in ["exp(x+y)", "exp(log(x)+y)", "𝑖^(1+a)"] {
            let input = vk_parse!(input).unwrap();
            assert_eq!(
                simplify_real(input.as_view()),
                previous(input.as_view()),
                "{input}"
            );
        }
    }

    #[test]
    fn gamma_cancellation_and_hidden_occurrences_keep_previous_results() {
        crate::Vakint::initialize_vakint_symbols();
        for input in [
            "EulerGamma*x+EulerGamma*y+x",
            "(1+EulerGamma)^3-EulerGamma^3-3*EulerGamma^2-3*EulerGamma",
            "opaque(EulerGamma,x)+x*EulerGamma",
            "x^EulerGamma+opaque(x^(2*EulerGamma))",
            "EulerGamma(x)+x",
            "log(x*exp(EulerGamma))-EulerGamma",
        ] {
            let input = vk_parse!(input).unwrap();
            assert!(input.contains_symbol(crate::symbols::S.euler_gamma));
            assert_eq!(
                simplify_real(input.as_view()),
                previous(input.as_view()),
                "{input}"
            );
        }
        let normalization = vk_parse!("exp(-ep*EulerGamma)*exp(ep*EulerGamma)*(1+ep)").unwrap();
        assert_eq!(
            simplify_real(normalization.as_view()),
            vk_parse!("1+ep").unwrap()
        );
        assert_eq!(
            simplify_real(normalization.as_view()),
            previous(normalization.as_view())
        );
    }

    #[test]
    fn floating_simplification_preserves_previous_rounding_and_structure() {
        crate::Vakint::initialize_vakint_symbols();
        for input in [
            "0.1*x+0.2*x^2+0.3*y",
            "(0.1+x)^3+0.2*x",
            "(0.1*x+0.2*x^2)/(1-x)+0.3*y",
            "opaque(0.5)*x+0.3*y",
            "0.1*EulerGamma*x+0.2*EulerGamma*x",
        ] {
            let input = vk_parse!(input).unwrap();
            assert_eq!(
                simplify_real(input.as_view()),
                previous(input.as_view()),
                "{input}"
            );
        }
    }
}
