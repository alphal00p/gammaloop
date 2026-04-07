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
    macro_rules! vk_parse {
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

    macro_rules! vk_symbol {
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

    pub(crate) use {vk_parse, vk_symbol};
}

#[macro_export]
macro_rules! vakint_parse {
    ($s: expr) => {{ symbolica::try_parse!($s, default_namespace = $crate::NAMESPACE) }};
    ($s: expr, $ns: expr) => {{ symbolica::try_parse!($sm, default_namespace = $ns) }};
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
    let mut res = replace_until_stable(
        input,
        &vk_parse!("x_^y_").unwrap().to_pattern(),
        &vk_parse!("1").unwrap().to_pattern(),
        Some(&Condition::from((crate::vk_symbol!("y_"), eq_condition(0)))),
        None,
    );
    res = replace_until_stable(
        res.as_view(),
        &vk_parse!("log(x__*exp(EulerGamma))").unwrap().to_pattern(),
        &vk_parse!("EulerGamma+log(x__)").unwrap().to_pattern(),
        None,
        None,
    );
    res = replace_until_stable(
        res.as_view(),
        &vk_parse!("exp(x_)^p_").unwrap().to_pattern(),
        &vk_parse!("exp(p_*x_)").unwrap().to_pattern(),
        None,
        None,
    );
    // res = replace_until_stable(
    //     res.as_view(),
    //     &vk_parse!("log(x__*𝜋)").unwrap().to_pattern(),
    //     &vk_parse!("log(x__)+log(𝜋)").unwrap().to_pattern(),
    //     None,
    //     None,
    // );
    // res = replace_until_stable(
    //     res.as_view(),
    //     &vk_parse!("log(x__*𝜋^p_)").unwrap().to_pattern(),
    //     &vk_parse!("log(x__)+p_*log(𝜋)").unwrap().to_pattern(),
    //     None,
    //     None,
    // );
    res = replace_until_stable(
        res.as_view(),
        &vk_parse!("log(exp(x_))").unwrap().to_pattern(),
        &vk_parse!("x_").unwrap().to_pattern(),
        None,
        None,
    );
    res = replace_until_stable(
        res.as_view(),
        &vk_parse!("exp(log(x_))").unwrap().to_pattern(),
        &vk_parse!("x_").unwrap().to_pattern(),
        None,
        None,
    );
    res = replace_until_stable(
        res.as_view(),
        &vk_parse!("log(x_^p_)").unwrap().to_pattern(),
        &vk_parse!("p_*log(x_)").unwrap().to_pattern(),
        None,
        None,
    );
    res = replace_until_stable(
        res.as_view(),
        &vk_parse!("exp(x_)*exp(y_)").unwrap().to_pattern(),
        &vk_parse!("exp(x_+y_)").unwrap().to_pattern(),
        None,
        None,
    );
    // Collect in EulerGamma to explicitly realise the MSbar cancellation
    res = res.collect::<i8>(vakint_parse!("EulerGamma").unwrap(), None, None);
    res
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
                truncated_float_re.set_prec(binary_prec);
                let mut truncated_float_im = fl_im.to_float();
                truncated_float_im.set_prec(binary_prec);
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
