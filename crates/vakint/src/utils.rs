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
                std::sync::LazyLock::new(|| {
                    symbolica::try_parse!($s, default_namespace = crate::NAMESPACE)
                });
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
            // A symbol with a display form is created here, at the point where the
            // parse below would create it, so the creation order is unchanged.
            let full_name = format!("{}", $s);
            let full_name = if full_name.contains("::") {
                full_name
            } else {
                format!("{}::{}", crate::NAMESPACE, full_name)
            };
            if let Some(decorated) = crate::display::create_decorated(&full_name) {
                decorated
            } else if format!("{}", $s).starts_with(format!("{}::", crate::NAMESPACE).as_str()) {
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

/// Whether every number in `atom` is an exact (complex) rational.
pub(crate) fn has_only_exact_numbers(atom: AtomView) -> bool {
    let mut exact = true;
    atom.visitor(&mut |part| {
        if let AtomView::Num(number) = part {
            exact &= matches!(
                number.get_coeff_view(),
                CoefficientView::Natural(..) | CoefficientView::Large(..)
            );
        }
        exact
    });
    exact
}

/// The left fold `0 + t_1 + ... + t_n` of normalized atoms under Symbolica's
/// `+`, bit for bit, with one n-way merge ([`Atom::add_many`]) per run of
/// terms instead of one renormalization of the growing sum per term.
///
/// Both merge a term into the sum by adding its number coefficient to that
/// of the term with the same non-number part. They differ only where the
/// order of these merges shows:
/// - a coefficient that becomes exactly 1 on a single sum factor, as in
///   `c*(a+b)`, leaves the sum `a+b`, which `+` flattens into the running
///   sum at once, so later `c'*(a+b)` no longer merge with it:
///   `(a+b)/2 + (a+b)/2 + (a+b)/2 = a+b+1/2*(a+b)`, where one merge of the
///   three terms gives `3/2*(a+b)`;
/// - `+` keeps an identical sum factored, `(a+b) + (a+b) = 2*(a+b)`, which
///   requires the next term to be a sum equal to the running sum;
/// - inexact numbers may round differently in another merge order.
///
/// A term that is a sum, a number times one sum ([`is_scaled_sum`]) or
/// holds an inexact number is therefore added with `+` itself, after the run
/// before it is merged, and once the running sum holds an inexact number
/// every later term is added with `+`. The non-number part of any other term
/// is not a single sum, so no merge with it forms a sum (a coefficient 1
/// leaves that non-number part, 0 removes the term), and its exact
/// coefficients add up to the same number in any order.
pub(crate) struct ExactSum {
    total: Atom,
    exact_total: bool,
    pending: Vec<Atom>,
}

impl ExactSum {
    pub(crate) fn new() -> Self {
        Self {
            total: Atom::Zero,
            exact_total: true,
            pending: Vec::new(),
        }
    }

    pub(crate) fn push(&mut self, term: Atom) {
        if self.exact_total
            && !matches!(term.as_view(), AtomView::Add(_))
            && !is_scaled_sum(term.as_view())
            && has_only_exact_numbers(term.as_view())
        {
            self.pending.push(term);
            return;
        }
        self.flush();
        self.exact_total &= has_only_exact_numbers(term.as_view());
        self.total += term;
    }

    fn flush(&mut self) {
        if self.pending.is_empty() {
            return;
        }
        let total = std::mem::take(&mut self.total);
        self.total = Atom::add_many(std::iter::once(total).chain(self.pending.drain(..)));
    }

    pub(crate) fn finish(mut self) -> Atom {
        self.flush();
        self.total
    }
}

/// Whether `term` is a number times one sum, `c*(a+b)`: a term that a merge
/// of coefficients can turn into a sum ([`ExactSum`]).
pub(crate) fn is_scaled_sum(term: AtomView) -> bool {
    let AtomView::Mul(product) = term else {
        return false;
    };
    let mut factors = product
        .iter()
        .filter(|factor| !matches!(factor, AtomView::Num(_)));
    matches!(
        (factors.next(), factors.next()),
        (Some(AtomView::Add(_)), None)
    )
}

/// The expanded summands that [`expand_merging`] holds unmerged, in bytes,
/// before it merges them into its running sum (at least the size of the
/// running sum, so every merge at least doubles what it merges).
const EXPAND_MERGING_PENDING_BYTES: usize = 64 << 20;

/// `expression.expand()`, with the summands of a top-level sum expanded one
/// at a time and merged into a running sum by one n-way merge
/// ([`Atom::add_many`]) per batch of pending summands, so the unmerged
/// expansion of the whole sum is never held at once.
///
/// `expand` of a sum normalizes the concatenated unmerged expansions of its
/// summands once; here each summand's expansion is normalized and the
/// normalized pieces are merged. Both sort the same terms canonically and
/// add the number coefficients of the terms with the same non-number part,
/// which with exact numbers gives the same coefficients in any grouping (a
/// coefficient 1 leaves that non-number part, 0 removes the term). They
/// could differ only where a merged term becomes a sum, which needs a term
/// that is a number times one sum ([`is_scaled_sum`]) and which an
/// expansion never leaves; and with inexact numbers, which may round
/// differently. In both cases the whole expression is expanded at once
/// instead.
pub(crate) fn expand_merging(expression: AtomView) -> Atom {
    expand_merging_with(expression, EXPAND_MERGING_PENDING_BYTES)
}

/// [`expand_merging`] with a merge once more than `pending_bound` bytes
/// (and more than the running sum) are pending.
fn expand_merging_with(expression: AtomView, pending_bound: usize) -> Atom {
    let AtomView::Add(sum) = expression else {
        return expression.expand();
    };
    if !has_only_exact_numbers(expression) {
        return expression.expand();
    }
    let mut total = Atom::Zero;
    let mut total_bytes = 0;
    let mut pending = Vec::<Atom>::new();
    let mut pending_bytes = 0;
    let merge = |total: Atom, pending: &mut Vec<Atom>| {
        Atom::add_many(std::iter::once(total).chain(pending.drain(..)))
    };
    for summand in sum.iter() {
        let expanded = summand.expand();
        let has_scaled_sum = match expanded.as_view() {
            AtomView::Add(terms) => terms.iter().any(is_scaled_sum),
            term => is_scaled_sum(term),
        };
        if has_scaled_sum {
            return expression.expand();
        }
        pending_bytes += expanded.as_view().get_byte_size();
        pending.push(expanded);
        if pending_bytes > pending_bound.max(total_bytes) {
            total = merge(total, &mut pending);
            total_bytes = total.as_view().get_byte_size();
            pending_bytes = 0;
        }
    }
    merge(total, &mut pending)
}

#[cfg(test)]
mod expand_merging_tests {
    use symbolica::atom::AtomCore;

    use super::{expand_merging, expand_merging_with};
    use crate::utils::vakint_macros::vk_parse;

    fn assert_same(source: &str) {
        let expression = vk_parse!(source).unwrap();
        let expected = expression.expand();
        for actual in [
            expand_merging(expression.as_view()),
            expand_merging_with(expression.as_view(), 0),
        ] {
            assert!(
                actual == expected
                    && actual.to_canonical_string() == expected.to_canonical_string(),
                "{source}: expand {expected}, merged {actual}"
            );
        }
    }

    #[test]
    fn expand_merging_equals_expand() {
        for source in [
            "0",
            "x",
            "(x+y)^2",
            "(x+y)^2-(x-y)^2",
            "(x+y)^2-x^2-2*x*y-y^2",
            "(a+b)*(c+d)+(a-b)*(c-d)+2*a*c",
            "(a+b)/2+(a+b)/2",
            "x*(a+b)/2+x*(a+b)/2+(a+b)*x",
            "(1/2+1𝑖)*(x+y)^2+(1/2-1𝑖)*(x+y)^2-(x+y)^2",
            "f(x+y)*(x+1)+f(x+y)*(1-x)",
            "(x+y)^-1*(x+1)+(x+y)^-1*(x-1)+2^(1/2)*(x+y)^2+2^(1/2)*x",
            "(x^(1/3)+y)^3+x^(2/3)*(y+1)",
            "1.5*(x+y)+(x+y)^2",
            "(dot(k(1),p(1))+dot(k(2),p(1)))^3-dot(k(1),p(1))^3+3*(x+dot(k(1),k(2)))^2",
        ] {
            assert_same(source);
        }
    }

    /// Random sums of products of powers of sums, with cancellations and
    /// merges across summands.
    #[test]
    fn expand_merging_equals_expand_on_random_sums() {
        let pool = "x y -x 2 1/3 -5/7 2𝑖 (x+y) (x-y) (1+x) (y+z) (x+y)^2 (x-y)^2 (1+x)^3 \
            f(x+y) (x+y)^-1 x^(1/2) 2^(1/2) (2^(1/2)*x+1) (x+y)/2 (x^2+2*x*y+y^2) z^-2";
        let factors = pool.split_whitespace().collect::<Vec<_>>();
        let mut state = 0x9e37_79b9_7f4a_7c15_u64;
        let mut below = |n: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        };
        for _ in 0..3000 {
            let mut sum = String::new();
            for summand in 0..1 + below(8) {
                if summand > 0 {
                    sum.push_str(if below(3) == 0 { "-" } else { "+" });
                }
                for factor in 0..1 + below(4) {
                    if factor > 0 {
                        sum.push('*');
                    }
                    sum.push_str(factors[below(factors.len())]);
                }
            }
            assert_same(&sum);
        }
    }
}

#[cfg(test)]
mod exact_sum_tests {
    use symbolica::atom::{Atom, AtomCore};

    use super::ExactSum;
    use crate::utils::vakint_macros::vk_parse;

    fn fold(terms: &[Atom]) -> Atom {
        let mut total = Atom::zero();
        for term in terms {
            total += term.clone();
        }
        total
    }

    fn merged(terms: &[Atom]) -> Atom {
        let mut sum = ExactSum::new();
        for term in terms {
            sum.push(term.clone());
        }
        sum.finish()
    }

    fn assert_same(sources: &[&str]) {
        let terms = sources
            .iter()
            .map(|source| vk_parse!(source).unwrap())
            .collect::<Vec<_>>();
        let expected = fold(&terms);
        let actual = merged(&terms);
        assert_eq!(
            actual.to_canonical_string(),
            expected.to_canonical_string(),
            "{sources:?}"
        );
        assert_eq!(actual, expected, "{sources:?}");
    }

    #[test]
    fn exact_sum_reproduces_the_left_fold() {
        let cases: &[&[&str]] = &[
            &[],
            &["0"],
            &["x"],
            &["x", "y", "-x", "z^2", "2*y"],
            &["3", "1/2", "-7/2"],
            &["x", "-x"],
            &["(a+b)/2", "a", "(a+b)/2"],
            &["(a+b)/2", "(a+b)/2", "-a"],
            &["2*(a+b)", "-2*(a+b)", "c"],
            // Merged terms that become a sum are flattened.
            &["3*(a+b)", "-2*(a+b)", "c"],
            &["c", "3*(a+b)", "-2*(a+b)"],
            &["2*(a+b)*x", "-(a+b)*x", "y"],
            &["x", "y", "-x", "-y"],
            // A coefficient that becomes 1 leaves the sum, which `+` flattens:
            // later multiples of the sum no longer merge with it.
            &["(a+b)/2", "(a+b)/2", "(a+b)/2"],
            &["(a+b)/3", "2*(a+b)/3", "(a+b)/3"],
            &["(a+b)/2", "(a+b)/2", "-(a+b)/2"],
            &["x", "(a+b)/2", "y", "(a+b)/2", "x", "(a+b)/2", "-y"],
            &["(1/2+1𝑖)*(a+b)", "(1/2-1𝑖)*(a+b)", "2*(a+b)", "c"],
            &["(x+(a+b)/2)/2", "(x+(a+b)/2)/2", "(a+b)/2", "x", "(a+b)/2"],
            // With another factor the merged coefficient never leaves a sum.
            &["x*(a+b)/2", "x*(a+b)/2", "x*(a+b)/2"],
            &["(a+b)^2/2", "(a+b)^2/2", "(a+b)^2/2"],
            // `+` keeps an identical running sum factored.
            &["a+b", "a+b", "c"],
            &["a", "b", "a+b", "c"],
            &["a", "b", "a+b", "-2*a", "-2*b"],
            &["a+b", "c", "a+b+c", "d"],
            &["a+b", "-a-b", "x"],
            &[
                "dot(k(1),k(1))*nf",
                "dot(k(1),k(2))*(1+ep)",
                "-dot(k(1),k(1))*nf",
            ],
            &["1.5*x", "2*x", "y"],
            &["x", "1.5*x", "y", "a+b", "a+b", "3*y"],
            &["1𝑖*x", "2*x", "-1𝑖*x"],
            &["f(a+b)", "f(a+b)", "g(1,2)*x^-1", "x^-1*g(1,2)"],
        ];
        for case in cases {
            assert_same(case);
        }
        // A long run with many merges, sums in between and a factored repeat.
        let mut long = Vec::new();
        for i in 0..200 {
            long.push(format!("{}*x^{}*y^{}", i % 7 - 3, i % 5, i % 3));
            if i % 37 == 0 {
                long.push(format!("x^{i}+y"));
            }
        }
        long.push("x^0+y".to_string());
        long.push("x^0+y".to_string());
        let long = long.iter().map(String::as_str).collect::<Vec<_>>();
        assert_same(&long);
    }

    /// Random sequences of terms that merge with each other and into sums:
    /// sums, multiples of sums (also nested and complex), roots, functions
    /// and cancellations.
    #[test]
    fn exact_sum_reproduces_the_left_fold_on_random_sequences() {
        let pool = "x -x 2*x a b -a -b c 0 1 -1 1/2 x^2 -x^2 x*y a+b -a-b a+b+c c+x -c-x \
            (a+b)/2 -(a+b)/2 2*(a+b) -(a+b) 3*(a+b) x*(a+b) -x*(a+b) (1/3)*(c+x) -(1/3)*(c+x) \
            (2/3)*(c+x) (1/3)*(x+c)*y (2/3)*(x+c)*y (1/2+1𝑖)*(a+b) (1/2-1𝑖)*(a+b) (1/2)*(x+y) \
            x+(a+b)/2 (x+(a+b)/2)/2 -(x+(a+b)/2)/2 (1/2)*(x+y+(a+b)/2) x+y+(a+b)/2 (a+b)/2+c \
            2^(1/2) -2^(1/2) x*2^(1/2) 3^(1/2) 2*3^(1/2) f(a+b) -f(a+b) (a+b)^-1 -(a+b)^-1 \
            2𝑖*x -2𝑖*x (a+b)^2 -(a+b)^2 a^2+2*a*b+b^2";
        let terms = pool
            .split_whitespace()
            .map(|source| vk_parse!(source).unwrap())
            .collect::<Vec<_>>();
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        let mut below = |n: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        };
        for case in 0..20_000 {
            let length = 1 + below(10);
            let sequence = (0..length)
                .map(|_| terms[below(terms.len())].clone())
                .collect::<Vec<_>>();
            let expected = fold(&sequence);
            let actual = merged(&sequence);
            assert!(
                actual == expected
                    && actual.to_canonical_string() == expected.to_canonical_string(),
                "case {case}: {:?}: fold {expected}, merged {actual}",
                sequence.iter().map(Atom::to_string).collect::<Vec<_>>()
            );
        }
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
