use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use linnet::half_edge::involution::EdgeIndex;
use schemars::{JsonSchema, json_schema};
use serde::{Deserialize, Serialize};
use spenso::shadowing::symbolica_utils::SpensoPrintSettings;
use std::{
    collections::{BTreeMap, HashSet},
    fmt::{Display, Write},
    ops::{Deref, DerefMut},
    sync::{
        LazyLock, Mutex,
        atomic::{AtomicUsize, Ordering},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Indeterminate, Symbol},
    coefficient::{Coefficient, CoefficientView},
    domains::{
        algebraic_number::AlgebraicExtension,
        finite_field::{FiniteFieldCore, PrimeIteratorU64, Zp64},
        float::{Complex, FloatLike},
        integer::IntegerRing,
        rational::{FractionField, Q, Rational},
    },
    function,
    id::{ReplaceWith, Replacement},
    parse,
    poly::{polynomial::PolynomialRing, series::SeriesDepth},
    printer::{PrintMode, PrintOptions, PrintState},
    symbol,
};

use crate::GammaLoopContext;

use super::{GS, W_};

const GAMMALOOP_ALIAS_SYMBOL: &str = "gammalooprs::alias";

static GAMMALOOP_ALIAS_COUNTER: AtomicUsize = AtomicUsize::new(0);
static GAMMALOOP_ALIAS_STORE: LazyLock<Mutex<BTreeMap<Atom, Atom>>> =
    LazyLock::new(|| Mutex::new(BTreeMap::new()));

pub(crate) fn alias_subexpressions_into_global_store(
    atom: Atom,
    min_byte_size: Option<usize>,
) -> Atom {
    let aliased = atom.alias_subexpressions(|subexpr, _, _| {
        alias_subexpression_is_allowed(subexpr, min_byte_size)
            .then(|| gammaloop_alias_atom(subexpr))
    });

    if !aliased.get_aliases().is_empty() {
        let mut store = GAMMALOOP_ALIAS_STORE
            .lock()
            .expect("GammaLoop alias store lock poisoned");
        for (alias, body) in aliased.get_aliases() {
            store.insert(alias.clone(), body.clone());
        }
    }

    aliased.get_root().clone()
}

pub(crate) fn get_all_aliases(atom: AtomView<'_>) -> Vec<Atom> {
    let mut aliases = Vec::new();
    let mut seen_aliases = HashSet::new();
    collect_all_aliases(atom, &mut seen_aliases, &mut aliases);
    aliases
}

pub(crate) fn get_alias_expanded_byte_size(atom: AtomView<'_>) -> usize {
    if let Some(alias_body) = lookup_gammaloop_alias(atom) {
        return get_alias_expanded_byte_size(alias_body.as_atom_view());
    }

    match atom {
        AtomView::Num(_) | AtomView::Var(_) => atom.get_byte_size(),
        AtomView::Fun(fun) => expanded_container_byte_size(atom, fun.iter()),
        AtomView::Pow(pow) => expanded_container_byte_size(atom, pow.iter()),
        AtomView::Mul(mul) => expanded_container_byte_size(atom, mul.iter()),
        AtomView::Add(add) => expanded_container_byte_size(atom, add.iter()),
    }
}

pub(crate) fn lookup_gammaloop_alias(atom: AtomView<'_>) -> Option<Atom> {
    is_gammaloop_alias_atom(atom).then(|| {
        GAMMALOOP_ALIAS_STORE
            .lock()
            .expect("GammaLoop alias store lock poisoned")
            .get::<[u8]>(atom.get_data())
            .cloned()
    })?
}

pub(crate) fn inline_gammaloop_aliases(mut atom: Atom) -> Atom {
    loop {
        let expanded = atom.replace_map(|subexpr, _, out| {
            if let Some(alias_body) = lookup_gammaloop_alias(subexpr) {
                **out = alias_body;
            }
        });
        if expanded == atom {
            return atom;
        }
        atom = expanded;
    }
}

fn expanded_container_byte_size<'a>(
    atom: AtomView<'a>,
    children: impl IntoIterator<Item = AtomView<'a>>,
) -> usize {
    children
        .into_iter()
        .fold(atom.get_byte_size(), |byte_size, child| {
            byte_size - child.get_byte_size() + get_alias_expanded_byte_size(child)
        })
}

fn collect_all_aliases(
    atom: AtomView<'_>,
    seen_aliases: &mut HashSet<Atom>,
    aliases: &mut Vec<Atom>,
) {
    if let Some(alias_body) = lookup_gammaloop_alias(atom) {
        let alias_atom = atom.to_owned();
        if seen_aliases.insert(alias_atom) {
            collect_all_aliases(alias_body.as_atom_view(), seen_aliases, aliases);
            aliases.push(alias_body);
        }
        return;
    }

    match atom {
        AtomView::Num(_) | AtomView::Var(_) => {}
        AtomView::Fun(fun) => {
            for child in fun {
                collect_all_aliases(child, seen_aliases, aliases);
            }
        }
        AtomView::Pow(pow) => {
            for child in pow {
                collect_all_aliases(child, seen_aliases, aliases);
            }
        }
        AtomView::Mul(mul) => {
            for child in mul {
                collect_all_aliases(child, seen_aliases, aliases);
            }
        }
        AtomView::Add(add) => {
            for child in add {
                collect_all_aliases(child, seen_aliases, aliases);
            }
        }
    }
}

pub(crate) fn collect_gammaloop_alias_entries(atom: AtomView<'_>) -> Vec<(Atom, Atom)> {
    let mut aliases = Vec::new();
    let mut seen_aliases = HashSet::new();
    collect_gammaloop_alias_entries_impl(atom, &mut seen_aliases, &mut aliases);
    aliases
}

fn collect_gammaloop_alias_entries_impl(
    atom: AtomView<'_>,
    seen_aliases: &mut HashSet<Atom>,
    aliases: &mut Vec<(Atom, Atom)>,
) {
    if let Some(alias_body) = lookup_gammaloop_alias(atom) {
        let alias_atom = atom.to_owned();
        if seen_aliases.insert(alias_atom.clone()) {
            collect_gammaloop_alias_entries_impl(alias_body.as_atom_view(), seen_aliases, aliases);
            aliases.push((alias_atom, alias_body));
        }
        return;
    }

    match atom {
        AtomView::Num(_) | AtomView::Var(_) => {}
        AtomView::Fun(fun) => {
            for child in fun {
                collect_gammaloop_alias_entries_impl(child, seen_aliases, aliases);
            }
        }
        AtomView::Pow(pow) => {
            for child in pow {
                collect_gammaloop_alias_entries_impl(child, seen_aliases, aliases);
            }
        }
        AtomView::Mul(mul) => {
            for child in mul {
                collect_gammaloop_alias_entries_impl(child, seen_aliases, aliases);
            }
        }
        AtomView::Add(add) => {
            for child in add {
                collect_gammaloop_alias_entries_impl(child, seen_aliases, aliases);
            }
        }
    }
}

fn alias_subexpression_is_allowed(subexpr: AtomView<'_>, min_byte_size: Option<usize>) -> bool {
    match subexpr {
        AtomView::Num(_) | AtomView::Var(_) | AtomView::Pow(_) => false,
        AtomView::Fun(fun) if is_structural_index_function(fun.get_symbol()) => false,
        _ if is_gammaloop_alias_atom(subexpr) => false,
        _ => min_byte_size.is_none_or(|min_byte_size| subexpr.get_byte_size() >= min_byte_size),
    }
}

fn is_structural_index_function(symbol: Symbol) -> bool {
    matches!(
        symbol.get_stripped_name(),
        "cind" | "mink" | "aind" | "uind" | "dind"
    )
}

fn gammaloop_alias_atom(subexpr: AtomView<'_>) -> Atom {
    let alias_index = GAMMALOOP_ALIAS_COUNTER.fetch_add(1, Ordering::Relaxed);
    let mut args = subexpr
        .get_all_indeterminates(true)
        .into_iter()
        .map(|arg| arg.to_owned())
        .filter(|arg| {
            !is_gammaloop_alias_atom(arg.as_atom_view())
                && !matches!(
                    arg.as_atom_view(),
                    AtomView::Fun(fun) if is_structural_index_function(fun.get_symbol())
                )
        })
        .collect::<Vec<_>>();
    args.sort();

    let mut alias = FunctionBuilder::new(symbol!(GAMMALOOP_ALIAS_SYMBOL));
    alias = alias.add_arg(alias_index);
    for arg in args {
        alias = alias.add_arg(arg);
    }
    alias.finish()
}

pub(crate) fn is_gammaloop_alias_atom(atom: AtomView<'_>) -> bool {
    let AtomView::Fun(fun) = atom else {
        return false;
    };
    fun.get_symbol().get_name() == GAMMALOOP_ALIAS_SYMBOL && fun.get_nargs() >= 1
}

pub(crate) fn add_numeric_constant_to_fn_map(
    fn_map: &mut symbolica::evaluate::FunctionMap,
    name: Atom,
    value: impl Into<Coefficient>,
) -> std::result::Result<(), String> {
    let value = value.into();
    if name.is_zero() {
        return if value.is_zero() {
            Ok(())
        } else {
            Err(format!("Cannot bind the literal zero atom to {value}"))
        };
    }

    let indeterminate = Indeterminate::try_from(name.clone())?;
    fn_map.add_function(indeterminate, Vec::<Symbol>::new(), Atom::num(value))
}

pub static Q_I: LazyLock<AlgebraicExtension<FractionField<IntegerRing>>> =
    LazyLock::new(|| AlgebraicExtension::new_complex(Q));
static RAW_UFO_MOMENTUM: LazyLock<Symbol> = LazyLock::new(|| symbol!("UFO::P"));
static RAW_UFO_PSLASH: LazyLock<Symbol> = LazyLock::new(|| symbol!("UFO::PSlash"));

pub static COMPLEXRATPOLYFIELD: LazyLock<
    FractionField<PolynomialRing<AlgebraicExtension<FractionField<IntegerRing>>, u16>>,
> = LazyLock::new(|| FractionField::new(PolynomialRing::<_, u16>::new(Q_I.clone())));

pub static LOGPRINTOPTS: PrintOptions = PrintOptions {
    hide_all_namespaces: true,
    color_namespace: false,
    color_builtin_symbols: false,
    color_top_level_sum: false,
    terms_on_new_line: false,
    print_ring: false,
    include_attributes: false,
    symmetric_representation_for_finite_field: false,
    explicit_rational_polynomial: false,
    number_thousands_separator: None,
    multiplication_operator: '*',
    double_star_for_exponentiation: false,
    num_exp_as_superscript: false,
    mode: PrintMode::Symbolica,
    precision: None,
    pretty_matrix: false,
    max_terms: None,
    custom_print_mode: None,
    hide_namespace: Some("gammalooprs"),
    ..PrintOptions::new()
};

pub trait IsNeg {
    fn is_negative(&self) -> bool;
}

impl<A: AtomCore> IsNeg for A {
    // fn is_one(&self)
    fn is_negative(&self) -> bool {
        match self.as_atom_view() {
            AtomView::Num(a) => match a.get_coeff_view() {
                CoefficientView::FiniteField(_, _) => false,
                CoefficientView::Float(re, im) => {
                    if im.is_zero() {
                        re.to_float().is_negative()
                    } else {
                        false
                    }
                }
                CoefficientView::Indeterminate => false,
                CoefficientView::Natural(n_re, _d_re, n_im, _de_im) => {
                    if n_im == 0 {
                        n_re.is_negative()
                    } else {
                        false
                    }
                }
                CoefficientView::Infinity(_) => false,
                CoefficientView::Large(re, im) => {
                    let re = re.to_rat();
                    let im = im.to_rat();
                    if im.is_zero() {
                        re.is_negative()
                    } else {
                        false
                    }
                }
                CoefficientView::RationalPolynomial(_) => false,
            },
            AtomView::Mul(a) => {
                for term in a.iter() {
                    if term.is_negative() {
                        return true;
                    }
                }
                false
                // if let Some(first) = a.iter().next() {
                //     first.is_negative()
                // } else {
                //     false
                // }
            }
            _ => false,
        }
    }
}

#[test]
fn isneg() {
    let a = Atom::num(-3) * Atom::var(GS.emr_mom);

    assert!(a.is_negative());
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct TypstState {
    pub inproduct: bool,
    pub without_minus: bool,
    pub in_power: bool,
    pub in_function: bool,
}
pub trait TypstFormat {
    fn typst_string(&self) -> String {
        let mut s = String::new();

        let symbols = self.preable(&mut s).unwrap();
        s.push_str("$ ");
        self.fmt_output(
            &mut s,
            &symbols,
            TypstState {
                inproduct: false,
                without_minus: false,
                in_power: false,
                in_function: false,
            },
        )
        .unwrap();
        s.push_str(" $");
        s
    }

    fn fmt_output<W: std::fmt::Write>(
        &self,
        f: &mut W,
        symbols: &BTreeMap<Symbol, bool>,
        print_state: TypstState,
    ) -> Result<bool>;

    fn preable<W: std::fmt::Write>(&self, f: &mut W) -> Result<BTreeMap<Symbol, bool>>;
}

impl<A: AtomCore> TypstFormat for A {
    fn preable<W: std::fmt::Write>(&self, fmt: &mut W) -> Result<BTreeMap<Symbol, bool>> {
        let mut symbols = BTreeMap::new();
        self.visitor(&mut |a| match a {
            AtomView::Add(_) => true,
            AtomView::Mul(_) => true,
            AtomView::Fun(f) => {
                *(symbols.entry(f.get_symbol()).or_insert_with(|| true)) = true;
                true
            }
            AtomView::Num(_) => false,
            AtomView::Pow(_) => true,
            AtomView::Var(v) => {
                symbols.entry(v.get_symbol()).or_insert_with(|| false);
                false
            }
        });

        writeln!(
            fmt,
            r#"#let defaultfn(..args, name:"NA") = {{
let a = args.pos().map(v => $#v$).join($, $);
    if args.pos().len() != 0 {{
        $op(name)(#a)$
    }} else {{
        $#name$
    }}
}}
           "#
        )?;

        for (s, isfun) in &symbols {
            if let Some(p) = s.get_print_function()
                && let Some(s) = p(
                    if *isfun {
                        Atom::var(GS.is_function)
                    } else {
                        Atom::var(GS.is_symbol)
                    }
                    .as_view(),
                    &PrintOptions {
                        custom_print_mode: Some(("typst", 2)),
                        ..Default::default()
                    },
                    &PrintState::new(),
                )
            {
                writeln!(fmt, "{}", s).unwrap();
                continue;
            }

            if *isfun {
                writeln!(
                    fmt,
                    "#let {}-{}(..args)= defaultfn(..args,name:\"{}\")",
                    s.get_namespace(),
                    s.get_stripped_name(),
                    s.get_stripped_name()
                )
                .unwrap();
            } else {
                writeln!(
                    fmt,
                    "#let {}-{}= $\"{}\"$",
                    s.get_namespace(),
                    s.get_stripped_name(),
                    s.get_stripped_name()
                )
                .unwrap();
            }
        }

        Ok(symbols)
    }
    fn fmt_output<W: std::fmt::Write>(
        &self,
        fmt: &mut W,
        symbols: &BTreeMap<Symbol, bool>,
        mut print_state: TypstState,
    ) -> Result<bool> {
        match self.as_atom_view() {
            AtomView::Num(n) => match n.get_coeff_view() {
                CoefficientView::FiniteField(_, _) => Ok(true),
                CoefficientView::Float(re, im) => {
                    let re = re.to_float();
                    let im = im.to_float();

                    if im.is_fully_zero() {
                        if print_state.without_minus && re.is_negative() {
                            write!(fmt, "{}", -re)?;
                        } else {
                            write!(fmt, "{}", re)?;
                        }
                    } else {
                        if print_state.inproduct || print_state.in_power {
                            write!(fmt, "(")?;
                        }
                        write!(fmt, "{} + {} i", re, im)?;
                        if print_state.inproduct || print_state.in_power {
                            write!(fmt, ")")?;
                        }
                    }
                    Ok(true)
                }
                CoefficientView::Indeterminate => {
                    write!(fmt, "NA")?;
                    Ok(true)
                }
                CoefficientView::Natural(n_re, d_re, n_im, de_im) => {
                    if d_re == 1 && de_im == 1 {
                        if n_im == 0 {
                            if print_state.without_minus && n_re.is_negative() {
                                write!(fmt, "{}", -n_re)?;
                            } else {
                                write!(fmt, "{}", n_re)?;
                            }
                        } else if n_re == 0 {
                            if n_im == 1 {
                                write!(fmt, "i")?;
                            } else {
                                write!(fmt, "{} i", n_im)?;
                            }
                        } else {
                            if print_state.inproduct || print_state.in_power {
                                write!(fmt, "(")?;
                            }
                            write!(fmt, "{} + {} i", n_re, n_im)?;
                            if print_state.inproduct || print_state.in_power {
                                write!(fmt, ")")?;
                            }
                        }
                    } else if n_im == 0 {
                        if d_re == 1 {
                            if print_state.without_minus && n_re.is_negative() {
                                write!(fmt, "{}", -n_re)?;
                            } else {
                                write!(fmt, "{}", n_re)?;
                            }
                            return Ok(true);
                        }
                        if print_state.without_minus && n_re.is_negative() {
                            write!(fmt, "({})/({})", -n_re, d_re)?;
                        } else {
                            write!(fmt, "({})/({})", n_re, d_re)?;
                        }
                        return Ok(true);
                    } else {
                        if de_im == 1 {
                            if print_state.inproduct || print_state.in_power {
                                write!(fmt, "(")?;
                            }
                            write!(fmt, "({})/({}) + {} i", n_re, d_re, n_im)?;
                            if print_state.inproduct || print_state.in_power {
                                write!(fmt, ")")?;
                            }
                            return Ok(true);
                        }
                        if print_state.inproduct || print_state.in_power {
                            write!(fmt, "(")?;
                        }
                        write!(fmt, "({})/({})+ ({})/({})i", n_re, d_re, n_im, de_im)?;
                        if print_state.inproduct || print_state.in_power {
                            write!(fmt, ")")?;
                        }
                    }
                    Ok(true)
                }
                CoefficientView::Infinity(_) => {
                    write!(fmt, "oo")?;
                    Ok(true)
                }
                CoefficientView::Large(re, im) => {
                    let re = re.to_rat();
                    let im = im.to_rat();
                    if im.is_zero() {
                        if re.is_integer() {
                            write!(fmt, "{}", re.numerator())?;
                        } else {
                            write!(fmt, "({})/({})", re.numerator(), re.denominator())?;
                        }
                    } else if re.is_zero() {
                        if im.is_integer() {
                            if im.is_one() {
                                write!(fmt, "i")?;
                            } else {
                                write!(fmt, "{} i", im.numerator())?;
                            }
                        } else {
                            write!(fmt, "({})/({}) i", im.numerator(), im.denominator())?;
                        }
                    } else {
                        if print_state.inproduct || print_state.in_power {
                            write!(fmt, "(")?;
                        }
                        if re.is_integer() && im.is_integer() {
                            write!(fmt, "{} + {} i", re.numerator(), im.numerator())?;
                        } else if re.is_integer() {
                            write!(
                                fmt,
                                "{} + ({})/({}) i",
                                re.numerator(),
                                im.numerator(),
                                im.denominator()
                            )?;
                        } else if im.is_integer() {
                            write!(
                                fmt,
                                "({})/({}) + {} i",
                                re.numerator(),
                                re.denominator(),
                                im.numerator()
                            )?;
                        } else {
                            write!(
                                fmt,
                                "({})/({}) + ({})/({}) i",
                                re.numerator(),
                                re.denominator(),
                                im.numerator(),
                                im.denominator()
                            )?;
                        }
                        if print_state.inproduct || print_state.in_power {
                            write!(fmt, ")")?;
                        }
                    }
                    Ok(true)
                }
                CoefficientView::RationalPolynomial(a) => {
                    write!(fmt, "{}", a.deserialize())?;
                    Ok(true)
                }
            },

            AtomView::Var(v) => {
                if !print_state.in_function {
                    write!(fmt, "#")?;
                }
                write!(
                    fmt,
                    "{}-{}",
                    v.get_symbol().get_namespace(),
                    v.get_symbol().get_stripped_name()
                )?;
                Ok(true)
            }
            AtomView::Fun(f) => {
                if !print_state.in_function {
                    write!(fmt, "#")?;
                }
                write!(
                    fmt,
                    "{}-{}(",
                    f.get_symbol().get_namespace(),
                    f.get_symbol().get_stripped_name()
                )?;

                let nargs = f.get_nargs();

                for (j, i) in f.iter().enumerate() {
                    i.fmt_output(
                        fmt,
                        symbols,
                        TypstState {
                            inproduct: false,
                            without_minus: false,
                            in_power: false,
                            in_function: true,
                        },
                    )?;
                    if j + 1 != nargs {
                        write!(fmt, ", ")?;
                    }
                }
                write!(fmt, ")")?;

                Ok(true)
            }
            AtomView::Pow(p) => {
                let (base, exp) = p.get_base_exp();
                if print_state.in_function {
                    write!(fmt, "$")?;
                }
                if exp.is_negative() {
                    write!(fmt, "1/")?;
                    base.fmt_output(
                        fmt,
                        symbols,
                        TypstState {
                            inproduct: false,
                            without_minus: false,
                            in_power: true,
                            in_function: false,
                        },
                    )?;

                    if (-exp).is_one() {
                        return Ok(true);
                    }
                    write!(fmt, "^(")?;
                    exp.fmt_output(
                        fmt,
                        symbols,
                        TypstState {
                            inproduct: false,
                            without_minus: true,
                            in_power: false,
                            in_function: false,
                        },
                    )?;
                    write!(fmt, ")")?;
                } else {
                    base.fmt_output(
                        fmt,
                        symbols,
                        TypstState {
                            inproduct: false,
                            without_minus: false,
                            in_power: true,
                            in_function: false,
                        },
                    )?;
                    if (exp).is_one() {
                        return Ok(true);
                    }
                    write!(fmt, "^(")?;
                    exp.fmt_output(
                        fmt,
                        symbols,
                        TypstState {
                            inproduct: false,
                            without_minus: false,
                            in_power: false,
                            in_function: false,
                        },
                    )?;

                    write!(fmt, ")")?;
                }
                if print_state.in_function {
                    write!(fmt, "$")?;
                }
                Ok(true)
            }

            AtomView::Mul(t) => {
                print_state.inproduct = true;
                let infun = print_state.in_function;
                print_state.in_function = false;
                if infun {
                    write!(fmt, "$")?;
                }
                let mut numbuffer = String::new();
                let mut denombuffer = String::new();

                for term in t.iter() {
                    if let AtomView::Pow(p) = term.as_atom_view() {
                        let (base, exp) = p.get_base_exp();
                        if exp.is_negative() {
                            let state = TypstState {
                                inproduct: false,
                                without_minus: false,
                                in_power: true,
                                in_function: false,
                            };
                            // if (-exp).is_one() {
                            //     state.in_function = true;
                            // }
                            base.fmt_output(&mut denombuffer, symbols, state)?;

                            if (-exp).is_one() {
                                continue;
                            }
                            write!(&mut denombuffer, "^(")?;
                            exp.fmt_output(
                                &mut denombuffer,
                                symbols,
                                TypstState {
                                    inproduct: false,
                                    without_minus: true,
                                    in_power: false,
                                    in_function: false,
                                },
                            )?;

                            write!(&mut denombuffer, ")")?;
                            write!(&mut denombuffer, " ")?;

                            continue;
                        } else {
                            term.fmt_output(&mut numbuffer, symbols, print_state)?;
                        }
                    } else if let Ok(t) = i32::try_from(term)
                        && t == -1
                    {
                        if !print_state.without_minus {
                            write!(fmt, "-")?;
                        } else {
                            continue;
                        }
                    } else {
                        term.fmt_output(&mut numbuffer, symbols, print_state)?;
                    }
                    write!(&mut numbuffer, " ")?;
                }

                if !denombuffer.is_empty() {
                    write!(
                        fmt,
                        "({})/({})",
                        numbuffer.trim_end(),
                        denombuffer.trim_end()
                    )?;
                } else {
                    write!(fmt, "{}", numbuffer.trim_end())?;
                }

                if infun {
                    write!(fmt, "$")?;
                }

                Ok(true)
            }

            AtomView::Add(e) => {
                let infun = print_state.in_function;
                print_state.in_function = false;
                if infun {
                    write!(fmt, "$")?;
                }
                if print_state.inproduct || print_state.in_power {
                    write!(fmt, "(")?;
                }
                let mut first = true;
                for term in e.iter() {
                    if term.is_negative() {
                        write!(fmt, " - ")?;
                    } else if !first {
                        write!(fmt, " + ")?;
                    }
                    first = false;
                    print_state.without_minus = true;
                    term.fmt_output(fmt, symbols, print_state)?;
                }

                if print_state.inproduct || print_state.in_power {
                    write!(fmt, ")")?;
                }
                if infun {
                    write!(fmt, "$")?;
                }
                Ok(true)
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Default, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct StringSerializedAtom(pub Atom);

impl JsonSchema for StringSerializedAtom {
    fn schema_name() -> std::borrow::Cow<'static, str> {
        "ParseableAtom".into()
    }
    fn json_schema(_generator: &mut schemars::SchemaGenerator) -> schemars::Schema {
        json_schema!({
            "description": "An atom that is serialized as a string. Do not make any assumptions about the state",
            "type": ["string"]
        })
    }
}

impl Deref for StringSerializedAtom {
    type Target = Atom;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl DerefMut for StringSerializedAtom {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
impl Display for StringSerializedAtom {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Serialize for StringSerializedAtom {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.0.to_canonical_string().serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for StringSerializedAtom {
    fn deserialize<D>(deserializer: D) -> std::result::Result<StringSerializedAtom, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        Ok(StringSerializedAtom(parse!(String::deserialize(
            deserializer
        )?)))
    }
}

pub trait PrimeGenerate {
    fn prime_generate_int(sample_iterator: &mut PrimeIteratorU64) -> Atom;
    fn prime_generate_ff_complex(
        sample_iterator: &mut PrimeIteratorU64,
        finite_field: &Zp64,
    ) -> Atom;
    fn prime_generate_ff(sample_iterator: &mut PrimeIteratorU64, finite_field: &Zp64) -> Atom;
    fn prime_generate_rat_complex(sample_iterator: &mut PrimeIteratorU64) -> Atom;
    fn prime_generate_rat(sample_iterator: &mut PrimeIteratorU64) -> Atom;
}

impl PrimeGenerate for Atom {
    fn prime_generate_int(sample_iterator: &mut PrimeIteratorU64) -> Atom {
        Atom::num(sample_iterator.next().unwrap())
    }

    fn prime_generate_rat(sample_iterator: &mut PrimeIteratorU64) -> Atom {
        let a = Rational::from((
            sample_iterator.next().unwrap(),
            sample_iterator.next().unwrap(),
        ));
        Atom::num(a)
    }

    fn prime_generate_ff_complex(
        sample_iterator: &mut PrimeIteratorU64,
        finite_field: &Zp64,
    ) -> Atom {
        let _a = finite_field.to_element(sample_iterator.next().unwrap());
        let _b = finite_field.to_element(sample_iterator.next().unwrap());
        Atom::num(1)
        // Atom::i() * Atom::num(b) + Atom::num(a)
    }

    fn prime_generate_ff(sample_iterator: &mut PrimeIteratorU64, _finite_field: &Zp64) -> Atom {
        let a = Rational::from((
            sample_iterator.next().unwrap(),
            sample_iterator.next().unwrap(),
        ));
        let b = Rational::from((
            sample_iterator.next().unwrap(),
            sample_iterator.next().unwrap(),
        ));
        Atom::num(Complex::new(a, b))
    }

    fn prime_generate_rat_complex(sample_iterator: &mut PrimeIteratorU64) -> Atom {
        let a = Rational::from((
            sample_iterator.next().unwrap(),
            sample_iterator.next().unwrap(),
        ));
        let b = Rational::from((
            sample_iterator.next().unwrap(),
            sample_iterator.next().unwrap(),
        ));
        Atom::num(Complex::new(a, b))
    }
}

pub trait DOD {
    ///Rescales momentum of edge `eid`, and computes the leading scaling
    fn edge_dod(&self, eid: EdgeIndex) -> i32;

    ///Rescales all momenta, and computes the leading scaling
    fn all_dod(&self) -> i32;

    fn trailing_exponent(&self) -> i32;
}

impl DOD for Atom {
    fn edge_dod(&self, eid: EdgeIndex) -> i32 {
        self.as_view().edge_dod(eid)
    }

    fn all_dod(&self) -> i32 {
        self.as_view().all_dod()
    }

    fn trailing_exponent(&self) -> i32 {
        self.as_view().trailing_exponent()
    }
}

impl DOD for AtomView<'_> {
    fn edge_dod(&self, eid: EdgeIndex) -> i32 {
        self.replace(GS.emr_mom(eid, W_.a___))
            .with(GS.emr_mom(eid, W_.a___) / GS.rescale)
            .trailing_exponent()
    }

    fn all_dod(&self) -> i32 {
        self.replace(function!(GS.emr_mom, W_.a___))
            .with(function!(GS.emr_mom, W_.a___) / GS.rescale)
            .replace(function!(*RAW_UFO_MOMENTUM, W_.a___))
            .with(function!(*RAW_UFO_MOMENTUM, W_.a___) / GS.rescale)
            .replace(function!(*RAW_UFO_PSLASH, W_.a___))
            .with(function!(*RAW_UFO_PSLASH, W_.a___) / GS.rescale)
            .trailing_exponent()
    }

    fn trailing_exponent(&self) -> i32 {
        let series = self
            .series(GS.rescale, Atom::Zero, SeriesDepth::relative(1))
            .unwrap();
        // println!("{}", series);

        let dod = series.get_trailing_exponent();

        if dod.is_integer() {
            -(dod.numerator().to_i64().unwrap() as i32)
        } else {
            panic!("{dod} for {self}")
        }
    }
}

#[test]
fn test_dod() {
    let (e1, e2) = (EdgeIndex(1), EdgeIndex(2));

    let atom = (GS.emr_mom(e1, Atom::Zero) * GS.emr_mom(e2, Atom::Zero)
        + GS.emr_mom(e2, Atom::Zero))
        / (GS.emr_mom(e1, Atom::Zero) * GS.emr_mom(e1, Atom::Zero));

    let atom2 =
        Atom::num(1) / (GS.emr_mom(e1, Atom::Zero) * GS.emr_mom(e1, Atom::Zero) + parse!("m"));
    let atom3 =
        GS.emr_mom(e1, Atom::Zero) * GS.emr_mom(e2, Atom::Zero) + GS.emr_mom(e2, Atom::Zero);

    assert_eq!(-1, atom.edge_dod(e1));
    assert_eq!(1, atom.edge_dod(e2));
    assert_eq!(-2, atom2.edge_dod(e1));
    assert_eq!(1, atom3.edge_dod(e1));
    assert_eq!(1, atom3.edge_dod(e2));
}

pub trait CallSymbol<T> {
    fn f(&self, args: T) -> Atom;
}

impl<'a> CallSymbol<AtomOrView<'a>> for Symbol {
    fn f(&self, arg: AtomOrView<'a>) -> Atom {
        FunctionBuilder::new(*self).add_arg(arg).finish()
    }
}

impl<'a> CallSymbol<AtomView<'a>> for Symbol {
    fn f(&self, arg: AtomView<'a>) -> Atom {
        FunctionBuilder::new(*self).add_arg(arg).finish()
    }
}

impl CallSymbol<Atom> for Symbol {
    fn f(&self, arg: Atom) -> Atom {
        FunctionBuilder::new(*self).add_arg(arg).finish()
    }
}

impl CallSymbol<&Atom> for Symbol {
    fn f(&self, arg: &Atom) -> Atom {
        FunctionBuilder::new(*self).add_arg(arg).finish()
    }
}

impl CallSymbol<usize> for Symbol {
    fn f(&self, arg: usize) -> Atom {
        FunctionBuilder::new(*self).add_arg(arg).finish()
    }
}

impl CallSymbol<Symbol> for Symbol {
    fn f(&self, arg: Symbol) -> Atom {
        FunctionBuilder::new(*self).add_arg(arg).finish()
    }
}

impl<'a, I> CallSymbol<&'a [I]> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
{
    fn f(&self, args: &'a [I]) -> Atom {
        FunctionBuilder::new(*self).add_args(args).finish()
    }
}
impl<'a, I, const N: usize> CallSymbol<&'a [I; N]> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
{
    fn f(&self, args: &'a [I; N]) -> Atom {
        FunctionBuilder::new(*self).add_args(args).finish()
    }
}

impl<'a, 'b, I, J> CallSymbol<(&'a [I], &'b [J])> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
    &'b J: Into<AtomOrView<'b>>,
{
    fn f(&self, args: (&'a [I], &'b [J])) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.0)
            .add_args(args.1)
            .finish()
    }
}

impl<'a, 'b, I, J, const N: usize> CallSymbol<(&'a [I; N], &'b [J])> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
    &'b J: Into<AtomOrView<'b>>,
{
    fn f(&self, args: (&'a [I; N], &'b [J])) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.0.as_slice())
            .add_args(args.1)
            .finish()
    }
}

impl<'a, 'b, I, J, const N: usize, const M: usize> CallSymbol<(&'a [I; N], &'b [J; M])> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
    &'b J: Into<AtomOrView<'b>>,
{
    fn f(&self, args: (&'a [I; N], &'b [J; M])) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.0.as_slice())
            .add_args(args.1)
            .finish()
    }
}
impl<'a, 'b, 'c, I, J, K, const N: usize, const M: usize, const O: usize>
    CallSymbol<(&'a [I; N], &'b [J; M], &'c [K; O])> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
    &'b J: Into<AtomOrView<'b>>,
    &'c K: Into<AtomOrView<'c>>,
{
    fn f(&self, args: (&'a [I; N], &'b [J; M], &'c [K; O])) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.0.as_slice())
            .add_args(args.1)
            .add_args(args.2)
            .finish()
    }
}

impl<I, const N: usize> CallSymbol<[I; N]> for Symbol
where
    for<'a> &'a I: Into<AtomOrView<'a>>,
{
    fn f(&self, args: [I; N]) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.as_slice())
            .finish()
    }
}

impl<'b, I, J, const N: usize> CallSymbol<([I; N], &'b [J])> for Symbol
where
    for<'a> &'a I: Into<AtomOrView<'a>>,
    &'b J: Into<AtomOrView<'b>>,
{
    fn f(&self, args: ([I; N], &'b [J])) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.0.as_slice())
            .add_args(args.1)
            .finish()
    }
}

impl<'a, I, J, K> CallSymbol<(&'a [I], &'a [J], &'a [K])> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
    &'a J: Into<AtomOrView<'a>>,
    &'a K: Into<AtomOrView<'a>>,
{
    fn f(&self, args: (&'a [I], &'a [J], &'a [K])) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.0)
            .add_args(args.1)
            .add_args(args.2)
            .finish()
    }
}

pub trait LogPrint {
    fn log_print(&self, max_line_length: Option<usize>) -> String;
}

impl<A: AtomCore> LogPrint for A {
    fn log_print(&self, max: Option<usize>) -> String {
        let mut settings = SpensoPrintSettings::compact().nice_symbolica();
        settings.max_line_length = max;
        // settings.hide_all_namespaces = false;
        self.printer(settings).to_string()
    }
}

pub trait Replaces {
    fn replace_with<R: Into<ReplaceWith<'static>>>(&self, rhs: R) -> Replacement;
}
impl<A: AtomCore> Replaces for A {
    fn replace_with<R: Into<ReplaceWith<'static>>>(&self, rhs: R) -> Replacement {
        Replacement::new(self.to_pattern(), rhs)
    }
}

#[cfg(test)]
mod alias_accounting_tests {
    use symbolica::{atom::AtomCore, function, parse, symbol};

    use super::{
        alias_subexpressions_into_global_store, get_alias_expanded_byte_size, get_all_aliases,
        inline_gammaloop_aliases,
    };

    fn alias_storage_byte_size(aliases: &[symbolica::atom::Atom]) -> usize {
        aliases
            .iter()
            .map(|alias| alias.as_atom_view().get_byte_size())
            .sum()
    }

    #[test]
    fn alias_accounting_has_no_overhead_for_unaliased_atoms() {
        let atom = parse!("probe::f(x+y,z)^2 + probe::g(x)");

        assert!(get_all_aliases(atom.as_atom_view()).is_empty());
        assert_eq!(
            get_alias_expanded_byte_size(atom.as_atom_view()),
            atom.as_atom_view().get_byte_size()
        );
    }

    #[test]
    fn alias_accounting_counts_repeated_alias_expansion_and_storage_once() {
        let body = parse!("x+y+z");
        let atom = function!(symbol!("probe::f"), body.clone())
            + function!(symbol!("probe::g"), body.clone())
            + function!(symbol!("probe::h"), body.clone());
        let aliased = alias_subexpressions_into_global_store(atom, None);

        let aliases = get_all_aliases(aliased.as_atom_view());
        let root_byte_size = aliased.as_atom_view().get_byte_size();
        let alias_storage_byte_size = alias_storage_byte_size(&aliases);
        let expression_byte_size = root_byte_size + alias_storage_byte_size;
        let expanded_byte_size = get_alias_expanded_byte_size(aliased.as_atom_view());

        assert!(!aliases.is_empty());
        assert!(alias_storage_byte_size >= body.as_atom_view().get_byte_size());
        let materialized_expanded_byte_size = inline_gammaloop_aliases(aliased.clone())
            .as_atom_view()
            .get_byte_size();
        assert_eq!(expanded_byte_size, materialized_expanded_byte_size);
        assert!(expression_byte_size > 0);
    }

    #[test]
    fn alias_accounting_collects_nested_alias_storage_once() {
        let inner_body = parse!("x+y+z");
        let inner_alias = alias_subexpressions_into_global_store(
            function!(
                symbol!("probe::inner_holder"),
                inner_body.clone(),
                inner_body.clone()
            ),
            None,
        );
        let outer_body = function!(
            symbol!("probe::outer"),
            inner_alias.clone(),
            inner_alias.clone(),
            inner_alias
        );
        let atom = function!(
            symbol!("probe::top"),
            outer_body.clone(),
            outer_body.clone()
        );
        let aliased = alias_subexpressions_into_global_store(atom, None);

        let aliases = get_all_aliases(aliased.as_atom_view());
        let root_byte_size = aliased.as_atom_view().get_byte_size();
        let alias_storage_byte_size = alias_storage_byte_size(&aliases);
        let expression_byte_size = root_byte_size + alias_storage_byte_size;
        let expanded_byte_size = get_alias_expanded_byte_size(aliased.as_atom_view());

        assert!(aliases.len() >= 2);
        let materialized_expanded_byte_size = inline_gammaloop_aliases(aliased.clone())
            .as_atom_view()
            .get_byte_size();
        assert_eq!(expanded_byte_size, materialized_expanded_byte_size);
        assert!(expression_byte_size > root_byte_size);
    }

    #[test]
    fn alias_accounting_matches_faithful_compression_metric() {
        let body =
            parse!("probe::b0(x+y+z+w)+probe::b1(x+y+z+w)+probe::b2(x+y+z+w)+probe::b3(x+y+z+w)");
        let atom = function!(symbol!("probe::f0"), body.clone())
            + function!(symbol!("probe::f1"), body.clone())
            + function!(symbol!("probe::f2"), body.clone())
            + function!(symbol!("probe::f3"), body.clone())
            + function!(symbol!("probe::f4"), body.clone())
            + function!(symbol!("probe::f5"), body.clone())
            + function!(symbol!("probe::f6"), body.clone())
            + function!(symbol!("probe::f7"), body.clone());
        let aliased = alias_subexpressions_into_global_store(atom, None);

        let aliases = get_all_aliases(aliased.as_atom_view());
        let expression_byte_size =
            aliased.as_atom_view().get_byte_size() + alias_storage_byte_size(&aliases);
        let expanded_byte_size = get_alias_expanded_byte_size(aliased.as_atom_view());

        assert!(!aliases.is_empty());
        let materialized_expanded_byte_size = inline_gammaloop_aliases(aliased.clone())
            .as_atom_view()
            .get_byte_size();
        assert!(
            expanded_byte_size.abs_diff(materialized_expanded_byte_size) <= 1,
            "structural expanded byte size {expanded_byte_size} differs from materialized \
             expanded byte size {materialized_expanded_byte_size}"
        );
        assert!(expression_byte_size > 0);
    }
}
