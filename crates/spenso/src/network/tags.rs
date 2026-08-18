use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    coefficient::CoefficientView,
    printer::{PrintOptions, PrintState},
    symbol, tag,
};

use crate::{
    shadowing::symbolica_utils::{SpensoPrintBackend, SpensoPrintSettings},
    structure::abstract_index::AIND_SYMBOLS,
};

pub struct SpensoTags {
    pub broadcast: String,
    /// Marks rank-one tensor symbols whose final argument is the tensor slot.
    ///
    /// Symbols carrying this tag must not use representation slots in earlier
    /// arguments; rank-one shorthand rewrites assume that contract.
    pub rank1: String,
    pub rank1_: Symbol,
    /// Reserved chain wiring markers registered in Spenso's Symbolica namespace.
    pub chain_in: Symbol,
    pub chain_out: Symbol,
    pub chain: Symbol,
    pub trace: Symbol,
    pub upper: String,
    pub lower: String,
    pub bracket: Symbol,
    pub pure_scalar: Symbol,
    pub scalar: Symbol,
    pub tensor: String,
    pub tensor_: Symbol,
    pub index: String,
    pub representation: String,
    pub i_: Symbol,
    pub dot: Symbol,
    pub rep_: Symbol,
    pub self_dual: String,
    pub self_dual_: Symbol,
    pub dualizable: String,
    pub dualizable_: Symbol,
}

crate::symbolica_init_lazy_static! {
    pub static SPENSO_TAG, SPENSO_TAG_INNER: SpensoTags = SpensoTags::new;
}

pub fn scalar_store_alias(index: usize) -> Atom {
    FunctionBuilder::new(SPENSO_TAG.scalar)
        .add_arg(Atom::num(
            i64::try_from(index).expect("scalar alias index must fit in i64"),
        ))
        .finish()
}

pub fn scalar_store_alias_index(value: AtomView<'_>) -> Option<usize> {
    let AtomView::Fun(fun) = value else {
        return None;
    };
    if fun.get_symbol() != SPENSO_TAG.scalar || fun.get_nargs() != 1 {
        return None;
    }

    let AtomView::Num(index) = fun.iter().next()? else {
        return None;
    };
    match index.get_coeff_view() {
        CoefficientView::Natural(index, 1, 0, 1) => usize::try_from(index).ok(),
        _ => None,
    }
}

/// Builds Symbolica atoms from a symbol and an optional argument list.
///
/// Spenso wildcard heads use the convention that a bare head with no arguments
/// is a variable atom, while a head with arguments is a function atom. This
/// helper keeps that convention in one place.
pub trait SymbolAtomExt {
    fn atom_with_args<'a, A>(self, args: impl IntoIterator<Item = A>) -> Atom
    where
        A: Into<AtomOrView<'a>>;
}

impl SymbolAtomExt for Symbol {
    fn atom_with_args<'a, A>(self, args: impl IntoIterator<Item = A>) -> Atom
    where
        A: Into<AtomOrView<'a>>,
    {
        let mut function = FunctionBuilder::new(self);
        let mut has_args = false;
        for arg in args {
            has_args = true;
            function = function.add_arg(arg);
        }

        if has_args {
            function.finish()
        } else {
            Atom::var(self)
        }
    }
}

macro_rules! define_numbered_tag_family_methods {
    ($(
        $(#[$meta:meta])*
        $vis:vis fn $method:ident => $base_field:ident, $prefix:literal, $symbol_method:ident;
    )*) => {
        $(
            $(#[$meta])*
            $vis fn $method<'a, const N: usize, A: Into<AtomOrView<'a>>>(
                &self,
                args: impl IntoIterator<Item = A>,
            ) -> Atom {
                let symbol = if N == 0 {
                    self.$base_field
                } else {
                    self.$symbol_method(&format!("{}{}_", $prefix, N))
                };
                symbol.atom_with_args(args)
            }
        )*
    };
}

macro_rules! define_numbered_tag_macros {
    ($d:tt; $($macro_name:ident => $method:ident;)*) => {
        $(
            #[macro_export]
            macro_rules! $macro_name {
                ($d n:literal; $d($d arg:expr),* $d(,)?) => {
                    $crate::network::tags::SPENSO_TAG.$method::<$d n, _>(
                        vec![$d($crate::shadowing::IntoAtom::into_atom($d arg)),*],
                    )
                };
                ($d n:literal $d(;)?) => {
                    $crate::network::tags::SPENSO_TAG.$method::<$d n, symbolica::atom::Atom>(
                        std::iter::empty::<symbolica::atom::Atom>(),
                    )
                };
            }
        )*
    };
}

define_numbered_tag_macros!($;
    rank1_ => rank1_;
    tensor_ => tensor_;
    rep_ => rep_;
    self_dual_ => self_dual_;
    dualizable_ => dualizable_;
    dualizable_dual_ => dualizable_dual_;
);

/// Creates a tensor-head symbol tagged with Spenso's generic tensor tag.
///
/// This expands `symbolica::symbol!` at the call site, so the symbol keeps the
/// caller's crate namespace while automatically receiving the Spenso tag. Any
/// Symbolica attributes and settings such as `print = ...` are forwarded, while
/// `tag`/`tags` remain owned by this macro so the tensor tag cannot be skipped.
#[macro_export]
macro_rules! tensor_symbol {
    ($name:ident) => {
        $crate::tensor_symbol!(stringify!($name))
    };
    ($name:ident; $($attr:ident),*) => {
        $crate::tensor_symbol!(stringify!($name); $($attr),*)
    };
    ($name:ident, $($setting:ident = $value:expr),*) => {
        $crate::tensor_symbol!(stringify!($name), $($setting = $value),*)
    };
    ($name:ident; $($attr:ident),+; $($setting:ident = $value:expr),*) => {
        $crate::tensor_symbol!(stringify!($name); $($attr),+; $($setting = $value),*)
    };
    ($id:expr) => {
        symbolica::symbol!($id, tag = &$crate::network::tags::SPENSO_TAG.tensor)
    };
    ($id:expr, tag = $tag:expr $(, $($rest:tt)*)?) => {
        compile_error!("tensor_symbol! owns the Spenso tensor tag; do not pass tag = ...")
    };
    ($id:expr, tags = $tags:expr $(, $($rest:tt)*)?) => {
        compile_error!("tensor_symbol! owns the Spenso tensor tag; do not pass tags = ...")
    };
    ($id:expr, $($setting:ident = $value:expr),*) => {
        symbolica::symbol!(
            $id,
            tag = &$crate::network::tags::SPENSO_TAG.tensor,
            $($setting = $value),*
        )
    };
    ($id:expr; $($attr:ident),*) => {
        symbolica::symbol!($id; $($attr),*; tag = &$crate::network::tags::SPENSO_TAG.tensor)
    };
    ($id:expr; $($attr:ident),+; tag = $tag:expr $(, $($rest:tt)*)?) => {
        compile_error!("tensor_symbol! owns the Spenso tensor tag; do not pass tag = ...")
    };
    ($id:expr; $($attr:ident),+; tags = $tags:expr $(, $($rest:tt)*)?) => {
        compile_error!("tensor_symbol! owns the Spenso tensor tag; do not pass tags = ...")
    };
    ($id:expr; $($attr:ident),+; $($setting:ident = $value:expr),*) => {
        symbolica::symbol!(
            $id;
            $($attr),+;
            tag = &$crate::network::tags::SPENSO_TAG.tensor,
            $($setting = $value),*
        )
    };
}

/// Creates a tensor-head symbol tagged as a Spenso vector.
///
/// This is the rank-one tensor head constructor. It expands
/// `symbolica::symbol!` at the call site, so `vector_symbol!(p)` gets the
/// caller's namespace plus the `tensor` and `rank1` tags.
#[macro_export]
macro_rules! vector_symbol {
    ($name:ident) => {
        symbolica::symbol!(
            stringify!($name),
            tags = [
                &$crate::network::tags::SPENSO_TAG.tensor,
                &$crate::network::tags::SPENSO_TAG.rank1
            ]
        )
    };
    ($name:ident, $($setting:ident = $value:expr),* $(,)?) => {
        symbolica::symbol!(
            stringify!($name),
            tags = [
                &$crate::network::tags::SPENSO_TAG.tensor,
                &$crate::network::tags::SPENSO_TAG.rank1
            ],
            $($setting = $value),*
        )
    };
    ($name:literal) => {
        symbolica::symbol!(
            $name,
            tags = [
                &$crate::network::tags::SPENSO_TAG.tensor,
                &$crate::network::tags::SPENSO_TAG.rank1
            ]
        )
    };
    ($name:literal, $($setting:ident = $value:expr),* $(,)?) => {
        symbolica::symbol!(
            $name,
            tags = [
                &$crate::network::tags::SPENSO_TAG.tensor,
                &$crate::network::tags::SPENSO_TAG.rank1
            ],
            $($setting = $value),*
        )
    };
}

/// Creates a representation symbol tagged with Spenso's representation tag.
///
/// This expands `symbolica::symbol!` at the call site and adds only the generic
/// representation tag.
#[macro_export]
macro_rules! representation_symbol {
    ($name:ident) => {
        symbolica::symbol!(
            stringify!($name),
            tag = &$crate::network::tags::SPENSO_TAG.representation
        )
    };
    ($name:literal) => {
        symbolica::symbol!(
            $name,
            tag = &$crate::network::tags::SPENSO_TAG.representation
        )
    };
}

/// Creates a self-dual representation symbol.
///
/// This expands `symbolica::symbol!` at the call site and adds the
/// `representation` and `self_dual` tags.
#[macro_export]
macro_rules! self_dual_symbol {
    ($name:ident) => {
        symbolica::symbol!(
            stringify!($name),
            tags = [
                &$crate::network::tags::SPENSO_TAG.representation,
                &$crate::network::tags::SPENSO_TAG.self_dual
            ]
        )
    };
    ($name:literal) => {
        symbolica::symbol!(
            $name,
            tags = [
                &$crate::network::tags::SPENSO_TAG.representation,
                &$crate::network::tags::SPENSO_TAG.self_dual
            ]
        )
    };
}

/// Creates a dualizable representation symbol.
///
/// This expands `symbolica::symbol!` at the call site and adds the
/// `representation` and `dualizable` tags.
#[macro_export]
macro_rules! dualizable_symbol {
    ($name:ident) => {
        symbolica::symbol!(
            stringify!($name),
            tags = [
                &$crate::network::tags::SPENSO_TAG.representation,
                &$crate::network::tags::SPENSO_TAG.dualizable
            ]
        )
    };
    ($name:literal) => {
        symbolica::symbol!(
            $name,
            tags = [
                &$crate::network::tags::SPENSO_TAG.representation,
                &$crate::network::tags::SPENSO_TAG.dualizable
            ]
        )
    };
}

/// Creates an abstract-index symbol tagged with Spenso's index tag.
///
/// This expands `symbolica::symbol!` at the call site and adds the Spenso
/// index tag.
#[macro_export]
macro_rules! index_symbol {
    ($name:ident) => {
        symbolica::symbol!(
            stringify!($name),
            tag = &$crate::network::tags::SPENSO_TAG.index
        )
    };
    ($name:literal) => {
        symbolica::symbol!($name, tag = &$crate::network::tags::SPENSO_TAG.index)
    };
}

/// Creates a function symbol tagged with Spenso's broadcast tag.
///
/// This expands `symbolica::symbol!` at the call site and adds the Spenso
/// broadcast tag.
#[macro_export]
macro_rules! broadcast_symbol {
    ($name:ident) => {
        symbolica::symbol!(
            stringify!($name),
            tag = &$crate::network::tags::SPENSO_TAG.broadcast
        )
    };
    ($name:literal) => {
        symbolica::symbol!($name, tag = &$crate::network::tags::SPENSO_TAG.broadcast)
    };
}

impl SpensoTags {
    fn print_chain_factor(value: AtomView<'_>, opt: &PrintOptions) -> Option<String> {
        fn without_visible_placeholders(value: AtomView<'_>, opt: &PrintOptions) -> Option<Atom> {
            let mut output = String::new();
            value.format(&mut output, opt, PrintState::new()).ok()?;
            let leaks_placeholder =
                [SPENSO_TAG.chain_in, SPENSO_TAG.chain_out]
                    .into_iter()
                    .any(|placeholder| {
                        let mut needle = String::new();
                        Atom::var(placeholder)
                            .as_view()
                            .format(&mut needle, opt, PrintState::new())
                            .is_ok()
                            && output.match_indices(&needle).any(|(start, _)| {
                                let end = start + needle.len();
                                let boundary = |character: char| {
                                    !character.is_alphanumeric() && character != '_'
                                };
                                output[..start].chars().next_back().is_none_or(boundary)
                                    && output[end..].chars().next().is_none_or(boundary)
                            })
                    });
            if !leaks_placeholder {
                return Some(value.to_owned());
            }

            match value {
                AtomView::Add(add) => add.iter().try_fold(Atom::Zero, |sum, term| {
                    Some(sum + without_visible_placeholders(term, opt)?)
                }),
                AtomView::Mul(mul) => mul.iter().try_fold(Atom::num(1), |product, factor| {
                    Some(product * without_visible_placeholders(factor, opt)?)
                }),
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    Some(without_visible_placeholders(base, opt)?.pow(exponent.to_owned()))
                }
                AtomView::Fun(function) => {
                    let arguments = function
                        .iter()
                        .filter(|argument| {
                            !matches!(argument, AtomView::Var(var) if var.get_symbol() == SPENSO_TAG.chain_in || var.get_symbol() == SPENSO_TAG.chain_out)
                        })
                        .map(|argument| without_visible_placeholders(argument, opt))
                        .collect::<Option<Vec<_>>>()?;
                    Some(if arguments.is_empty() {
                        Atom::var(function.get_symbol())
                    } else {
                        FunctionBuilder::new(function.get_symbol())
                            .add_args(arguments)
                            .finish()
                    })
                }
                _ => None,
            }
        }

        let sanitized = without_visible_placeholders(value, opt)?;
        let mut output = String::new();
        sanitized
            .as_view()
            .format(&mut output, opt, PrintState::new())
            .ok()?;
        Some(output)
    }

    fn print_dot(a: AtomView<'_>, opt: &PrintOptions, _state: &PrintState) -> Option<String> {
        let resolved = SpensoPrintSettings::resolve(opt)?;
        let parens = resolved.presentation.parens;
        let AtomView::Fun(f) = a else {
            return None;
        };

        if f.get_nargs() != 2 {
            return None;
        }
        let mut argitem = f.iter();
        let a = argitem.next()?;
        let b = argitem.next()?;

        let mut out = String::new();
        if parens {
            out.push('(');
        }
        a.format(&mut out, opt, PrintState::new()).ok()?;
        out.push_str(match resolved.backend {
            SpensoPrintBackend::Plain => ".",
            SpensoPrintBackend::Typst => " dot ",
            SpensoPrintBackend::Latex => r"\cdot ",
        });
        b.format(&mut out, opt, PrintState::new()).ok()?;
        if parens {
            out.push(')');
        }
        Some(out)
    }

    fn new() -> Self {
        let broadcast = tag!("broadcast");
        let upper = tag!("upper");
        let lower = tag!("lower");
        let rank1 = tag!("rank1");
        let tensor = tag!("tensor");
        let index = tag!("index");
        let representation = tag!("representation");
        let self_dual = tag!("self_dual");
        let dualizable = tag!("dualizable");
        Self {
            chain_in: symbol!("in"),
            chain_out: symbol!("out"),
            chain: symbol!(
                "chain";Linear;
                print = |a, opt, _state| {
                    let resolved = SpensoPrintSettings::resolve(opt)?;
                    let parens = resolved.presentation.parens;
                    let AtomView::Fun(f) = a else {
                        return None;
                    };
                    if f.get_nargs() < 2 {
                        return None;
                    }

                    let mut args = f.iter();
                    let in_index = args.next()?;
                    let out_index = args.next()?;

                    let mut s = String::new();
                    if !matches!(in_index, AtomView::Var(v) if v.get_symbol() == SPENSO_TAG.chain_in) {
                            in_index.format(&mut s, opt, PrintState::new()).unwrap();
                    }
                    if parens {
                        s.push('[');
                    }
                    for a in args {
                        s.push_str(&Self::print_chain_factor(a, opt)?);
                    }
                    if parens {
                        s.push(']');
                    }
                    if !matches!(out_index, AtomView::Var(v) if v.get_symbol() == SPENSO_TAG.chain_out) {
                            out_index.format(&mut s, opt, PrintState::new()).unwrap();
                    }
                    Some(s)
                }
            ),
            trace: symbol!(
                "trace";Linear;
                print = |a, opt, _state| {
                    let resolved = SpensoPrintSettings::resolve(opt)?;
                    let SpensoPrintSettings {
                        parens, with_dim, ..
                    } = resolved.presentation;
                    let AtomView::Fun(f) = a else {
                        return None;
                    };
                    if !(1..=2).contains(&f.get_nargs()) {
                        return None;
                    }

                    let mut args = f.iter();
                    let rep = args.next()?;
                    let mut s = match resolved.backend {
                        SpensoPrintBackend::Typst => r#"op("Tr")"#,
                        SpensoPrintBackend::Latex => r#"\operatorname{Tr}"#,
                        SpensoPrintBackend::Plain => "Tr",
                    }
                    .to_string();
                            if with_dim {
                                rep.format(&mut s, opt, PrintState::new()).unwrap();
                            }
                            if parens {
                                s.push('(');
                            }
                            match args.next() {
                                None => s.push('1'),
                                // Non-empty traces store factors under the canonical cyclic head.
                                Some(AtomView::Fun(cyclic))
                                    if cyclic.get_symbol() == *crate::shadowing::CYCLIC =>
                                {
                                    for factor in cyclic.iter() {
                                        s.push_str(&Self::print_chain_factor(factor, opt)?);
                                    }
                                }
                                Some(_) => return None,
                            }

                            if parens {
                                s.push(')');
                            }
                            Some(s)
                }
            ),
            rank1_: symbol!("rank1_", tags = [&tensor, &rank1]),
            bracket: symbol!(
                "bracket",
                print = |a, opt, state| {
                    SpensoPrintSettings::resolve(opt)?;
                    let AtomView::Fun(f) = a else {
                        return None;
                    };
                    if f.get_nargs() == 0 {
                        return None;
                    }

                    let parenthesize = state.in_exp || state.in_exp_base;
                    let mut product_state = *state;
                    product_state.in_product = true;
                    product_state.in_sum = false;
                    product_state.level += 1;

                    let mut output = String::new();
                    if parenthesize {
                        output.push('(');
                    }
                    for (position, factor) in f.iter().enumerate() {
                        if position > 0 {
                            if opt.mode.is_latex() {
                                output.push(' ');
                            } else {
                                output.push(opt.multiplication_operator);
                            }
                        }
                        factor.format(&mut output, opt, product_state).ok()?;
                    }
                    if parenthesize {
                        output.push(')');
                    }
                    Some(output)
                }
            ),
            pure_scalar: symbol!("pure_scalar"),
            scalar: symbol!("scalar"),
            dot: symbol!("dot";Symmetric,Linear; print = Self::print_dot),
            tensor_: symbol!("tensor_", tag = tensor),
            i_: symbol!("i_", tag = &index),
            rep_: symbol!("rep_", tag = &representation),
            self_dual_: symbol!("self_dual_", tags = [&representation, &self_dual]),
            dualizable_: symbol!("dualizable_", tags = [&representation, &dualizable]),
            broadcast,
            upper,
            lower,
            rank1,
            tensor,
            index,
            representation,
            self_dual,
            dualizable,
        }
    }

    define_numbered_tag_family_methods! {
        pub fn rank1_ => rank1_, "rank1", rank_one_tensor_symbol;
        pub fn rep_ => rep_, "rep", representation_symbol;
    }

    pub fn tensor_symbol(&self, name: &str) -> Symbol {
        symbol!(name, tag = &self.tensor)
    }

    pub fn representation_symbol(&self, name: &str) -> Symbol {
        symbol!(name, tag = &self.representation)
    }

    pub fn self_dual_symbol(&self, name: &str) -> Symbol {
        symbol!(name, tags = [&self.representation, &self.self_dual])
    }

    pub fn dualizable_symbol(&self, name: &str) -> Symbol {
        symbol!(name, tags = [&self.representation, &self.dualizable])
    }

    pub fn rank_one_tensor_symbol(&self, name: &str) -> Symbol {
        symbol!(name, tags = [&self.tensor, &self.rank1])
    }

    define_numbered_tag_family_methods! {
        pub fn tensor_ => tensor_, "tensor", tensor_symbol;
    }

    pub fn chain<'a, 'b, 'c, A, B, F>(
        &self,
        start: A,
        end: B,
        factors: impl IntoIterator<Item = F>,
    ) -> Atom
    where
        A: Into<AtomOrView<'a>>,
        B: Into<AtomOrView<'b>>,
        F: Into<AtomOrView<'c>>,
    {
        let mut f = FunctionBuilder::new(self.chain).add_arg(start).add_arg(end);
        for factor in factors {
            f = f.add_arg(factor);
        }
        f.finish()
    }

    pub fn trace<'a, 'b, R, F>(&self, rep: R, factors: impl IntoIterator<Item = F>) -> Atom
    where
        R: Into<AtomOrView<'a>>,
        F: Into<AtomOrView<'b>>,
    {
        let mut f = FunctionBuilder::new(self.trace).add_arg(rep);
        for factor in factors {
            f = f.add_arg(factor);
        }
        f.finish()
    }

    pub fn reverse_flip_factor(&self, factor: AtomView<'_>) -> Atom {
        let tmp = symbol!("spenso::chain_flip_tmp");
        factor
            .to_owned()
            .replace(self.chain_in)
            .with(tmp)
            .replace(self.chain_out)
            .with(self.chain_in)
            .replace(tmp)
            .with(self.chain_out)
    }

    pub fn reverse_flip_factors(&self, factors: impl IntoIterator<Item = Atom>) -> Vec<Atom> {
        let mut factors = factors.into_iter().collect::<Vec<_>>();
        factors.reverse();
        factors
            .into_iter()
            .map(|factor| self.reverse_flip_factor(factor.as_view()))
            .collect()
    }

    define_numbered_tag_family_methods! {
        pub fn self_dual_ => self_dual_, "self_dual", self_dual_symbol;
        pub fn dualizable_ => dualizable_, "dualizable", dualizable_symbol;
    }

    pub fn dualizable_dual_<'a, const N: usize, A: Into<AtomOrView<'a>>>(
        &self,
        args: impl IntoIterator<Item = A>,
    ) -> Atom {
        AIND_SYMBOLS.dual(self.dualizable_::<N, A>(args))
    }
}

#[cfg(test)]
mod tests {
    use symbolica::{
        atom::{Atom, AtomCore, AtomView, FunctionBuilder},
        printer::PrintOptions,
        symbol,
    };

    use crate::{cyclic, shadowing::symbolica_utils::SpensoPrintSettings};

    use super::{SPENSO_TAG, SymbolAtomExt};

    #[test]
    fn numbered_wildcard_macros_build_variables_without_args() {
        let expr = rank1_!(0);
        let AtomView::Var(var) = expr.as_view() else {
            panic!("empty wildcard head should be a variable");
        };

        assert_eq!(var.get_symbol(), SPENSO_TAG.rank1_);
    }

    #[test]
    fn numbered_wildcard_macros_build_functions_with_args() {
        let expr = rank1_!(
            1;
            Atom::var(symbol!("a___")),
            rep_!(2; Atom::var(symbol!("d_")))
        );

        let AtomView::Fun(fun) = expr.as_view() else {
            panic!("wildcard head with args should be a function");
        };

        assert_eq!(
            fun.get_symbol(),
            SPENSO_TAG.rank_one_tensor_symbol("rank11_")
        );
        assert_eq!(fun.get_nargs(), 2);
    }

    #[test]
    fn numbered_representation_families_use_their_own_prefixes() {
        let self_dual = self_dual_!(1; Atom::var(symbol!("d_")));
        let dualizable = dualizable_!(1; Atom::var(symbol!("d_")));

        let AtomView::Fun(self_dual) = self_dual.as_view() else {
            panic!("self-dual wildcard should be a function");
        };
        let AtomView::Fun(dualizable) = dualizable.as_view() else {
            panic!("dualizable wildcard should be a function");
        };

        assert_eq!(
            self_dual.get_symbol(),
            SPENSO_TAG.self_dual_symbol("self_dual1_")
        );
        assert_eq!(
            dualizable.get_symbol(),
            SPENSO_TAG.dualizable_symbol("dualizable1_")
        );
    }

    #[test]
    fn symbol_atom_ext_uses_variable_for_empty_args() {
        let symbol = SPENSO_TAG.tensor_symbol("empty_tensor_pattern");

        assert_eq!(
            symbol.atom_with_args(std::iter::empty::<Atom>()),
            Atom::var(symbol)
        );
        assert!(matches!(
            tensor_!(0; Atom::var(symbol!("a___"))).as_view(),
            AtomView::Fun(_)
        ));
    }

    #[test]
    fn trace_uses_typst_operator_only_in_typst_print_mode() {
        let trace = SPENSO_TAG.trace(Atom::var(symbol!("rep")), [cyclic!(Atom::num(1))]);

        assert_eq!(
            trace
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            r#"op("Tr")(1)"#
        );
        assert_eq!(
            trace
                .printer(SpensoPrintSettings::compact().nice_symbolica())
                .to_string(),
            "Tr(1)"
        );
        assert_eq!(
            trace
                .printer(SpensoPrintSettings::typst().nice_symbolica())
                .to_string(),
            "Tr(1)"
        );
    }

    #[test]
    fn chain_hides_namespaced_endpoint_placeholders() {
        assert_eq!(SPENSO_TAG.chain_in.get_namespace(), "spenso");
        assert_eq!(SPENSO_TAG.chain_out.get_namespace(), "spenso");

        let chain = SPENSO_TAG.chain(
            Atom::var(SPENSO_TAG.chain_in),
            Atom::var(SPENSO_TAG.chain_out),
            [Atom::var(SPENSO_TAG.tensor_symbol("factor"))],
        );

        let compact = chain
            .printer(SpensoPrintSettings::compact().nice_symbolica())
            .to_string();
        let typst = chain.printer(PrintOptions::typst()).to_string();

        assert_eq!(compact, "[factor]");
        assert_eq!(typst, "[\"factor\"]");
        assert!(!compact.contains("in"));
        assert!(!compact.contains("out"));
    }

    #[test]
    fn reverse_flip_swaps_namespaced_endpoint_placeholders() {
        let factor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("flip_factor"))
            .add_arg(Atom::var(SPENSO_TAG.chain_in))
            .add_arg(Atom::var(SPENSO_TAG.chain_out))
            .finish();
        let flipped = SPENSO_TAG.reverse_flip_factor(factor.as_view());
        let AtomView::Fun(flipped) = flipped.as_view() else {
            panic!("expected a tensor factor")
        };
        let arguments = flipped.iter().collect::<Vec<_>>();

        assert!(
            matches!(arguments[0], AtomView::Var(var) if var.get_symbol() == SPENSO_TAG.chain_out)
        );
        assert!(
            matches!(arguments[1], AtomView::Var(var) if var.get_symbol() == SPENSO_TAG.chain_in)
        );
    }

    #[test]
    fn bracket_prints_as_an_ordered_product() {
        let bracket = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(Atom::var(symbol!("z")))
            .add_arg(Atom::var(symbol!("a")))
            .finish();
        let mut compact = SpensoPrintSettings::compact().nice_symbolica();
        compact.color_builtin_symbols = false;

        assert_eq!(bracket.printer(compact).to_string(), "z·a");
        assert_eq!(bracket.printer(PrintOptions::typst()).to_string(), "z a");
    }
}
