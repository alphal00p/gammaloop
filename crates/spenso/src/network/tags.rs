use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    printer::PrintState,
    symbol, tag,
};

use crate::{
    shadowing::symbolica_utils::SpensoPrintSettings, structure::abstract_index::AIND_SYMBOLS,
};

pub struct SpensoTags {
    pub broadcast: String,
    /// Marks rank-one tensor symbols whose final argument is the tensor slot.
    ///
    /// Symbols carrying this tag must not use representation slots in earlier
    /// arguments; rank-one shorthand rewrites assume that contract.
    pub rank1: String,
    pub rank1_: Symbol,
    pub chain_in: Symbol,
    pub chain_out: Symbol,
    pub chain: Symbol,
    pub trace: Symbol,
    pub upper: String,
    pub lower: String,
    pub bracket: Symbol,
    pub pure_scalar: Symbol,
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

pub static SPENSO_TAG: std::sync::LazyLock<SpensoTags> = std::sync::LazyLock::new(SpensoTags::new);

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
    ($name:literal) => {
        symbolica::symbol!(
            $name,
            tags = [
                &$crate::network::tags::SPENSO_TAG.tensor,
                &$crate::network::tags::SPENSO_TAG.rank1
            ]
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
                print = |a, opt| {
                    match opt.custom_print_mode {
                        Some(("spenso", i)) => {
                            let SpensoPrintSettings { parens, .. } = SpensoPrintSettings::from(i);

                            let AtomView::Fun(f) = a else {
                                return None;
                            };

                            let mut args = f.iter();

                            let in_index = args.next().unwrap();
                            let out_index = args.next().unwrap();

                            let mut s = String::new();
                            in_index.format(&mut s, opt, PrintState::new()).unwrap();
                            if parens {
                                s.push('[');
                            }
                            for a in args {
                                a.format(&mut s, opt, PrintState::new()).unwrap();
                            }
                            if parens {
                                s.push(']');
                            }
                            out_index.format(&mut s, opt, PrintState::new()).unwrap();
                            Some(s)
                        }
                        _ => None,
                    }
                }
            ),
            trace: symbol!(
                "trace";Linear;
                print = |a, opt| {
                    match opt.custom_print_mode {
                        Some(("spenso", i)) => {
                            let SpensoPrintSettings {
                                parens, with_dim, ..
                            } = SpensoPrintSettings::from(i);

                            let AtomView::Fun(f) = a else {
                                return None;
                            };

                            let mut args = f.iter();

                            let rep = args.next().unwrap();

                            let mut s = String::from("Tr");
                            if with_dim {
                                rep.format(&mut s, opt, PrintState::new()).unwrap();
                            }
                            if parens {
                                s.push('(');
                            }
                            for a in args {
                                a.format(&mut s, opt, PrintState::new()).unwrap();
                            }
                            if parens {
                                s.push(')');
                            }
                            Some(s)
                        }
                        _ => None,
                    }
                }
            ),
            rank1_: symbol!("rank1_", tags = [&tensor, &rank1]),
            bracket: symbol!("bracket"),
            pure_scalar: symbol!("pure_scalar"),
            dot: symbol!("dot";Symmetric,Linear; print = |a,opt|{
                    match opt.custom_print_mode {
                        Some(("spenso",i))=>{
                            let SpensoPrintSettings{
                                parens,
                                with_dim,..
                            } = SpensoPrintSettings::from(i);


                    let AtomView::Fun(f) = a else {
                        return None;
                    };

                    if f.get_nargs() != 3 {
                        return None;
                    }
                    let mut args = f.iter();

                    let a = args.next().unwrap();
                    let b = args.next().unwrap();
                    let c = args.next().unwrap();

                    fn is_rep(view:AtomView<'_>)->bool{
                        match view {
                            AtomView::Fun(f) if f.get_symbol().has_tag(&SPENSO_TAG.upper) => true,
                            AtomView::Var(s) if s.get_symbol().has_tag(&SPENSO_TAG.upper) => true,
                            _=>false
                        }
                    }

                    let (a,b,c) = if is_rep(a) && !is_rep(b) && !is_rep(c) {
                        (a,b,c)
                    } else if is_rep(b) && !is_rep(a) && !is_rep(c) {
                        (b,c,a)
                    } else if is_rep(c) && !is_rep(a) && !is_rep(b) {
                        (c,a,b)
                    } else { return None};

                    let mut s = String::new();
                    if parens {
                        s.push('(');
                    }
                    b.format(&mut s, opt,PrintState::new()).unwrap();
                    s.push('.');
                    if with_dim {a.format(&mut s, opt, PrintState::new()).unwrap();
                        s.push('.');
                    }
                    c.format(&mut s, opt,PrintState::new()).unwrap();
                    if parens {
                        s.push(')');
                    }
                    Some(s)

                },
                _=>None
            }



            }),
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
        atom::{Atom, AtomView},
        symbol,
    };

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
}
