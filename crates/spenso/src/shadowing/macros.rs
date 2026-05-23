#![allow(unused_macros)]

/// Creates a Symbolica symbol from an identifier.
///
/// This is a small shorthand for `symbolica::symbol!(stringify!(name))`.
/// It is useful in expression-building code where many symbolic abstract
/// indices are needed.
#[macro_export]
macro_rules! s {
    ($name:ident) => {
        symbolica::symbol!(stringify!($name))
    };
}

/// Creates an abstract-index slot from a representation and an index.
///
/// The slot index is fixed to Spenso's symbolic `AbstractIndex`, which avoids
/// type-inference ambiguity when passing `Symbol`s. The second argument can be
/// an identifier, expanded through [`s!`], an integer index, or an arbitrary
/// expression convertible into `AbstractIndex`.
///
/// # Examples
///
/// ```ignore
/// use spenso::{s, slot};
/// use spenso::structure::representation::{Minkowski, RepName};
///
/// let mink4 = Minkowski {}.new_rep(4);
/// let mu = slot!(mink4, mu);
/// let nu = slot!(mink4, s!(nu));
/// let one = slot!(mink4, 1);
/// ```
#[macro_export]
macro_rules! slot {
    ($rep:expr, $index:ident) => {
        ($rep).slot::<$crate::structure::abstract_index::AbstractIndex, _>($crate::s!($index))
    };
    ($rep:expr, $index:expr) => {
        ($rep).slot::<$crate::structure::abstract_index::AbstractIndex, _>($index)
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! spenso_rep_atom {
    ($rep:expr; $dim:ident $(,)?) => {{
        let rep = $crate::structure::representation::RepName::new_rep(&$rep, $crate::s!($dim));
        rep.to_symbolic(std::iter::empty::<symbolica::atom::Atom>())
    }};
    ($rep:expr; $dim:expr $(,)?) => {{
        let rep = $crate::structure::representation::RepName::new_rep(&$rep, $dim);
        rep.to_symbolic(std::iter::empty::<symbolica::atom::Atom>())
    }};
    ($rep:expr; $dim:ident, $index:ident $(,)?) => {{
        let rep = $crate::structure::representation::RepName::new_rep(&$rep, $crate::s!($dim));
        rep.to_symbolic([$crate::shadowing::IntoAtom::into_atom($crate::s!($index))])
    }};
    ($rep:expr; $dim:ident, $index:expr $(,)?) => {{
        let rep = $crate::structure::representation::RepName::new_rep(&$rep, $crate::s!($dim));
        rep.to_symbolic([$crate::shadowing::IntoAtom::into_atom($index)])
    }};
    ($rep:expr; $dim:expr, $index:ident $(,)?) => {{
        let rep = $crate::structure::representation::RepName::new_rep(&$rep, $dim);
        rep.to_symbolic([$crate::shadowing::IntoAtom::into_atom($crate::s!($index))])
    }};
    ($rep:expr; $dim:expr, $index:expr $(,)?) => {{
        let rep = $crate::structure::representation::RepName::new_rep(&$rep, $dim);
        rep.to_symbolic([$crate::shadowing::IntoAtom::into_atom($index)])
    }};
}

/// Builds a stripped or indexed Minkowski representation atom.
#[macro_export]
macro_rules! mink {
    ($($args:tt)*) => {
        $crate::spenso_rep_atom!($crate::structure::representation::Minkowski {}; $($args)*)
    };
}

/// Builds a stripped or indexed Euclidean representation atom.
#[macro_export]
macro_rules! euc {
    ($($args:tt)*) => {
        $crate::spenso_rep_atom!($crate::structure::representation::Euclidean {}; $($args)*)
    };
}

/// Builds a stripped or indexed Lorentz representation atom.
#[macro_export]
macro_rules! lor {
    ($($args:tt)*) => {
        $crate::spenso_rep_atom!($crate::structure::representation::Lorentz {}; $($args)*)
    };
}

/// Builds a symbolic tensor function from a name and tensor arguments.
///
/// Identifier heads are created with [`tensor_symbol!`], so they keep the
/// caller's symbol namespace and carry the generic Spenso tensor tag. Pass an
/// explicit `Symbol` expression when the head was built elsewhere.
///
/// Arguments are converted through [`IntoAtom`], so callers can mix scalar
/// literal arguments, atoms, slots, and stripped representations. Passing a
/// representation, for example `tensor!(p, rep)`, emits compact Schoonschip
/// syntax; passing a slot, for example `tensor!(p, slot!(rep, i))`, emits an
/// explicitly indexed tensor.
///
/// # Examples
///
/// ```ignore
/// use spenso::{p, q, slot, tensor};
///
/// let compact = tensor!(k, 1, mink4);
/// let indexed = p!(1, slot!(mink4, mu)) * q!(2, slot!(mink4, mu));
/// ```
#[macro_export]
macro_rules! tensor {
    ($name:ident $(,)?) => {
        symbolica::atom::FunctionBuilder::new($crate::tensor_symbol!($name)).finish()
    };
    ($name:ident, $($arg:expr),+ $(,)?) => {{
        let mut tensor = symbolica::atom::FunctionBuilder::new($crate::tensor_symbol!($name));
        $(
            tensor = tensor.add_arg($crate::shadowing::IntoAtom::into_atom($arg));
        )*
        tensor.finish()
    }};
    ($name:expr $(,)?) => {
        symbolica::atom::FunctionBuilder::new($name).finish()
    };
    ($name:expr, $($arg:expr),+ $(,)?) => {{
        let mut tensor = symbolica::atom::FunctionBuilder::new($name);
        $(
            tensor = tensor.add_arg($crate::shadowing::IntoAtom::into_atom($arg));
        )*
        tensor.finish()
    }};
}

/// Builds a symbolic vector function from a name and tensor arguments.
///
/// The head is created with [`vector_symbol!`], so `vector!(p, ...)` is a
/// rank-one tensor in the caller's symbol namespace.
#[macro_export]
macro_rules! vector {
    ($name:ident $(,)?) => {
        symbolica::atom::FunctionBuilder::new($crate::vector_symbol!($name)).finish()
    };
    ($name:ident, $($arg:expr),+ $(,)?) => {{
        let mut tensor = symbolica::atom::FunctionBuilder::new($crate::vector_symbol!($name));
        $(
            tensor = tensor.add_arg($crate::shadowing::IntoAtom::into_atom($arg));
        )*
        tensor.finish()
    }};
    ($name:expr $(,)?) => {
        symbolica::atom::FunctionBuilder::new($name).finish()
    };
    ($name:expr, $($arg:expr),+ $(,)?) => {{
        let mut tensor = symbolica::atom::FunctionBuilder::new($name);
        $(
            tensor = tensor.add_arg($crate::shadowing::IntoAtom::into_atom($arg));
        )*
        tensor.finish()
    }};
}

/// Builds a symbolic vector `p(...)`.
///
/// This is a convenience wrapper around [`vector!`] with a rank-one tensor head.
#[macro_export]
macro_rules! p {
    () => {
        $crate::vector!(p)
    };
    ($($arg:expr),+ $(,)?) => {
        $crate::vector!(p, $($arg),+)
    };
}

/// Builds a symbolic vector `q(...)`.
///
/// This is a convenience wrapper around [`vector!`] with a rank-one tensor head.
#[macro_export]
macro_rules! q {
    () => {
        $crate::vector!(q)
    };
    ($($arg:expr),+ $(,)?) => {
        $crate::vector!(q, $($arg),+)
    };
}

/// Builds a metric tensor atom `g(a,b)`.
#[macro_export]
macro_rules! metric {
    ($a:expr, $b:expr $(,)?) => {
        symbolica::function!(
            $crate::network::library::symbolic::ETS.metric,
            $crate::shadowing::IntoAtom::into_atom($a),
            $crate::shadowing::IntoAtom::into_atom($b)
        )
    };
}

/// Builds a metric tensor atom `g(a,b)`.
#[macro_export]
macro_rules! g {
    ($a:expr, $b:expr $(,)?) => {
        $crate::metric!($a, $b)
    };
}

/// Builds a two-argument compact dot-product atom.
#[macro_export]
macro_rules! dot {
    ($a:expr, $b:expr $(,)?) => {
        symbolica::function!(
            $crate::network::tags::SPENSO_TAG.dot,
            $crate::shadowing::IntoAtom::into_atom($a),
            $crate::shadowing::IntoAtom::into_atom($b)
        )
    };
}

/// Forces an expression through the parser's pure-scalar boundary.
#[macro_export]
macro_rules! pure_scalar {
    ($expr:expr $(,)?) => {
        symbolica::function!(
            $crate::network::tags::SPENSO_TAG.pure_scalar,
            $crate::shadowing::IntoAtom::into_atom($expr)
        )
    };
}

/// Builds a parser grouping over network factors.
#[macro_export]
macro_rules! bracket {
    ($($expr:expr),* $(,)?) => {
        symbolica::function!(
            $crate::network::tags::SPENSO_TAG.bracket,
            $($crate::shadowing::IntoAtom::into_atom($expr)),*
        )
    };
}

/// Bundles structural slots inside one tensor argument.
#[macro_export]
macro_rules! aind {
    ($($slot:expr),* $(,)?) => {
        symbolica::function!(
            $crate::structure::abstract_index::AIND_SYMBOLS.aind,
            $($crate::shadowing::IntoAtom::into_atom($slot)),*
        )
    };
}

/// Wraps a slot or representation atom in the dual-index marker.
#[macro_export]
macro_rules! dind {
    ($slot:expr $(,)?) => {
        symbolica::function!(
            $crate::structure::abstract_index::AIND_SYMBOLS.dind,
            $crate::shadowing::IntoAtom::into_atom($slot)
        )
    };
}

/// Builds a generic chain or trace factor with explicit `in` and `out` markers.
///
/// Use the identifier tokens `in` and `out` in the argument list to insert the
/// chain placeholder symbols. Other arguments are converted through `IntoAtom`.
#[macro_export]
macro_rules! chain_factor {
    ($name:ident $(,)?) => {{
        let factor = symbolica::atom::FunctionBuilder::new($crate::tensor_symbol!($name));
        $crate::chain_factor!(@finish factor)
    }};
    ($name:ident, $($args:tt)+) => {{
        let factor = symbolica::atom::FunctionBuilder::new($crate::tensor_symbol!($name));
        $crate::chain_factor!(@finish factor, $($args)+)
    }};
    ($name:expr; $($args:tt)*) => {{
        let factor = symbolica::atom::FunctionBuilder::new($name);
        $crate::chain_factor!(@finish factor $(, $args)*)
    }};
    (@finish $factor:ident) => {
        $factor.finish()
    };
    (@finish $factor:ident, in $(, $rest:tt)*) => {{
        let $factor = $factor.add_arg(symbolica::atom::Atom::var($crate::network::tags::SPENSO_TAG.chain_in));
        $crate::chain_factor!(@finish $factor $(, $rest)*)
    }};
    (@finish $factor:ident, out $(, $rest:tt)*) => {{
        let $factor = $factor.add_arg(symbolica::atom::Atom::var($crate::network::tags::SPENSO_TAG.chain_out));
        $crate::chain_factor!(@finish $factor $(, $rest)*)
    }};
    (@finish $factor:ident, $arg:expr $(, $rest:tt)*) => {{
        let $factor = $factor.add_arg($crate::shadowing::IntoAtom::into_atom($arg));
        $crate::chain_factor!(@finish $factor $(, $rest)*)
    }};
}

/// Builds an inert normalized symmetrizer over chain/trace factors.
///
/// The resulting `sym(...)` is canonicalized through Symbolica's symmetric
/// function attribute. It expands only when explicitly passed through
/// [`ProjectorExpander::expand_projectors`].
///
/// # Examples
///
/// ```ignore
/// use spenso::{chain, sym};
///
/// let projected = chain!(i, j, sym!(a, b, c));
/// let from_vec = sym!(; vec![a, b, c]);
/// ```
#[macro_export]
macro_rules! sym {
    (; $factors:expr $(,)?) => {
        $crate::shadowing::sym(($factors).into_iter())
    };
    ($($factor:expr),* $(,)?) => {
        $crate::shadowing::sym(vec![$($crate::shadowing::IntoAtom::into_atom($factor)),*])
    };
}

/// Builds an inert normalized antisymmetrizer over chain/trace factors.
///
/// The resulting `antisym(...)` is canonicalized through Symbolica's
/// antisymmetric function attribute. It expands only when explicitly passed
/// through [`ProjectorExpander::expand_projectors`].
///
/// # Examples
///
/// ```ignore
/// use spenso::{antisym, chain};
///
/// let commutator_part = chain!(i, j, antisym!(a, b));
/// let from_vec = antisym!(; vec![a, b, c]);
/// ```
#[macro_export]
macro_rules! antisym {
    (; $factors:expr $(,)?) => {
        $crate::shadowing::antisym(($factors).into_iter())
    };
    ($($factor:expr),* $(,)?) => {
        $crate::shadowing::antisym(vec![$($crate::shadowing::IntoAtom::into_atom($factor)),*])
    };
}

/// Builds an inert normalized cyclic symmetrizer over chain/trace factors.
///
/// The resulting `cyclic(...)` is canonicalized through Symbolica's
/// cyclesymmetric function attribute. Inside traces it is the compact syntax
/// for cyclic trace equivalence, for example `trace(rep, cyclic(a,b,c))`.
///
/// # Examples
///
/// ```ignore
/// use spenso::{cyclic, trace};
///
/// let cyclic_trace = trace!(rep, cyclic!(a, b, c));
/// let from_vec = cyclic!(; vec![a, b, c]);
/// ```
#[macro_export]
macro_rules! cyclic {
    (; $factors:expr $(,)?) => {
        $crate::shadowing::cyclic(($factors).into_iter())
    };
    ($($factor:expr),* $(,)?) => {
        $crate::shadowing::cyclic(vec![$($crate::shadowing::IntoAtom::into_atom($factor)),*])
    };
}

/// Builds a symbolic open chain with a start slot, an end slot, and any number
/// of factors.
///
/// Arguments are converted through [`IntoAtom`], so callers can pass Spenso
/// slots, atoms, atom views, or symbols. This macro is only the variadic surface
/// for [`SPENSO_TAG.chain`](crate::network::tags::SPENSO_TAG); it does not
/// simplify or normalize the expression.
///
/// # Examples
///
/// ```ignore
/// use spenso::{chain, slot};
///
/// let expr = chain!(
///     slot!(bis4, a),
///     slot!(bis4, b),
///     gamma_mu,
///     gamma_nu,
/// );
///
/// let factors = vec![gamma_mu, gamma_nu];
/// let expr = chain!(slot!(bis4, a), slot!(bis4, b); factors);
/// ```
#[macro_export]
macro_rules! chain {
    ($start:expr, $end:expr $(,)?) => {
        $crate::network::tags::SPENSO_TAG.chain(
            $crate::shadowing::IntoAtom::into_atom($start),
            $crate::shadowing::IntoAtom::into_atom($end),
            std::iter::empty::<symbolica::atom::Atom>(),
        )
    };
    ($start:expr, $end:expr; $factors:expr $(,)?) => {
        $crate::network::tags::SPENSO_TAG.chain(
            $crate::shadowing::IntoAtom::into_atom($start),
            $crate::shadowing::IntoAtom::into_atom($end),
            ($factors)
                .into_iter()
                .map($crate::shadowing::IntoAtom::into_atom),
        )
    };
    ($start:expr, $end:expr $(, $factor:expr)+; $factors:expr $(,)?) => {
        $crate::network::tags::SPENSO_TAG.chain(
            $crate::shadowing::IntoAtom::into_atom($start),
            $crate::shadowing::IntoAtom::into_atom($end),
            vec![$($crate::shadowing::IntoAtom::into_atom($factor)),*]
                .into_iter()
                .chain(($factors).into_iter().map($crate::shadowing::IntoAtom::into_atom)),
        )
    };
    ($start:expr, $end:expr $(, $factor:expr)* $(,)?) => {
        $crate::network::tags::SPENSO_TAG.chain(
            $crate::shadowing::IntoAtom::into_atom($start),
            $crate::shadowing::IntoAtom::into_atom($end),
            vec![$($crate::shadowing::IntoAtom::into_atom($factor)),*],
        )
    };
}

/// Builds a symbolic trace with a representation and any number of factors.
///
/// The representation and factors are converted through [`IntoAtom`]. Empty
/// traces are emitted as `trace(rep)`. Non-empty traces are emitted in the
/// canonical cyclic form `trace(rep, cyclic(factors...))`.
///
/// # Examples
///
/// ```ignore
/// use spenso::{cyclic, trace};
///
/// let expr = trace!(
///     &cof_nc,
///     color_t_a,
///     color_t_b,
/// );
/// assert_eq!(expr, trace!(&cof_nc, cyclic!(color_t_a, color_t_b)));
///
/// let factors = vec![color_t_a, color_t_b];
/// let expr = trace!(&cof_nc; factors);
/// ```
#[macro_export]
macro_rules! trace {
    ($rep:expr $(,)?) => {
        $crate::shadowing::trace(
            $crate::shadowing::IntoAtom::into_atom($rep),
            std::iter::empty::<symbolica::atom::Atom>(),
        )
    };
    ($rep:expr; $factors:expr $(,)?) => {
        $crate::shadowing::trace(
            $crate::shadowing::IntoAtom::into_atom($rep),
            ($factors)
                .into_iter()
                .map($crate::shadowing::IntoAtom::into_atom),
        )
    };
    ($rep:expr $(, $factor:expr)+; $factors:expr $(,)?) => {
        $crate::shadowing::trace(
            $crate::shadowing::IntoAtom::into_atom($rep),
            vec![$($crate::shadowing::IntoAtom::into_atom($factor)),*]
                .into_iter()
                .chain(($factors).into_iter().map($crate::shadowing::IntoAtom::into_atom)),
        )
    };
    ($rep:expr $(, $factor:expr)* $(,)?) => {
        $crate::shadowing::trace(
            $crate::shadowing::IntoAtom::into_atom($rep),
            vec![$($crate::shadowing::IntoAtom::into_atom($factor)),*],
        )
    };
}

/// Builds a symbolic trace with a single symmetric projector payload.
///
/// This emits `trace(rep, sym(factors...))`. Use `trace!` for ordinary cyclic
/// closed traces.
///
/// # Examples
///
/// ```ignore
/// use spenso::trace_sym;
///
/// let invariant = trace_sym!(&cof_nc, color_t_a, color_t_b, color_t_c);
/// let from_vec = trace_sym!(&cof_nc; factors);
/// ```
#[macro_export]
macro_rules! trace_sym {
    ($rep:expr; $factors:expr $(,)?) => {
        $crate::shadowing::trace_sym(
            $crate::shadowing::IntoAtom::into_atom($rep),
            ($factors)
                .into_iter()
                .map($crate::shadowing::IntoAtom::into_atom),
        )
    };
    ($rep:expr $(, $factor:expr)+; $factors:expr $(,)?) => {
        $crate::shadowing::trace_sym(
            $crate::shadowing::IntoAtom::into_atom($rep),
            vec![$($crate::shadowing::IntoAtom::into_atom($factor)),*]
                .into_iter()
                .chain(($factors).into_iter().map($crate::shadowing::IntoAtom::into_atom)),
        )
    };
    ($rep:expr $(, $factor:expr)+ $(,)?) => {
        $crate::shadowing::trace_sym(
            $crate::shadowing::IntoAtom::into_atom($rep),
            vec![$($crate::shadowing::IntoAtom::into_atom($factor)),*],
        )
    };
}

#[macro_export]
macro_rules!  symbol_set {
        // Identifier symbols with an explicit namespace.
        ($struct_name:ident, $static_name:ident, namespace = $namespace:literal; $($char:ident)*) => {
            #[allow(dead_code, non_snake_case)]
            pub struct $struct_name {
                $(pub $char: ::symbolica::atom::Symbol,)*
            }

            pub static $static_name: ::std::sync::LazyLock<$struct_name> = ::std::sync::LazyLock::new(|| $struct_name {
                $($char: ::symbolica::symbol!(concat!($namespace, "::", stringify!($char))),)*
            });
        };

        // Identifier symbols with explicit struct and static names
        ($struct_name:ident, $static_name:ident; $($char:ident)*) => {
            #[allow(dead_code, non_snake_case)]
            pub struct $struct_name {
                $(pub $char: ::symbolica::atom::Symbol,)*
            }

            pub static $static_name: ::std::sync::LazyLock<$struct_name> = ::std::sync::LazyLock::new(|| $struct_name {
                $($char: ::symbolica::symbol!(stringify!($char)),)*
            });
        };

        // String literals as individual statics with an explicit namespace.
        (statics, namespace = $namespace:literal; $($field:ident : $string:literal),* $(,)?) => {
            $(
                #[allow(non_upper_case_globals)]
                pub static $field: ::std::sync::LazyLock<::symbolica::atom::Symbol> =
                    ::std::sync::LazyLock::new(|| ::symbolica::symbol!(concat!($namespace, "::", $string)));
            )*
        };

        // String literals as individual statics
        (statics; $($field:ident : $string:literal),* $(,)?) => {
            $(
                #[allow(non_upper_case_globals)]
                pub static $field: ::std::sync::LazyLock<::symbolica::atom::Symbol> =
                    ::std::sync::LazyLock::new(|| ::symbolica::symbol!($string));
            )*
        };

        // String literals grouped in a struct with an explicit namespace.
        ($struct_name:ident, $static_name:ident, namespace = $namespace:literal; $($field:ident : $string:literal),* $(,)?) => {
            #[allow(dead_code, non_snake_case)]
            pub struct $struct_name {
                $(pub $field: ::symbolica::atom::Symbol,)*
            }

            pub static $static_name: ::std::sync::LazyLock<$struct_name> = ::std::sync::LazyLock::new(|| $struct_name {
                $($field: ::symbolica::symbol!(concat!($namespace, "::", $string)),)*
            });
        };

        // String literals grouped in a struct
        ($struct_name:ident, $static_name:ident; $($field:ident : $string:literal),* $(,)?) => {
            #[allow(dead_code, non_snake_case)]
            pub struct $struct_name {
                $(pub $field: ::symbolica::atom::Symbol,)*
            }

            pub static $static_name: ::std::sync::LazyLock<$struct_name> = ::std::sync::LazyLock::new(|| $struct_name {
                $($field: ::symbolica::symbol!($string),)*
            });
        };
    }
