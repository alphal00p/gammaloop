/// Builds an identity atom between two indices.
///
/// The arguments are parsed as Symbolica literals, preserving the historical
/// `id!(i, j)` shorthand used in idenso tests and examples.
#[macro_export]
macro_rules! id {
    ($i: expr, $j: expr) => {{
        let i = symbolica::parse_lit!($i);
        let j = symbolica::parse_lit!($j);
        $crate::dirac::id_atom(i, j)
    }};
}

/// Builds a gamma matrix.
///
/// With one argument, this builds a chain factor using the placeholder indices
/// `in` and `out`; use this form only as a factor inside `chain!` or `trace!`.
/// With three arguments, this builds the ordinary gamma tensor with explicit
/// spinor endpoints and a Lorentz slot.
///
/// Arguments are converted through `spenso::symbolica_atom::IntoAtom`, so they
/// can be typed slots, atoms, or atom views.
///
/// # Examples
///
/// ```ignore
/// use idenso::gamma;
/// use spenso::{chain, slot};
///
/// let factor = gamma!(slot!(mink4, mu));
/// let chain_expr = chain!(slot!(bis4, a), slot!(bis4, b), factor);
/// let default_tensor = gamma!(mu, a, b);
/// let indexed_tensor = gamma!(1, 2, 3);
/// let pattern_tensor = gamma!(RS.a__, RS.b__, RS.c__);
/// let pattern_tensor_with_slots = gamma!(RS.a__, [RS.d_, RS.i_], [RS.d_, RS.j_]);
/// let mixed_tensor = gamma!(mu, slot!(bis_d, a), 1);
/// let explicit_tensor = gamma!(slot!(mink_d, mu), slot!(bis_d, a), slot!(bis_d, b));
/// ```
#[macro_export]
macro_rules! gamma {
    ($mu:expr) => {
        symbolica::atom::FunctionBuilder::new($crate::dirac::AGS.gamma)
            .add_arg(symbolica::atom::Atom::var(spenso::network::tags::SPENSO_TAG.chain_in))
            .add_arg(symbolica::atom::Atom::var(spenso::network::tags::SPENSO_TAG.chain_out))
            .add_arg(spenso::symbolica_atom::IntoAtom::into_atom($mu))
            .finish()
    };
    ($base:ident . $mu:ident, $($rest:tt)*) => {
        $crate::gamma!(@tensor spenso::structure::representation::RepName::to_symbolic(
            &spenso::structure::representation::Minkowski {},
            [$base.$mu],
        ); $($rest)*)
    };
    ($mu:ident, $($rest:tt)*) => {{
        let mink = spenso::structure::representation::RepName::new_rep(
            &spenso::structure::representation::Minkowski {},
            4,
        );
        $crate::gamma!(@tensor spenso::slot!(mink, $mu); $($rest)*)
    }};
    ($mu:literal, $($rest:tt)*) => {{
        let mink = spenso::structure::representation::RepName::new_rep(
            &spenso::structure::representation::Minkowski {},
            4,
        );
        $crate::gamma!(@tensor spenso::slot!(mink, $mu); $($rest)*)
    }};
    ($mu:expr, $($rest:tt)*) => {
        $crate::gamma!(@tensor $mu; $($rest)*)
    };
    (@tensor $mu:expr; [$($i:expr),+ $(,)?], $($rest:tt)*) => {
        $crate::gamma!(@tensor2 $mu, spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::Bispinor {},
            [$($i),+],
        ); $($rest)*)
    };
    (@tensor $mu:expr; $base:ident . $i:ident, $($rest:tt)*) => {
        $crate::gamma!(@tensor2 $mu, spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::Bispinor {},
            [$base.$i],
        ); $($rest)*)
    };
    (@tensor $mu:expr; $i:ident, $($rest:tt)*) => {{
        let bis = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::Bispinor {},
            4,
        );
        $crate::gamma!(@tensor2 $mu, spenso::slot!(bis, $i); $($rest)*)
    }};
    (@tensor $mu:expr; $i:literal, $($rest:tt)*) => {{
        let bis = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::Bispinor {},
            4,
        );
        $crate::gamma!(@tensor2 $mu, spenso::slot!(bis, $i); $($rest)*)
    }};
    (@tensor $mu:expr; $i:expr, $($rest:tt)*) => {
        $crate::gamma!(@tensor2 $mu, $i; $($rest)*)
    };
    (@tensor2 $mu:expr, $i:expr; [$($j:expr),+ $(,)?]) => {
        $crate::gamma!(@tensor_done $mu, $i, spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::Bispinor {},
            [$($j),+],
        ))
    };
    (@tensor2 $mu:expr, $i:expr; $base:ident . $j:ident) => {
        $crate::gamma!(@tensor_done $mu, $i, spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::Bispinor {},
            [$base.$j],
        ))
    };
    (@tensor2 $mu:expr, $i:expr; $j:ident) => {{
        let bis = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::Bispinor {},
            4,
        );
        $crate::gamma!(@tensor_done $mu, $i, spenso::slot!(bis, $j))
    }};
    (@tensor2 $mu:expr, $i:expr; $j:literal) => {{
        let bis = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::Bispinor {},
            4,
        );
        $crate::gamma!(@tensor_done $mu, $i, spenso::slot!(bis, $j))
    }};
    (@tensor2 $mu:expr, $i:expr; $j:expr) => {
        $crate::gamma!(@tensor_done $mu, $i, $j)
    };
    (@tensor_done $mu:expr, $i:expr, $j:expr) => {
        symbolica::atom::FunctionBuilder::new($crate::dirac::AGS.gamma)
            .add_arg(spenso::symbolica_atom::IntoAtom::into_atom($i))
            .add_arg(spenso::symbolica_atom::IntoAtom::into_atom($j))
            .add_arg(spenso::symbolica_atom::IntoAtom::into_atom($mu))
            .finish()
    };
}

/// Builds a gamma-zero factor or tensor.
///
/// With no arguments, this builds a chain factor using the placeholder indices
/// `in` and `out`; use this form only as a factor inside `chain!` or `trace!`.
/// With two arguments, this builds the ordinary gamma-zero tensor between
/// explicit spinor endpoints.
///
/// Endpoints can be typed slots, atoms, atom views, default four-dimensional
/// bispinor identifiers/literals, pattern variables such as `RS.i__`, or
/// bracketed bispinor argument lists such as `[RS.d_, RS.i_]`.
///
/// # Examples
///
/// ```ignore
/// use idenso::gamma0;
///
/// let factor = gamma0!();
/// let default_tensor = gamma0!(i, j);
/// let pattern_tensor = gamma0!(RS.i__, RS.j__);
/// let dimensioned_pattern_tensor = gamma0!([RS.d_, RS.i_], [RS.d_, RS.j_]);
/// ```
#[macro_export]
macro_rules! gamma0 {
    () => {
        symbolica::atom::FunctionBuilder::new($crate::dirac::AGS.gamma0)
            .add_arg(symbolica::atom::Atom::var(spenso::network::tags::SPENSO_TAG.chain_in))
            .add_arg(symbolica::atom::Atom::var(spenso::network::tags::SPENSO_TAG.chain_out))
            .finish()
    };
    ($base:ident . $i:ident, $($rest:tt)*) => {
        $crate::gamma0!(@tensor spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::Bispinor {},
            [$base.$i],
        ); $($rest)*)
    };
    ([$($i:expr),+ $(,)?], $($rest:tt)*) => {
        $crate::gamma0!(@tensor spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::Bispinor {},
            [$($i),+],
        ); $($rest)*)
    };
    ($i:ident, $($rest:tt)*) => {{
        let bis = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::Bispinor {},
            4,
        );
        $crate::gamma0!(@tensor spenso::slot!(bis, $i); $($rest)*)
    }};
    ($i:literal, $($rest:tt)*) => {{
        let bis = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::Bispinor {},
            4,
        );
        $crate::gamma0!(@tensor spenso::slot!(bis, $i); $($rest)*)
    }};
    ($i:expr, $($rest:tt)*) => {
        $crate::gamma0!(@tensor $i; $($rest)*)
    };
    (@tensor $i:expr; $base:ident . $j:ident) => {
        $crate::gamma0!(@tensor_done $i, spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::Bispinor {},
            [$base.$j],
        ))
    };
    (@tensor $i:expr; [$($j:expr),+ $(,)?]) => {
        $crate::gamma0!(@tensor_done $i, spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::Bispinor {},
            [$($j),+],
        ))
    };
    (@tensor $i:expr; $j:ident) => {{
        let bis = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::Bispinor {},
            4,
        );
        $crate::gamma0!(@tensor_done $i, spenso::slot!(bis, $j))
    }};
    (@tensor $i:expr; $j:literal) => {{
        let bis = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::Bispinor {},
            4,
        );
        $crate::gamma0!(@tensor_done $i, spenso::slot!(bis, $j))
    }};
    (@tensor $i:expr; $j:expr) => {
        $crate::gamma0!(@tensor_done $i, $j)
    };
    (@tensor_done $i:expr, $j:expr) => {
        symbolica::atom::FunctionBuilder::new($crate::dirac::AGS.gamma0)
            .add_arg(spenso::symbolica_atom::IntoAtom::into_atom($i))
            .add_arg(spenso::symbolica_atom::IntoAtom::into_atom($j))
            .finish()
    };
}

/// Builds a gamma-five factor or tensor.
///
/// With no arguments, this builds a chain factor using the placeholder indices
/// `in` and `out`; the surrounding chain or trace owns the physical endpoints.
/// With two arguments, this builds the ordinary gamma-five tensor between
/// explicit spinor endpoints.
///
/// Pattern variables such as `RS.i__` are wrapped as bispinor arguments.
#[macro_export]
macro_rules! gamma5 {
    () => {
        symbolica::atom::FunctionBuilder::new($crate::dirac::AGS.gamma5)
            .add_arg(symbolica::atom::Atom::var(spenso::network::tags::SPENSO_TAG.chain_in))
            .add_arg(symbolica::atom::Atom::var(spenso::network::tags::SPENSO_TAG.chain_out))
            .finish()
    };
    ($base:ident . $i:ident, $($rest:tt)*) => {
        $crate::gamma5!(@tensor spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::Bispinor {},
            [$base.$i],
        ); $($rest)*)
    };
    ($i:expr, $($rest:tt)*) => {
        $crate::gamma5!(@tensor $i; $($rest)*)
    };
    (@tensor $i:expr; $base:ident . $j:ident) => {
        $crate::gamma5!(@tensor_done $i, spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::Bispinor {},
            [$base.$j],
        ))
    };
    (@tensor $i:expr; $j:expr) => {
        $crate::gamma5!(@tensor_done $i, $j)
    };
    (@tensor_done $i:expr, $j:expr) => {
        symbolica::atom::FunctionBuilder::new($crate::dirac::AGS.gamma5)
            .add_arg(spenso::symbolica_atom::IntoAtom::into_atom($i))
            .add_arg(spenso::symbolica_atom::IntoAtom::into_atom($j))
            .finish()
    };
}
