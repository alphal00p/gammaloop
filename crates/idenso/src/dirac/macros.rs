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
/// let mixed_tensor = gamma!(mu, slot!(bis_d, a), 1);
/// let explicit_tensor = gamma!(slot!(mink_d, mu), slot!(bis_d, a), slot!(bis_d, b));
/// ```
#[macro_export]
macro_rules! gamma {
    ($mu:expr) => {
        $crate::dirac::AGS.chain_gamma(spenso::symbolica_atom::IntoAtom::into_atom($mu))
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

/// Builds a gamma-five factor for use inside `chain!` or `trace!`.
///
/// The factor uses the chain placeholder indices `in` and `out`; the surrounding
/// chain or trace owns the physical endpoints.
#[macro_export]
macro_rules! gamma5 {
    () => {
        $crate::dirac::AGS.chain_gamma5()
    };
}
