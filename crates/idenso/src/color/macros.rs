/// Builds a color-generator atom.
///
/// With one argument, this builds a fundamental generator factor for use inside
/// `chain!` or `trace!`. The factor uses the chain placeholder indices `in` and
/// `out`; the surrounding chain or trace owns the physical endpoints.
///
/// With three arguments, this builds an explicit generator `t(a,i,j)`. The first
/// argument is an adjoint color index, the second is a fundamental color index,
/// and the third is an anti-fundamental color index.
///
/// Arguments can be typed slots, atoms, atom views, default-dimension
/// identifiers/literals, pattern variables such as `RS.a__`, or bracketed
/// representation argument lists such as `[CS.adj_, RS.a_]`.
///
/// # Examples
///
/// ```ignore
/// use idenso::color_t;
/// use spenso::{slot, trace};
///
/// let factor = color_t!(slot!(coad_na, a));
/// let expr = trace!(&cof_nc, factor);
/// let default_factor = color_t!(a);
/// let default_tensor = color_t!(a, i, j);
/// let pattern_tensor = color_t!(RS.a__, RS.i__, RS.j__);
/// let dimensioned_pattern = color_t!([CS.adj_, RS.a_], [CS.nc_, RS.i_], [CS.nc_, RS.j_]);
/// ```
#[macro_export]
macro_rules! color_t {
    ([$($a:expr),+ $(,)?], $($rest:tt)+) => {
        $crate::color_t!(@tensor spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::ColorAdjoint {},
            [$($a),+],
        ); $($rest)*)
    };
    ($base:ident . $a:ident, $($rest:tt)+) => {
        $crate::color_t!(@tensor spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::ColorAdjoint {},
            [$base.$a],
        ); $($rest)*)
    };
    ($a:ident, $($rest:tt)+) => {{
        let coad = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::ColorAdjoint {},
            $crate::color::CS.na,
        );
        $crate::color_t!(@tensor spenso::slot!(coad, $a); $($rest)*)
    }};
    ($a:literal, $($rest:tt)+) => {{
        let coad = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::ColorAdjoint {},
            $crate::color::CS.na,
        );
        $crate::color_t!(@tensor spenso::slot!(coad, $a); $($rest)*)
    }};
    ($a:expr, $($rest:tt)+) => {
        $crate::color_t!(@tensor $a; $($rest)*)
    };
    ([$($a:expr),+ $(,)?] $(,)?) => {
        $crate::color::CS.chain_t(spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::ColorAdjoint {},
            [$($a),+],
        ))
    };
    ($base:ident . $a:ident $(,)?) => {
        $crate::color::CS.chain_t(spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::ColorAdjoint {},
            [$base.$a],
        ))
    };
    ($a:ident $(,)?) => {{
        let coad = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::ColorAdjoint {},
            $crate::color::CS.na,
        );
        $crate::color::CS.chain_t(spenso::shadowing::IntoAtom::into_atom(spenso::slot!(coad, $a)))
    }};
    ($a:literal $(,)?) => {{
        let coad = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::ColorAdjoint {},
            $crate::color::CS.na,
        );
        $crate::color::CS.chain_t(spenso::shadowing::IntoAtom::into_atom(spenso::slot!(coad, $a)))
    }};
    ($a:expr $(,)?) => {
        $crate::color::CS.chain_t(spenso::shadowing::IntoAtom::into_atom($a))
    };
    (@tensor $a:expr; [$($i:expr),+ $(,)?], $($rest:tt)*) => {
        $crate::color_t!(@tensor2 $a, spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::ColorFundamental {},
            [$($i),+],
        ); $($rest)*)
    };
    (@tensor $a:expr; $base:ident . $i:ident, $($rest:tt)*) => {
        $crate::color_t!(@tensor2 $a, spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::ColorFundamental {},
            [$base.$i],
        ); $($rest)*)
    };
    (@tensor $a:expr; $i:ident, $($rest:tt)*) => {{
        let cof = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::ColorFundamental {},
            $crate::color::CS.nc,
        );
        $crate::color_t!(@tensor2 $a, spenso::slot!(cof, $i); $($rest)*)
    }};
    (@tensor $a:expr; $i:literal, $($rest:tt)*) => {{
        let cof = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::ColorFundamental {},
            $crate::color::CS.nc,
        );
        $crate::color_t!(@tensor2 $a, spenso::slot!(cof, $i); $($rest)*)
    }};
    (@tensor $a:expr; $i:expr, $($rest:tt)*) => {
        $crate::color_t!(@tensor2 $a, $i; $($rest)*)
    };
    (@tensor2 $a:expr, $i:expr; [$($j:expr),+ $(,)?] $(,)?) => {
        $crate::color_t!(@tensor_done $a, $i, spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::ColorAntiFundamental {},
            [$($j),+],
        ))
    };
    (@tensor2 $a:expr, $i:expr; $base:ident . $j:ident $(,)?) => {
        $crate::color_t!(@tensor_done $a, $i, spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::ColorAntiFundamental {},
            [$base.$j],
        ))
    };
    (@tensor2 $a:expr, $i:expr; $j:ident $(,)?) => {{
        let coaf = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::ColorAntiFundamental {},
            $crate::color::CS.nc,
        );
        $crate::color_t!(@tensor_done $a, $i, spenso::slot!(coaf, $j))
    }};
    (@tensor2 $a:expr, $i:expr; $j:literal $(,)?) => {{
        let coaf = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::ColorAntiFundamental {},
            $crate::color::CS.nc,
        );
        $crate::color_t!(@tensor_done $a, $i, spenso::slot!(coaf, $j))
    }};
    (@tensor2 $a:expr, $i:expr; $j:expr $(,)?) => {
        $crate::color_t!(@tensor_done $a, $i, $j)
    };
    (@tensor_done $a:expr, $i:expr, $j:expr) => {
        $crate::color::CS.explicit_t(
            spenso::shadowing::IntoAtom::into_atom($a),
            spenso::shadowing::IntoAtom::into_atom($i),
            spenso::shadowing::IntoAtom::into_atom($j),
        )
    };
}

/// Builds an adjoint color structure-constant atom `f(a,b,c)`.
///
/// Arguments can be typed slots, atoms, atom views, default-dimension adjoint
/// identifiers/literals, pattern variables such as `RS.a__`, or bracketed
/// adjoint representation argument lists such as `[CS.adj_, RS.a_]`.
///
/// # Examples
///
/// ```ignore
/// use idenso::color_f;
/// use spenso::slot;
///
/// let expr = color_f!(slot!(coad_na, a), slot!(coad_na, b), slot!(coad_na, c));
/// let default_tensor = color_f!(a, b, c);
/// let pattern_tensor = color_f!(RS.a__, RS.b__, RS.c__);
/// let dimensioned_pattern = color_f!([CS.adj_, RS.a_], [CS.adj_, RS.b_], [CS.adj_, RS.c_]);
/// ```
#[macro_export]
macro_rules! color_f {
    ([$($a:expr),+ $(,)?], $($rest:tt)+) => {
        $crate::color_f!(@tensor spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::ColorAdjoint {},
            [$($a),+],
        ); $($rest)*)
    };
    ($base:ident . $a:ident, $($rest:tt)+) => {
        $crate::color_f!(@tensor spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::ColorAdjoint {},
            [$base.$a],
        ); $($rest)*)
    };
    ($a:ident, $($rest:tt)+) => {{
        let coad = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::ColorAdjoint {},
            $crate::color::CS.na,
        );
        $crate::color_f!(@tensor spenso::slot!(coad, $a); $($rest)*)
    }};
    ($a:literal, $($rest:tt)+) => {{
        let coad = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::ColorAdjoint {},
            $crate::color::CS.na,
        );
        $crate::color_f!(@tensor spenso::slot!(coad, $a); $($rest)*)
    }};
    ($a:expr, $($rest:tt)+) => {
        $crate::color_f!(@tensor $a; $($rest)*)
    };
    (@tensor $a:expr; [$($b:expr),+ $(,)?], $($rest:tt)*) => {
        $crate::color_f!(@tensor2 $a, spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::ColorAdjoint {},
            [$($b),+],
        ); $($rest)*)
    };
    (@tensor $a:expr; $base:ident . $b:ident, $($rest:tt)*) => {
        $crate::color_f!(@tensor2 $a, spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::ColorAdjoint {},
            [$base.$b],
        ); $($rest)*)
    };
    (@tensor $a:expr; $b:ident, $($rest:tt)*) => {{
        let coad = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::ColorAdjoint {},
            $crate::color::CS.na,
        );
        $crate::color_f!(@tensor2 $a, spenso::slot!(coad, $b); $($rest)*)
    }};
    (@tensor $a:expr; $b:literal, $($rest:tt)*) => {{
        let coad = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::ColorAdjoint {},
            $crate::color::CS.na,
        );
        $crate::color_f!(@tensor2 $a, spenso::slot!(coad, $b); $($rest)*)
    }};
    (@tensor $a:expr; $b:expr, $($rest:tt)*) => {
        $crate::color_f!(@tensor2 $a, $b; $($rest)*)
    };
    (@tensor2 $a:expr, $b:expr; [$($c:expr),+ $(,)?] $(,)?) => {
        $crate::color_f!(@tensor_done $a, $b, spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::ColorAdjoint {},
            [$($c),+],
        ))
    };
    (@tensor2 $a:expr, $b:expr; $base:ident . $c:ident $(,)?) => {
        $crate::color_f!(@tensor_done $a, $b, spenso::structure::representation::RepName::to_symbolic(
            &$crate::representations::ColorAdjoint {},
            [$base.$c],
        ))
    };
    (@tensor2 $a:expr, $b:expr; $c:ident $(,)?) => {{
        let coad = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::ColorAdjoint {},
            $crate::color::CS.na,
        );
        $crate::color_f!(@tensor_done $a, $b, spenso::slot!(coad, $c))
    }};
    (@tensor2 $a:expr, $b:expr; $c:literal $(,)?) => {{
        let coad = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::ColorAdjoint {},
            $crate::color::CS.na,
        );
        $crate::color_f!(@tensor_done $a, $b, spenso::slot!(coad, $c))
    }};
    (@tensor2 $a:expr, $b:expr; $c:expr $(,)?) => {
        $crate::color_f!(@tensor_done $a, $b, $c)
    };
    (@tensor_done $a:expr, $b:expr, $c:expr) => {
        $crate::color::CS.structure_f(
            spenso::shadowing::IntoAtom::into_atom($a),
            spenso::shadowing::IntoAtom::into_atom($b),
            spenso::shadowing::IntoAtom::into_atom($c),
        )
    };
}

/// Alias for [`color_t!`].
#[macro_export]
macro_rules! t {
    ($($args:tt)*) => {
        $crate::color_t!($($args)*)
    };
}

/// Alias for [`color_f!`].
#[macro_export]
macro_rules! f {
    ($($args:tt)*) => {
        $crate::color_f!($($args)*)
    };
}

/// Builds a symmetric color invariant `d(rep,args...)`.
#[macro_export]
macro_rules! color_d {
    ($rep:expr $(, $arg:expr)+ $(,)?) => {
        $crate::color::CS.symmetric_d(
            spenso::shadowing::IntoAtom::into_atom($rep),
            vec![$(spenso::shadowing::IntoAtom::into_atom($arg)),+],
        )
    };
}

/// Builds a contracted symmetric color invariant `d33(left,right)`.
#[macro_export]
macro_rules! color_d33 {
    ($left:expr, $right:expr $(,)?) => {
        $crate::color::CS.d33(
            spenso::shadowing::IntoAtom::into_atom($left),
            spenso::shadowing::IntoAtom::into_atom($right),
        )
    };
}
