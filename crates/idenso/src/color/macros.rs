/// Builds a color-generator atom.
///
/// The one-argument form builds a fundamental generator factor for use inside
/// `chain!` or `trace!`. The adjoint-index argument is converted through
/// `spenso::symbolica_atom::IntoAtom`, so it can be a typed color-adjoint slot,
/// an atom, or an atom view. The factor uses the chain placeholder indices `in`
/// and `out`; the surrounding chain or trace owns the physical endpoints.
///
/// The three-argument form builds an explicit generator `t(a,i,j)`. The
/// `pattern:` form builds the legacy typed raw tensor pattern with explicit
/// symbolic fundamental and adjoint dimensions.
///
/// # Examples
///
/// ```ignore
/// use idenso::color_t;
/// use spenso::{slot, trace};
///
/// let factor = color_t!(slot!(coad_na, a));
/// let expr = trace!(&cof_nc, factor);
/// let raw = color_t!(a, i, j);
/// let typed_pattern = color_t!(pattern: nc_, adj_; a_, i_, j_);
/// ```
#[macro_export]
macro_rules! color_t {
    (pattern: $fundamental_dimension:expr, $adjoint_dimension:expr; $a:expr, $i:expr, $j:expr $(,)?) => {
        $crate::color::CS.t_pattern($fundamental_dimension, $adjoint_dimension, $a, $i, $j)
    };
    ($a:expr $(,)?) => {
        $crate::color::CS.chain_t(spenso::symbolica_atom::IntoAtom::into_atom($a))
    };
    ($a:expr, $i:expr, $j:expr $(,)?) => {
        $crate::color::CS.explicit_t(
            spenso::symbolica_atom::IntoAtom::into_atom($a),
            spenso::symbolica_atom::IntoAtom::into_atom($i),
            spenso::symbolica_atom::IntoAtom::into_atom($j),
        )
    };
}

/// Builds an adjoint color structure-constant atom `f(a,b,c)`.
///
/// Arguments are converted through `spenso::symbolica_atom::IntoAtom`, so typed
/// color-adjoint slots can be used directly. This avoids accidentally mixing
/// parsed symbols that print alike but are not identical atoms.
///
/// # Examples
///
/// ```ignore
/// use idenso::color_f;
/// use spenso::slot;
///
/// let expr = color_f!(slot!(coad_na, a), slot!(coad_na, b), slot!(coad_na, c));
/// ```
#[macro_export]
macro_rules! color_f {
    (pattern: $adjoint_dimension:expr; $a:expr, $b:expr, $c:expr $(,)?) => {
        $crate::color::CS.f_pattern($adjoint_dimension, $a, $b, $c)
    };
    ($a:expr, $b:expr, $c:expr $(,)?) => {
        $crate::color::CS.structure_f(
            spenso::symbolica_atom::IntoAtom::into_atom($a),
            spenso::symbolica_atom::IntoAtom::into_atom($b),
            spenso::symbolica_atom::IntoAtom::into_atom($c),
        )
    };
}

/// Alias for [`color_t!`].
#[macro_export]
macro_rules! t {
    (pattern: $fundamental_dimension:expr, $adjoint_dimension:expr; $a:expr, $i:expr, $j:expr $(,)?) => {
        $crate::color_t!(pattern: $fundamental_dimension, $adjoint_dimension; $a, $i, $j)
    };
    ($a:expr $(,)?) => {
        $crate::color_t!($a)
    };
    ($a:expr, $i:expr, $j:expr $(,)?) => {
        $crate::color_t!($a, $i, $j)
    };
}

/// Alias for [`color_f!`].
#[macro_export]
macro_rules! f {
    (pattern: $adjoint_dimension:expr; $a:expr, $b:expr, $c:expr $(,)?) => {
        $crate::color_f!(pattern: $adjoint_dimension; $a, $b, $c)
    };
    ($a:expr, $b:expr, $c:expr $(,)?) => {
        $crate::color_f!($a, $b, $c)
    };
}

/// Builds a symmetric color invariant `d(rep,args...)`.
#[macro_export]
macro_rules! color_d {
    ($rep:expr $(, $arg:expr)+ $(,)?) => {
        $crate::color::CS.symmetric_d(
            spenso::symbolica_atom::IntoAtom::into_atom($rep),
            vec![$(spenso::symbolica_atom::IntoAtom::into_atom($arg)),+],
        )
    };
}

/// Builds a contracted symmetric color invariant `d33(left,right)`.
#[macro_export]
macro_rules! color_d33 {
    ($left:expr, $right:expr $(,)?) => {
        $crate::color::CS.d33(
            spenso::symbolica_atom::IntoAtom::into_atom($left),
            spenso::symbolica_atom::IntoAtom::into_atom($right),
        )
    };
}
