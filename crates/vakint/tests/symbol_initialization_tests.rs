//! Public initialization must not parse a numeric literal as a Symbol.

use symbolica::atom::{Atom, AtomCore};
use symbolica::domains::float::{Complex, Float};
use vakint::{Vakint, VakintSettings, symbols::S, vakint_parse};

#[test]
fn public_initialization_preserves_namespaced_imaginary_placeholder() {
    Vakint::initialize_vakint_symbols();
    Vakint::initialize_vakint_symbols();
    assert_eq!(S.cmplx_i.get_namespace(), "vakint");
    let placeholder = Atom::var(S.cmplx_i);
    assert_eq!(vakint_parse!("vakint::𝑖").unwrap(), placeholder);

    let numeric_i = vakint_parse!("𝑖").unwrap();
    assert_ne!(numeric_i, placeholder);
    assert_eq!(numeric_i.pow(Atom::num(2)), Atom::num(-1));

    // The historical FORM-facing symbolic placeholder still becomes the
    // numerical imaginary unit only at the existing evaluator boundary.
    let settings = VakintSettings::default();
    let (_, complex) =
        Vakint::get_constants_map(&settings, &Default::default(), &Default::default(), None)
            .unwrap();
    let precision = settings.get_binary_precision();
    assert_eq!(
        complex.get(&placeholder),
        Some(&Complex::new(
            Float::with_val(precision, 0),
            Float::with_val(precision, 1),
        )),
    );
}
