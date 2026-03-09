use spenso::structure::representation::{LibraryRep, RepName};
use spenso_macros::SimpleRepresentation;

#[rustfmt::skip]
#[derive(
    SimpleRepresentation)]
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Hash,
    PartialOrd,
    Ord,
    Default,
)]
#[representation(name = "spf", dual_name = "SpinAntiFundamental")] // Specify the dual name
pub struct SpinFundamental {}
// This will now generate:
// #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
// pub struct SpinAntiFundamental {} // <-- Custom name used!
// and all impls referencing SpinAntiFundamental instead of DualSpinFundamental

// // --- Usage with Default Dual Name (still works) ---

#[rustfmt::skip]
#[derive(SimpleRepresentation)]
#[representation(name = "vec")] // No dual_name specified
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub struct Vector {}
// // This will generate:
// // #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
// // pub struct DualVector {} // <-- Default name used
// // and impls for Vector and DualVector

// // --- Self-Dual Usage (unchanged) ---
#[derive(
    SimpleRepresentation, Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default,
)]
#[representation(name = "sc", self_dual)]
pub struct Scalar {}

#[test]
fn reps() {
    let sf = SpinFundamental {};

    // let sf_dual = sf.dual().dual();
    // Use the custom type name here!
    let sf_dual: SpinAntiFundamental = sf.dual();
    let v = Vector {};
    let v_dual: DualVector = v.dual(); // Default name still works

    assert_eq!(format!("{}", sf), "spf🠑");
    assert_eq!(format!("{}", sf_dual), "spf🠓"); // Display uses Rep, name is unchanged
    assert_eq!(format!("{}", v), "vec🠑");
    assert_eq!(format!("{}", v_dual), "vec🠓");

    // Check RepName::Dual associated type (implicitly checked by sf.dual() type)
    assert!(
        std::any::TypeId::of::<<SpinFundamental as RepName>::Dual>()
            == std::any::TypeId::of::<SpinAntiFundamental>()
    );
    assert!(
        std::any::TypeId::of::<<Vector as RepName>::Dual>() == std::any::TypeId::of::<DualVector>()
    );

    assert_eq!(sf.dual(), sf_dual);
    assert_eq!(sf_dual.dual(), sf); // Dual of custom dual is base
    assert_eq!(v.dual(), v_dual);
    assert_eq!(v_dual.dual(), v);

    // Check TryFrom<Rep> uses the correct type
    let rep_sf: LibraryRep = sf.into();
    let rep_sf_dual: LibraryRep = sf_dual.into();

    assert_eq!(SpinFundamental::try_from(rep_sf).unwrap(), sf);
    // Use the custom type for TryFrom
    assert_eq!(SpinAntiFundamental::try_from(rep_sf_dual).unwrap(), sf_dual);
    assert!(SpinAntiFundamental::try_from(rep_sf).is_err());

    // Check base() method on custom dual type
    assert_eq!(sf_dual.base(), sf);

    let sc = Scalar::default();
    let sc_dual = sc.dual();
    assert_eq!(sc_dual.base(), sc);

    assert_eq!(format!("{}", sc), "sc");

    // assert_eq!(format!("{}", sf.dual().dual().to_symbolic([])), "spf🠑");
}
