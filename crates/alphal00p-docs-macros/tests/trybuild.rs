#[test]
fn documentation_attributes_compile_and_diagnose() {
    let cases = trybuild::TestCases::new();
    cases.pass("tests/ui/pass.rs");
    cases.pass("tests/ui/source-backed-pass.rs");
    cases.compile_fail("tests/ui/wrong-target.rs");
    cases.compile_fail("tests/ui/macro-wrong-target.rs");
    cases.compile_fail("tests/ui/unknown-format.rs");
    cases.compile_fail("tests/ui/owner-mismatch.rs");
    cases.compile_fail("tests/ui/owner-invalid.rs");
    cases.compile_fail("tests/ui/owner-wrong-attribute.rs");
}
