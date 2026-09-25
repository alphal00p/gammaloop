#[test]
fn documentation_attributes_compile_and_diagnose() {
    #[alphal00p_docs::ty]
    enum Mode {
        /// The production mode remains documented.
        Production,
        #[cfg(test)]
        /// A probe retained only in unit-test builds.
        TestProbe,
    }
    let mode = __alphal00p_docs_ty_Mode();
    assert_eq!(
        mode.members
            .iter()
            .map(|member| member.name.as_str())
            .collect::<Vec<_>>(),
        ["Production"]
    );
    let _ = (Mode::Production, Mode::TestProbe);

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
