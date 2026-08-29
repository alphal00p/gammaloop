use vakint::{
    RustRedOptions, TensorReductionError, TensorReductionMode, Vakint, VakintSettings, vakint_parse,
};

#[test]
fn tensor_reducer_defaults_to_the_form_backend() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings::default();
    let reducer = vakint.tensor_reducer(&settings);

    assert_eq!(reducer.selected_mode(), &TensorReductionMode::Form);
    assert_eq!(
        reducer
            .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
            .selected_mode(),
        &TensorReductionMode::RustRed(RustRedOptions::new())
    );
}

#[test]
fn default_builder_executes_the_legacy_form_backend() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings::default();
    let input = vakint_parse!("k(1,1)*k(1,2)*topo(I1L(muvsq,1))").unwrap();

    let legacy = vakint.tensor_reduce(&settings, input.as_view()).unwrap();
    let through_default = vakint
        .tensor_reducer(&settings)
        .reduce(input.as_view())
        .unwrap();

    assert_eq!(through_default, legacy);
}

#[test]
fn rustred_mode_stops_before_form_until_the_core_service_exists() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        ..VakintSettings::default()
    };
    let input = vakint_parse!("k(1,mu)*topo(I1L(muvsq,1))").unwrap();

    let error = vakint
        .tensor_reducer(&settings)
        .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
        .reduce(input.as_view())
        .unwrap_err();

    assert!(matches!(error, TensorReductionError::RustRedUnavailable));
}
