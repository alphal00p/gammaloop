mod test_utils;
use symbolica::atom::Atom;
use test_utils::{TestVakint, compare_output, get_vakint};
use vakint::vakint_parse;
use vakint::{RustRedOptions, TensorReductionMode, VakintSettings};

fn tensor_oracle_settings() -> VakintSettings {
    VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    }
}

fn rustred_tensor_oracle_settings() -> VakintSettings {
    VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        ..tensor_oracle_settings()
    }
}

fn reduction_1l_a_input() -> Atom {
    vakint_parse!(
        "(k(1,1)*k(1,2)+k(1,3)*p(1,3))*topo(\
            prop(1,edge(1,1),k(1),muvsq,1)\
        )"
    )
    .unwrap()
}

fn reduction_1l_a_expected() -> Atom {
    vakint_parse!(
        "(\
            -(2*ε-4)^-1*dot(k(1),k(1))*g(1,2)\
        )*topo(I1L(muvsq,1))"
    )
    .unwrap()
}

fn reduction_1l_b_input() -> Atom {
    vakint_parse!(
        "((k(1,1)*k(1,2))^2*g(1,2)+k(1,3)*p(1,3)+k(1,1)*k(1,2)*p(2,1)*p(3,2))*topo(\
            prop(1,edge(1,1),k(1),muvsq,1)\
        )"
    )
    .unwrap()
}

fn reduction_1l_b_expected() -> Atom {
    vakint_parse!(
        "(\
            dot(k(1),k(1))^2*g(1,2)-(2*ε-4)^-1*dot(p(2),p(3))*dot(k(1),k(1))\
        )*topo(I1L(muvsq,1))"
    )
    .unwrap()
}

fn reduction_2l_a_input() -> Atom {
    vakint_parse!(
        "(\
            (k(1,1)*k(2,2))^2*g(1,2)+k(2,3)*p(1,3)+k(1,1)*k(2,2)*p(2,1)*p(3,2)\
        )*topo(I2L(mUVsq,1,2,1))"
    )
    .unwrap()
}

fn reduction_2l_a_expected() -> Atom {
    vakint_parse!(
        "(\
            -(2*ε-4)^-1*dot(p(2),p(3))*dot(k(1),k(2))+dot(k(1),k(1))*dot(k(2),k(2))*g(1,2)\
        )*topo(I2L(mUVsq,1,2,1))"
    )
    .unwrap()
}

fn rustred_reduce_to_short_form(vakint: &TestVakint, input: Atom) -> Atom {
    let reduced = vakint
        .vakint
        .tensor_reducer(&vakint.settings)
        .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
        .reduce(input.as_view())
        .unwrap();
    vakint.to_canonical(reduced.as_view(), true).unwrap()
}

#[test]
fn default_tensor_builder_executes_the_legacy_form_backend() {
    let vakint = get_vakint(VakintSettings::default());
    let input = vakint_parse!("k(1,1)*k(1,2)*topo(I1L(muvsq,1))").unwrap();

    let legacy = vakint.tensor_reduce(input.as_view()).unwrap();
    let through_default = vakint
        .vakint
        .tensor_reducer(&vakint.settings)
        .reduce(input.as_view())
        .unwrap();

    assert_eq!(through_default, legacy);
}

#[test_log::test]
fn test_reduction_1l_a() {
    let vakint = get_vakint(tensor_oracle_settings());

    let integral = vakint
        .to_canonical(reduction_1l_a_input().as_view(), true)
        .unwrap();

    _ = compare_output(
        vakint
            .tensor_reduce(integral.as_view())
            .as_ref()
            .map(|a| a.as_view()),
        reduction_1l_a_expected(),
    );
}

#[test_log::test]
fn test_reduction_1l_b() {
    let vakint = get_vakint(tensor_oracle_settings());

    let integral = vakint
        .to_canonical(reduction_1l_b_input().as_view(), true)
        .unwrap();

    _ = compare_output(
        vakint
            .tensor_reduce(integral.as_view())
            .as_ref()
            .map(|a| a.as_view()),
        reduction_1l_b_expected(),
    );
}

#[test_log::test]
fn test_reduction_2l_a() {
    let vakint = get_vakint(tensor_oracle_settings());

    let integral = vakint
        .to_canonical(reduction_2l_a_input().as_view(), true)
        .unwrap();

    _ = compare_output(
        vakint
            .tensor_reduce(integral.as_view())
            .as_ref()
            .map(|a| a.as_view()),
        reduction_2l_a_expected(),
    );
}

#[test]
fn test_reduction_1l_a_rustred() {
    let vakint = get_vakint(rustred_tensor_oracle_settings());
    assert_eq!(
        rustred_reduce_to_short_form(&vakint, reduction_1l_a_input()),
        reduction_1l_a_expected()
    );
}

#[test]
fn test_reduction_2l_a_rustred() {
    let vakint = get_vakint(rustred_tensor_oracle_settings());
    let explicit_input = vakint
        .to_canonical(reduction_2l_a_input().as_view(), false)
        .unwrap();
    assert_eq!(
        rustred_reduce_to_short_form(&vakint, explicit_input),
        reduction_2l_a_expected()
    );
}

#[allow(dead_code)]
fn run_tensor_reduction_tests() {
    test_reduction_1l_a();
    test_reduction_1l_b();
    test_reduction_2l_a();
}
