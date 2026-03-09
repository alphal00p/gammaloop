mod test_utils;
use test_utils::{compare_output, get_vakint};
use vakint::VakintSettings;
use vakint::vakint_parse;

#[test_log::test]
fn test_reduction_1l_a() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    });

    let integral = vakint
        .to_canonical(
            vakint_parse!(
                "(k(1,1)*k(1,2)+k(1,3)*p(1,3))*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
            )"
            )
            .unwrap()
            .as_view(),
            true,
        )
        .unwrap();

    _ = compare_output(
        vakint
            .tensor_reduce(integral.as_view())
            .as_ref()
            .map(|a| a.as_view()),
        vakint_parse!(
            "(\
                -(2*ε-4)^-1*dot(k(1),k(1))*g(1,2)\
            )*topo(I1L(muvsq,1))"
        )
        .unwrap(),
    );
}

#[test_log::test]
fn test_reduction_1l_b() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    });

    let integral = vakint
        .to_canonical(
            vakint_parse!(
                "((k(1,1)*k(1,2))^2*g(1,2)+k(1,3)*p(1,3)+k(1,1)*k(1,2)*p(2,1)*p(3,2))*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
            )"
            )
            .unwrap()
            .as_view(),
            true,
        )
        .unwrap();

    _ = compare_output(
        vakint
            .tensor_reduce(integral.as_view())
            .as_ref()
            .map(|a| a.as_view()),
        vakint_parse!(
            "(\
                dot(k(1),k(1))^2*g(1,2)-(2*ε-4)^-1*dot(p(2),p(3))*dot(k(1),k(1))\
            )*topo(I1L(muvsq,1))"
        )
        .unwrap(),
    );
}

#[test_log::test]
fn test_reduction_2l_a() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    });

    let integral = vakint
        .to_canonical(
            vakint_parse!(
                "(\
                    (k(1,1)*k(2,2))^2*g(1,2)+k(2,3)*p(1,3)+k(1,1)*k(2,2)*p(2,1)*p(3,2)\
                )*topo(I2L(mUVsq,1,2,1))"
            )
            .unwrap()
            .as_view(),
            true,
        )
        .unwrap();

    _ = compare_output(
        vakint
            .tensor_reduce(integral.as_view())
            .as_ref()
            .map(|a| a.as_view()),
        vakint_parse!(
            "(\
                -(2*ε-4)^-1*dot(p(2),p(3))*dot(k(1),k(2))+dot(k(1),k(1))*dot(k(2),k(2))*g(1,2)\
            )*topo(I2L(mUVsq,1,2,1))"
        )
        .unwrap(),
    );
}

#[allow(dead_code)]
fn run_tensor_reduction_tests() {
    test_reduction_1l_a();
    test_reduction_1l_b();
    test_reduction_2l_a();
}
