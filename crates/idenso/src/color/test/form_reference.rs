use super::*;
use crate::coad;
use insta::assert_snapshot;
use spenso::{antisym, g, mink, shadowing::IntoAtom};
use symbolica::{function, symbol};

macro_rules! fco {
    ($r:ident, $a:tt, $b:tt, $c:tt) => {
        color_f!(
            slot!($r.coad_da, $a),
            slot!($r.coad_da, $b),
            slot!($r.coad_da, $c),
        )
    };
}

macro_rules! tco {
    ($r:ident, $a:tt) => {
        color_t!(slot!($r.coad_da, $a))
    };
}

#[test]
fn one_generator_trace_vanishes() {
    test_initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, i))),
        default_namespace = "spenso"
    );

    // Strict representation dimensions are atomic. Restore the FORM dimension
    // after the typed operation, before comparing with the original reference.
    let dimension = parse_lit!(Nc ^ 2 - 1, default_namespace = "spenso");
    let atomic_dimension = TestReps::new().coad_da.dim.to_symbolic();
    let admitted = expr
        .replace(dimension.to_pattern())
        .with(atomic_dimension.to_pattern());
    assert_eq!(
        admitted
            .replace(atomic_dimension.to_pattern())
            .with(dimension.to_pattern()),
        expr
    );
    let result = crate::tensor::SymbolicTensor::infer(admitted)
        .unwrap()
        .simplify_color(crate::color::ColorSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression()
        .replace(atomic_dimension.to_pattern())
        .with(dimension.to_pattern());

    assert!(result.is_zero());
}

#[test]
fn two_generator_trace_normalizes_to_tr_metric() {
    test_initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j)))
            * t(coad(Nc ^ 2 - 1, b), cof(Nc, j), dind(cof(Nc, i))),
        default_namespace = "spenso"
    );

    // Strict representation dimensions are atomic. Restore the FORM dimension
    // after the typed operation, before comparing with the original reference.
    let dimension = parse_lit!(Nc ^ 2 - 1, default_namespace = "spenso");
    let atomic_dimension = TestReps::new().coad_da.dim.to_symbolic();
    let admitted = expr
        .replace(dimension.to_pattern())
        .with(atomic_dimension.to_pattern());
    assert_eq!(
        admitted
            .replace(atomic_dimension.to_pattern())
            .with(dimension.to_pattern()),
        expr
    );
    let result = crate::tensor::SymbolicTensor::infer(admitted)
        .unwrap()
        .simplify_color(crate::color::ColorSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression()
        .replace(atomic_dimension.to_pattern())
        .with(dimension.to_pattern());

    assert_snapshot!(result.to_bare_ordered_string(), @"g(coad(-1+Nc^2,a),coad(-1+Nc^2,b))*idx(2,cof(Nc))");
}

#[test]
fn fierz_generator_contraction() {
    test_initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j)))
            * t(coad(Nc ^ 2 - 1, a), cof(Nc, k), dind(cof(Nc, l))),
        default_namespace = "spenso"
    );

    // Strict representation dimensions are atomic. Restore the FORM dimension
    // after the typed operation, before comparing with the original reference.
    let dimension = parse_lit!(Nc ^ 2 - 1, default_namespace = "spenso");
    let atomic_dimension = TestReps::new().coad_da.dim.to_symbolic();
    let admitted = expr
        .replace(dimension.to_pattern())
        .with(atomic_dimension.to_pattern());
    assert_eq!(
        admitted
            .replace(atomic_dimension.to_pattern())
            .with(dimension.to_pattern()),
        expr
    );
    let result = crate::tensor::SymbolicTensor::infer(admitted)
        .unwrap()
        .simplify_color(crate::color::ColorSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression()
        .replace(atomic_dimension.to_pattern())
        .with(dimension.to_pattern());

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((result.expand()).as_atom_view().to_owned()).unwrap().contract(crate::tensor::ContractionSettings::default().without_rank_one_tensors()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"-1*Nc^(-1)*g(cof(Nc,i),dind(cof(Nc,j)))*g(cof(Nc,k),dind(cof(Nc,l)))*idx(2,cof(Nc))+g(cof(Nc,i),dind(cof(Nc,l)))*g(cof(Nc,k),dind(cof(Nc,j)))*idx(2,cof(Nc))");
}

#[test]
fn separated_generator_casimir_shortcut() {
    test_initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j)))
            * t(coad(Nc ^ 2 - 1, b), cof(Nc, j), dind(cof(Nc, k)))
            * t(coad(Nc ^ 2 - 1, a), cof(Nc, k), dind(cof(Nc, l))),
        default_namespace = "spenso"
    );

    // Strict representation dimensions are atomic. Restore the FORM dimension
    // after the typed operation, before comparing with the original reference.
    let dimension = parse_lit!(Nc ^ 2 - 1, default_namespace = "spenso");
    let atomic_dimension = TestReps::new().coad_da.dim.to_symbolic();
    let admitted = expr
        .replace(dimension.to_pattern())
        .with(atomic_dimension.to_pattern());
    assert_eq!(
        admitted
            .replace(atomic_dimension.to_pattern())
            .with(dimension.to_pattern()),
        expr
    );
    let result = crate::tensor::SymbolicTensor::infer(admitted)
        .unwrap()
        .simplify_color(crate::color::ColorSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression()
        .replace(atomic_dimension.to_pattern())
        .with(dimension.to_pattern());

    use crate::shorthands::UndoShorthands;

    // Keep the original FORM identity, allowing the equivalent Fierz result
    // -T_R/N. Factor out exactly the same generator before reducing only its
    // scalar coefficient; the tensor numerator is never distributed.
    let expected = parse!(
        "-1/2*cas(2,coad(-1+Nc^2))*t(coad(-1+Nc^2,b),cof(Nc,i),dind(cof(Nc,l)))+cas(2,cof(Nc))*t(coad(-1+Nc^2,b),cof(Nc,i),dind(cof(Nc,l)))",
        default_namespace = "spenso"
    );
    let generator = parse_lit!(
        t(coad(Nc ^ 2 - 1, b), cof(Nc, i), dind(cof(Nc, l))),
        default_namespace = "spenso"
    );
    let coefficient = |expression: Atom| {
        // The graph's slot grammar still uses atomic dimensions when undoing
        // the chain. Restore the original dimension before scalar comparison.
        let factored = expression
            .replace(dimension.to_pattern())
            .with(atomic_dimension.to_pattern())
            .undo_chain::<AbstractIndex>()
            .unwrap()
            .replace(atomic_dimension.to_pattern())
            .with(dimension.to_pattern())
            .collect_factors();
        let scalar = factored
            .replace(generator.to_pattern())
            .with(Atom::num(1).to_pattern());
        assert_eq!(factored, &scalar * &generator);
        scalar.to_cof_dimension_invariants().together().cancel()
    };
    assert_eq!(coefficient(result), coefficient(expected));
}

#[test]
fn chain_one_generator_trace_vanishes() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), i),
        color_t!(slot!(r.coad_da, a)),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"0");
}

#[test]
fn chain_two_generator_trace_normalizes() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), i),
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, b)),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"g(coad(dA,a),coad(dA,b))*idx(2,cof(Nc))");
}

#[test]
fn three_generator_trace_terminal() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), i),
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, b)),
        color_t!(slot!(r.coad_da, c)),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"f(coad(dA,a),coad(dA,b),coad(dA,c))*idx(2,cof(Nc))*𝑖/2+trace(cof(Nc),sym(t(coad(dA,a),in,out),t(coad(dA,b),in,out),t(coad(dA,c),in,out)))");
}

#[test]
fn adjacent_generator_casimir_chain() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), k),
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, a)),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"cas(2,cof(Nc))*g(cof(Nc,i),dind(cof(Nc,k)))");
}

#[test]
fn separated_generator_casimir_trace() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, b)),
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, c)),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"-1/2*cas(2,coad(dA))*g(coad(dA,b),coad(dA,c))*idx(2,cof(Nc))+cas(2,cof(Nc))*g(coad(dA,b),coad(dA,c))*idx(2,cof(Nc))");
}

#[test]
fn two_f_loop_contracts_to_ca_metric() {
    test_initialize();
    let expr = parse_lit!(
        f(coad(dA, a), coad(dA, c), coad(dA, x)) * f(coad(dA, b), coad(dA, c), coad(dA, x)),
        default_namespace = "spenso"
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"cas(2,coad(dA))*g(coad(dA,a),coad(dA,b))");
}

#[test]
fn metric_contraction_respects_antisymmetric_normalization() {
    test_initialize();
    let coad8 = ColorAdjoint {}.new_rep(8);
    // Seed the symbol order from the three-loop graph that exposed the lost
    // permutation parity in the old sequence-wildcard contraction.
    let a = slot!(coad8, standalone_metric_a);
    let x = slot!(coad8, standalone_metric_x);
    let u = slot!(coad8, standalone_metric_u);
    let v = slot!(coad8, standalone_metric_v);
    let _ = crate::tensor::SymbolicTensor::infer(
        (color_f!(a, x, u) * color_f!(a, x, v))
            .as_atom_view()
            .to_owned(),
    )
    .unwrap()
    .contract(crate::tensor::ContractionSettings::default().without_rank_one_tensors())
    .unwrap()
    .resolved()
    .unwrap()
    .into_expression();
    let r = slot!(coad8, standalone_metric_r);
    let s = slot!(coad8, standalone_metric_s);
    let j = slot!(coad8, standalone_metric_j);

    let contracted = crate::tensor::SymbolicTensor::infer(
        (g!(r, s) * color_f!(u, j, r)).as_atom_view().to_owned(),
    )
    .unwrap()
    .contract(crate::tensor::ContractionSettings::default().without_rank_one_tensors())
    .unwrap()
    .resolved()
    .unwrap()
    .into_expression();
    // Re-normalizing f after direct slot substitution may change its displayed
    // argument order and coefficient, so compare with a freshly built target.
    assert_eq!(
        contracted.to_bare_ordered_string(),
        color_f!(u, j, s).to_bare_ordered_string(),
    );

    let closed = g!(u, v) * g!(r, s) * color_f!(u, j, r) * color_f!(v, j, s);
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((closed).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        Atom::num(8) * color_cas!(2, &coad8)
    );
}

#[test]
fn metric_contraction_only_replaces_the_immediate_slot() {
    test_initialize();
    let coad8 = ColorAdjoint {}.new_rep(8);
    let u = slot!(coad8, standalone_nested_metric_u);
    let j = slot!(coad8, standalone_nested_metric_j);
    let r = slot!(coad8, standalone_nested_metric_r);
    let s = slot!(coad8, standalone_nested_metric_s);
    let probe = symbol!("spenso::standalone_metric_probe");
    let structure = color_f!(u, j, r);
    let r_atom = r.into_atom();
    let s_atom = s.into_atom();

    // This raw owner regression intentionally leaves the outer wrapper untagged:
    // its nested tensor-looking metadata must not be treated as a second port.
    let source = g!(r, s) * function!(probe, &structure, r_atom);
    let contracted = crate::shorthands::schoonschip::SlotContraction::run(
        source.as_view(),
        true,
        false,
        &mut Vec::new(),
    );
    assert_eq!(contracted, function!(probe, structure, s_atom));
}

#[test]
fn three_f_loop_contracts_to_ca_f() {
    test_initialize();
    let expr = parse_lit!(
        f(coad(dA, a), coad(dA, b), coad(dA, e))
            * f(coad(dA, b), coad(dA, c), coad(dA, f_))
            * f(coad(dA, c), coad(dA, a), coad(dA, g_)),
        default_namespace = "spenso"
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"1/2*cas(2,coad(dA))*f(coad(dA,e),coad(dA,f_),coad(dA,g_))");
}

#[test]
fn three_f_loop_preserves_antisymmetric_orientation() {
    test_initialize();
    let coad8 = ColorAdjoint {}.new_rep(8);
    // Seed the symbol order that made the unoriented triangle reduction pick
    // the wrong sign in the three-loop graph.
    let r = slot!(coad8, standalone_triangle_r);
    let s = slot!(coad8, standalone_triangle_s);
    let x = slot!(coad8, standalone_triangle_x);
    let z = slot!(coad8, standalone_triangle_z);
    let _ = g!(r, s) * color_f!(x, z, r);
    let u = slot!(coad8, standalone_triangle_u);
    let y = slot!(coad8, standalone_triangle_y);
    let triangle = color_f!(u, x, r) * color_f!(y, z, r) * color_f!(u, z, s);
    let closing_structure = color_f!(s, x, y);
    let expected = Atom::num(-1) * color_cas!(2, &coad8) / Atom::num(2) * color_f!(s, x, y);

    assert_eq!(
        crate::tensor::SymbolicTensor::infer((triangle).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        expected
    );
    // Closing the loop before or after triangle reduction must keep the same
    // antisymmetric orientation.
    assert_eq!(
        crate::tensor::SymbolicTensor::infer(
            (&triangle * &closing_structure).as_atom_view().to_owned()
        )
        .unwrap()
        .simplify_color(crate::color::ColorSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression(),
        crate::tensor::SymbolicTensor::infer(
            (crate::tensor::SymbolicTensor::infer((triangle).as_atom_view().to_owned())
                .unwrap()
                .simplify_color(crate::color::ColorSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression()
                * closing_structure)
                .as_atom_view()
                .to_owned()
        )
        .unwrap()
        .simplify_color(crate::color::ColorSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression(),
    );
}

#[test]
fn six_f_k33_odd_automorphism_simplifies_to_zero() {
    test_initialize();
    let r = TestReps::new();
    let contraction = fco!(r, a, b, c)
        * fco!(r, d, e, f)
        * fco!(r, h, i, j)
        * fco!(r, a, d, h)
        * fco!(r, b, e, i)
        * fco!(r, c, f, j);
    // Exchanging the first two K3,3 vertices is a dummy relabeling, but it
    // reverses three f tensors. Thus C=-C and the contraction is zero.
    let odd_relabeling = fco!(r, d, e, f)
        * fco!(r, a, b, c)
        * fco!(r, h, i, j)
        * fco!(r, d, a, h)
        * fco!(r, e, b, i)
        * fco!(r, f, c, j);

    assert_eq!(odd_relabeling, -&contraction);
    // Color simplification applies signed tensor canonicalization only to its
    // extracted color factor.
    assert!(
        crate::tensor::SymbolicTensor::infer((contraction).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
            .is_zero()
    );
}

#[test]
fn six_f_k33_with_structured_indices_simplifies_to_zero() {
    test_initialize();
    let cooking = crate::CookSettings::reversible()
        .with_index_payload_filter(None)
        .with_output_tags(["idenso::canonical_color_index"]);
    let hedge = symbol!("spenso::hedge");
    let structured_slot = |index| coad!(8, function!(hedge, Atom::num(index as i64)));
    let a = structured_slot(1);
    let b = structured_slot(2);
    let c = structured_slot(3);
    let d = structured_slot(4);
    let e = structured_slot(5);
    let f = structured_slot(6);
    let h = structured_slot(7);
    let i = structured_slot(8);
    let j = structured_slot(9);
    let contraction = color_f!(&a, &b, &c)
        * color_f!(&d, &e, &f)
        * color_f!(&h, &i, &j)
        * color_f!(&a, &d, &h)
        * color_f!(&b, &e, &i)
        * color_f!(&c, &f, &j);

    assert!(
        crate::tensor::SymbolicTensor::infer(cooking.cook_indices(contraction.as_view()))
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
            .is_zero()
    );
    // Pruning a closed zero factor must restore every surviving structured
    // label, including free ports on the remaining color tensor.
    let surviving = color_f!(
        structured_slot(21),
        structured_slot(22),
        structured_slot(23)
    );
    let input = (contraction + Atom::one()) * &surviving;
    let simplified = crate::tensor::SymbolicTensor::infer(cooking.cook_indices(input.as_view()))
        .unwrap()
        .simplify_color(crate::color::ColorSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression();
    let simplified = cooking.uncook(simplified.as_view());
    assert_eq!(simplified, surviving);
    assert!(simplified.contains_symbol(hedge));
}

#[test]
fn signed_pruning_reenters_color_simplification() {
    use spenso::structure::slot::IsAbstractSlot;

    test_initialize();
    let r = TestReps::new();
    let odd = fco!(r, a, b, c)
        * fco!(r, d, e, f)
        * fco!(r, h, i, j)
        * fco!(r, a, d, h)
        * fco!(r, b, e, i)
        * fco!(r, c, f, j);
    let vanishing_factor = fco!(r, u, v, w);
    let surviving_factor = fco!(r, x, y, z);
    // An algebraically vanishing coefficient does not excuse incompatible
    // source ports: validate before signed-zero pruning.
    let incompatible = (&odd * &vanishing_factor + &surviving_factor) * &surviving_factor;
    assert!(crate::tensor::SymbolicTensor::infer(incompatible).is_err());
    let vanishing_factor = vanishing_factor
        .replace(slot!(r.coad_da, u).to_atom().to_pattern())
        .with(slot!(r.coad_da, x).to_atom().to_pattern())
        .replace(slot!(r.coad_da, v).to_atom().to_pattern())
        .with(slot!(r.coad_da, y).to_atom().to_pattern())
        .replace(slot!(r.coad_da, w).to_atom().to_pattern())
        .with(slot!(r.coad_da, z).to_atom().to_pattern());
    assert_eq!(vanishing_factor, surviving_factor);
    let expected = crate::tensor::SymbolicTensor::infer(
        (surviving_factor.clone().pow(2)).as_atom_view().to_owned(),
    )
    .unwrap()
    .simplify_color(crate::color::ColorSimplifySettings::default())
    .unwrap()
    .resolved()
    .unwrap()
    .into_expression();

    assert_eq!(
        crate::tensor::SymbolicTensor::infer(
            ((odd * vanishing_factor + &surviving_factor) * surviving_factor)
                .as_atom_view()
                .to_owned()
        )
        .unwrap()
        .simplify_color(crate::color::ColorSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression(),
        expected
    );
}

#[test]
fn color_simplification_does_not_canonicalize_lorentz_tensors() {
    test_initialize();
    let r = TestReps::new();
    let a = mink!(4, color_only_canonicalization_a);
    let b = mink!(4, color_only_canonicalization_b);
    let c = mink!(4, color_only_canonicalization_c);
    let odd_lorentz =
        antisym!(a.clone(), b.clone()) * antisym!(b.clone(), c.clone()) * antisym!(c, a);
    let color = fco!(r, x, y, z);
    let expected = crate::tensor::SymbolicTensor::infer((color).as_atom_view().to_owned())
        .unwrap()
        .simplify_color(crate::color::ColorSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression()
        * &odd_lorentz;

    assert!(!expected.is_zero());
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((color * &odd_lorentz).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        expected
    );
    assert!(
        odd_lorentz
            .canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .expect("test expression should canonicalize")
            .is_zero()
    );
}

#[test]
fn two_f_even_automorphism_survives_canonicalization() {
    test_initialize();
    let r = TestReps::new();
    let contraction = fco!(r, a, b, c) * fco!(r, a, b, c);

    // Exchanging two contracted edges reverses both f tensors, so the total
    // graph-automorphism sign is positive.
    assert!(
        !contraction
            .canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .expect("test expression should canonicalize")
            .is_zero()
    );
}

#[test]
fn four_f_closed_ghost_topology_preserves_antisymmetric_orientation() {
    test_initialize();
    let expr = parse_lit!(
        -1 / 8
            * g(coad(8, hedge(0)), coad(8, hedge(1)))
            * f(coad(8, hedge(3)), coad(8, hedge(0)), coad(8, hedge(5)))
            * f(coad(8, hedge(3)), coad(8, hedge(7)), coad(8, hedge(9)))
            * f(coad(8, hedge(7)), coad(8, hedge(11)), coad(8, hedge(1)))
            * f(coad(8, hedge(9)), coad(8, hedge(5)), coad(8, hedge(11))),
        default_namespace = "spenso"
    );

    // The FORM fixture labels indices with hedge(n); typed admission uses the
    // existing reversible index encoding and restores it at the Atom boundary.
    let cooking = crate::CookSettings::indices().with_mode(crate::CookMode::ReversibleEncoding);
    let admitted = cooking.try_cook_indices(expr.as_view()).unwrap();
    assert_eq!(cooking.uncook(admitted.as_view()), expr);
    let result = crate::tensor::SymbolicTensor::infer(admitted)
        .unwrap()
        .simplify_color(crate::color::ColorSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression();
    let result = cooking.uncook(result.as_view());

    assert_snapshot!(result.to_bare_ordered_string(), @"(cas(2,coad(8)))^2*1/2");
}

#[test]
fn factored_four_f_closed_ghost_topology_preserves_antisymmetric_orientation() {
    test_initialize();
    let expr = parse_lit!(
        -1 / 8
            * g(coad(8, hedge(0)), coad(8, hedge(1)))
            * (5 / 16
                * x
                * f(coad(8, hedge(3)), coad(8, hedge(0)), coad(8, hedge(5)))
                * f(coad(8, hedge(3)), coad(8, hedge(7)), coad(8, hedge(9)))
                * f(coad(8, hedge(7)), coad(8, hedge(11)), coad(8, hedge(1)))
                * f(coad(8, hedge(9)), coad(8, hedge(5)), coad(8, hedge(11)))
                + 3 / 8
                    * y
                    * f(coad(8, hedge(3)), coad(8, hedge(0)), coad(8, hedge(5)))
                    * f(coad(8, hedge(3)), coad(8, hedge(7)), coad(8, hedge(9)))
                    * f(coad(8, hedge(7)), coad(8, hedge(11)), coad(8, hedge(1)))
                    * f(coad(8, hedge(9)), coad(8, hedge(5)), coad(8, hedge(11)))),
        default_namespace = "spenso"
    );

    // The FORM fixture labels indices with hedge(n); typed admission uses the
    // existing reversible index encoding and restores it at the Atom boundary.
    let cooking = crate::CookSettings::indices().with_mode(crate::CookMode::ReversibleEncoding);
    let admitted = cooking.try_cook_indices(expr.as_view()).unwrap();
    assert_eq!(cooking.uncook(admitted.as_view()), expr);
    let result = crate::tensor::SymbolicTensor::infer(admitted)
        .unwrap()
        .simplify_color(crate::color::ColorSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression();
    let result = cooking.uncook(result.as_view());

    assert_snapshot!(result.expand().to_bare_ordered_string(), @"(cas(2,coad(8)))^2*3/16*y+(cas(2,coad(8)))^2*5/32*x");
}

#[test]
fn mixed_trace_structure_contraction() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, b)),
        color_t!(slot!(r.coad_da, d)),
    ) * color_f!(
        slot!(r.coad_da, a),
        slot!(r.coad_da, b),
        slot!(r.coad_da, c),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"cas(2,coad(dA))*g(coad(dA,c),coad(dA,d))*idx(2,cof(Nc))*𝑖/2");
}

#[test]
fn symmetric_invariant_d33_partial_contraction() {
    test_initialize();
    let expr = parse_lit!(
        d(sym_x, coad(dA, a), coad(dA, b), coad(dA, c))
            * d(sym_y, coad(dA, a), coad(dA, b), coad(dA, d_)),
        default_namespace = "spenso"
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"dA^(-1)*g(coad(dA,c),coad(dA,d_))*gram(3,sym_x,sym_y)");
}

#[test]
fn four_generator_trace_terminal() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, b)),
        color_t!(slot!(r.coad_da, c)),
        color_t!(slot!(r.coad_da, d)),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"-1/6*f(coad(dA,a),coad(dA,c),coad(dA,x))*f(coad(dA,b),coad(dA,d),coad(dA,x))*idx(2,cof(Nc))+1/3*f(coad(dA,a),coad(dA,d),coad(dA,x))*f(coad(dA,b),coad(dA,c),coad(dA,x))*idx(2,cof(Nc))+f(coad(dA,a),coad(dA,b),coad(dA,x))*trace(cof(Nc),sym(t(coad(dA,c),in,out),t(coad(dA,d),in,out),t(coad(dA,x),in,out)))*𝑖/2+f(coad(dA,c),coad(dA,d),coad(dA,x))*trace(cof(Nc),sym(t(coad(dA,a),in,out),t(coad(dA,b),in,out),t(coad(dA,x),in,out)))*𝑖/2+trace(cof(Nc),sym(t(coad(dA,a),in,out),t(coad(dA,b),in,out),t(coad(dA,c),in,out),t(coad(dA,d),in,out)))");
}

// FORM's repository includes a valgrind-oriented size-5 port of color.h's
// tloop.frm check:
// https://github.com/form-dev/form/blob/fdf6af0ce520f7da2dfe0b0b61bf4cef396770be/check/extra/color.frm
//
// These ignored tests keep the upstream expressions and expected invariant
// families close to the implementation. They require higher symmetric
// invariant reductions than the current chain/trace simplifier provides.

fn form_color_tloop_q10(r: &TestReps) -> Atom {
    trace!(
        &r.cof_nc,
        tco!(r, j1),
        tco!(r, j2),
        tco!(r, j3),
        tco!(r, j4),
        tco!(r, j5),
        tco!(r, j1),
        tco!(r, j2),
        tco!(r, j3),
        tco!(r, j4),
        tco!(r, j5),
    )
}

fn form_color_tloop_g10(r: &TestReps) -> Atom {
    fco!(r, i1, i2, j1)
        * fco!(r, i2, i3, j2)
        * fco!(r, i3, i4, j3)
        * fco!(r, i4, i5, j4)
        * fco!(r, i5, i6, j5)
        * fco!(r, i6, i7, j1)
        * fco!(r, i7, i8, j2)
        * fco!(r, i8, i9, j3)
        * fco!(r, i9, i10, j4)
        * fco!(r, i10, i1, j5)
}

fn form_color_tloop_qq5(r: &TestReps) -> Atom {
    trace!(
        &r.cof_nc,
        tco!(r, j1),
        tco!(r, j2),
        tco!(r, j3),
        tco!(r, j4),
        tco!(r, j5),
    ) * trace!(
        &r.cof_nc,
        tco!(r, j1),
        tco!(r, j2),
        tco!(r, j3),
        tco!(r, j4),
        tco!(r, j5),
    )
}

fn form_color_tloop_qg5(r: &TestReps) -> Atom {
    trace!(
        &r.cof_nc,
        tco!(r, j1),
        tco!(r, j2),
        tco!(r, j3),
        tco!(r, j4),
        tco!(r, j5),
    ) * fco!(r, k1, k2, j1)
        * fco!(r, k2, k3, j2)
        * fco!(r, k3, k4, j3)
        * fco!(r, k4, k5, j4)
        * fco!(r, k5, k1, j5)
}

fn form_color_tloop_gg5(r: &TestReps) -> Atom {
    fco!(r, i1, i2, j1)
        * fco!(r, i2, i3, j2)
        * fco!(r, i3, i4, j3)
        * fco!(r, i4, i5, j4)
        * fco!(r, i5, i1, j5)
        * fco!(r, k1, k2, j1)
        * fco!(r, k2, k3, j2)
        * fco!(r, k3, k4, j3)
        * fco!(r, k4, k5, j4)
        * fco!(r, k5, k1, j5)
}

#[test]
#[ignore = "pending FORM color.frm size-5 Q10 invariant-family reduction"]
fn form_github_color_tloop_q10_size_5() {
    test_initialize();
    let r = TestReps::new();
    let expr = form_color_tloop_q10(&r);

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().expand().to_bare_ordered_string(), @"trace(cof(Nc),cyclic(t(coad(dA,j1),in,out),t(coad(dA,j2),in,out),t(coad(dA,j3),in,out),t(coad(dA,j4),in,out),t(coad(dA,j5),in,out),t(coad(dA,j1),in,out),t(coad(dA,j2),in,out),t(coad(dA,j3),in,out),t(coad(dA,j4),in,out),t(coad(dA,j5),in,out)))");
}

#[test]
#[ignore = "pending FORM color.frm size-5 G10 invariant-family reduction"]
fn form_github_color_tloop_g10_size_5() {
    test_initialize();
    let r = TestReps::new();
    let expr = form_color_tloop_g10(&r);

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().expand().to_bare_ordered_string(), @"f(coad(dA,i1),coad(dA,i10),coad(dA,j5))*f(coad(dA,i1),coad(dA,i2),coad(dA,j1))*f(coad(dA,i10),coad(dA,i9),coad(dA,j4))*f(coad(dA,i2),coad(dA,i3),coad(dA,j2))*f(coad(dA,i3),coad(dA,i4),coad(dA,j3))*f(coad(dA,i4),coad(dA,i5),coad(dA,j4))*f(coad(dA,i5),coad(dA,i6),coad(dA,j5))*f(coad(dA,i6),coad(dA,i7),coad(dA,j1))*f(coad(dA,i7),coad(dA,i8),coad(dA,j2))*f(coad(dA,i8),coad(dA,i9),coad(dA,j3))");
}

#[test]
#[ignore = "pending FORM color.frm size-5 QQ5 invariant-family reduction"]
fn form_github_color_tloop_qq5_size_5() {
    test_initialize();
    let r = TestReps::new();
    let expr = form_color_tloop_qq5(&r);

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().expand().to_bare_ordered_string(), @"(trace(cof(Nc),cyclic(t(coad(dA,j1),in,out),t(coad(dA,j2),in,out),t(coad(dA,j3),in,out),t(coad(dA,j4),in,out),t(coad(dA,j5),in,out))))^2");
}

#[test]
#[ignore = "pending FORM color.frm size-5 QG5 invariant-family reduction"]
fn form_github_color_tloop_qg5_size_5() {
    test_initialize();
    let r = TestReps::new();
    let expr = form_color_tloop_qg5(&r);

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().expand().to_bare_ordered_string(), @"-1*f(coad(dA,j1),coad(dA,k1),coad(dA,k2))*f(coad(dA,j2),coad(dA,k2),coad(dA,k3))*f(coad(dA,j3),coad(dA,k3),coad(dA,k4))*f(coad(dA,j4),coad(dA,k4),coad(dA,k5))*f(coad(dA,j5),coad(dA,k1),coad(dA,k5))*trace(cof(Nc),cyclic(t(coad(dA,j1),in,out),t(coad(dA,j2),in,out),t(coad(dA,j3),in,out),t(coad(dA,j4),in,out),t(coad(dA,j5),in,out)))");
}

#[test]
#[ignore = "pending FORM color.frm size-5 GG5 invariant-family reduction"]
fn form_github_color_tloop_gg5_size_5() {
    test_initialize();
    let r = TestReps::new();
    let expr = form_color_tloop_gg5(&r);

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().expand().to_bare_ordered_string(), @"f(coad(dA,i1),coad(dA,i2),coad(dA,j1))*f(coad(dA,i1),coad(dA,i5),coad(dA,j5))*f(coad(dA,i2),coad(dA,i3),coad(dA,j2))*f(coad(dA,i3),coad(dA,i4),coad(dA,j3))*f(coad(dA,i4),coad(dA,i5),coad(dA,j4))*f(coad(dA,j1),coad(dA,k1),coad(dA,k2))*f(coad(dA,j2),coad(dA,k2),coad(dA,k3))*f(coad(dA,j3),coad(dA,k3),coad(dA,k4))*f(coad(dA,j4),coad(dA,k4),coad(dA,k5))*f(coad(dA,j5),coad(dA,k1),coad(dA,k5))");
}

#[test]
#[ignore = "pending color.h simpli strategy over invariant environments"]
fn simpli_contracts_projected_f_pair() {
    test_initialize();
    let expr = parse_lit!(
        f(coad(dA, x), coad(dA, a), coad(dA, c))
            * f(coad(dA, x), coad(dA, b), coad(dA, c))
            * invariant_environment(a, b),
        default_namespace = "spenso"
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"CA*g(coad(dA,a),coad(dA,b))*invariant_environment(a,b)");
}

#[test]
#[ignore = "pending color.tar.gz tloop qloop size 3 family"]
fn tloop_qloop_size_3() {
    test_initialize();
    let result =
        parse_lit!(dA * TR * CF ^ 2 - 3 / 2 * dA * TR * CA * CF + 1 / 2 * dA * TR * CA ^ 2);

    assert_snapshot!(result.to_bare_ordered_string(), @"-3/2*CA*CF*dA*TR+1/2*CA^2*dA*TR+CF^2*dA*TR");
}

#[test]
#[ignore = "pending color.tar.gz tloop gloop size 3 family"]
fn tloop_gloop_size_3() {
    test_initialize();
    let result = parse_lit!(0);

    assert_snapshot!(result.to_bare_ordered_string(), @"0");
}

#[test]
#[ignore = "pending color.tar.gz tloop qqloop size 3 family"]
fn tloop_qqloop_size_3() {
    test_initialize();
    let result = parse_lit!(-1 / 4 * dA * TR ^ 2 * CA + d33(R1, R2));

    assert_snapshot!(result.to_bare_ordered_string(), @"-1/4*CA*dA*TR^2+d33(R1,R2)");
}

#[test]
#[ignore = "pending color.tar.gz tloop qgloop size 3 family"]
fn tloop_qgloop_size_3() {
    test_initialize();
    let result = parse_lit!(1𝑖 / 4 * dA * TR * CA ^ 2);

    assert_snapshot!(result.to_bare_ordered_string(), @"1𝑖/4*CA^2*dA*TR");
}

#[test]
#[ignore = "pending color.tar.gz tloop ggloop size 3 family"]
fn tloop_ggloop_size_3() {
    test_initialize();
    let result = parse_lit!(1 / 4 * dA * CA ^ 3);

    assert_snapshot!(result.to_bare_ordered_string(), @"1/4*CA^3*dA");
}

#[test]
#[ignore = "pending color.tar.gz fixed g14 family"]
fn tloop_g14() {
    test_initialize();
    let result = parse_lit!(
        1 / 648 * dA * CA ^ 7 - 8 / 15 * d444(A1, A2, A3) * CA + 16 / 9 * d644(A1, A2, A3)
    );

    assert_snapshot!(result.to_bare_ordered_string(), @"-8/15*CA*d444(A1,A2,A3)+1/648*CA^7*dA+16/9*d644(A1,A2,A3)");
}

#[test]
#[ignore = "pending color.tar.gz fixed fiveq family"]
fn tloop_fiveq() {
    test_initialize();
    let result = parse_lit!(
        1 / 192 * dA * TR
            ^ 5 * CA
            ^ 3 + 1 / 4 * d33(R1, R2) * TR
            ^ 3 * CA
            ^ 2 + 5 / 48 * d33(R1, R2) * d33(R3, R4) * dA
            ^ -1 * TR * CA + 5 / 48 * d33(R1, R3) * d33(R2, R4) * dA
            ^ -1 * TR * CA + 1 / 8 * d33(R1, R4) * d33(R2, R3) * dA
            ^ -1 * TR * CA + 1 / 16 * d44(R1, A1) * TR
            ^ 4 + 3 / 8 * d433(R3, R1, R2) * TR
            ^ 2 * CA + 1 / 2 * d3333(R1, R2, R3, R4) * TR * CA + d43333a(R5, R2, R1, R4, R3)
    );

    assert_snapshot!(result.to_bare_ordered_string(), @"1/16*TR^4*d44(R1,A1)+1/192*CA^3*dA*TR^5+1/2*CA*TR*d3333(R1,R2,R3,R4)+1/4*CA^2*TR^3*d33(R1,R2)+1/8*CA*dA^(-1)*TR*d33(R1,R4)*d33(R2,R3)+3/8*CA*TR^2*d433(R3,R1,R2)+5/48*CA*dA^(-1)*TR*d33(R1,R2)*d33(R3,R4)+5/48*CA*dA^(-1)*TR*d33(R1,R3)*d33(R2,R4)+d43333a(R5,R2,R1,R4,R3)");
}

#[test]
#[ignore = "pending color.tar.gz su.frm F3F3 case"]
fn su_f3f3() {
    test_initialize();
    let result = parse_lit!(2 * a ^ 3 * NF ^ -1 - 5 / 2 * a ^ 3 * NF + 1 / 2 * a ^ 3 * NF ^ 3);

    assert_snapshot!(result.to_bare_ordered_string(), @"-5/2*NF*a^3+1/2*NF^3*a^3+2*NF^(-1)*a^3");
}

#[test]
#[ignore = "pending color.tar.gz su.frm F4F4 case"]
fn su_f4f4() {
    test_initialize();
    let result = parse_lit!(
        -3 * a
            ^ 4 * nf
            ^ 2 * NF
            ^ -2 + 4 * a
            ^ 4 * nf
            ^ 2 - 7 / 6 * a
            ^ 4 * nf
            ^ 2 * NF
            ^ 2 + 1 / 6 * a
            ^ 4 * nf
            ^ 2 * NF
            ^ 4
    );

    assert_snapshot!(result.to_bare_ordered_string(), @"-3*NF^(-2)*a^4*nf^2+-7/6*NF^2*a^4*nf^2+1/6*NF^4*a^4*nf^2+4*a^4*nf^2");
}

#[test]
#[ignore = "pending color.tar.gz su.frm F4A4 case"]
fn su_f4a4() {
    test_initialize();
    let result = parse_lit!(
        -2 * a ^ 4 * nf * NF + 5 / 3 * a ^ 4 * nf * NF ^ 3 + 1 / 3 * a ^ 4 * nf * NF ^ 5
    );

    assert_snapshot!(result.to_bare_ordered_string(), @"-2*NF*a^4*nf+1/3*NF^5*a^4*nf+5/3*NF^3*a^4*nf");
}

#[test]
#[ignore = "pending color.tar.gz su.frm A4A4 case"]
fn su_a4a4() {
    test_initialize();
    let result =
        parse_lit!(-24 * a ^ 4 * NF ^ 2 + 70 / 3 * a ^ 4 * NF ^ 4 + 2 / 3 * a ^ 4 * NF ^ 6);

    assert_snapshot!(result.to_bare_ordered_string(), @"-24*NF^2*a^4+2/3*NF^6*a^4+70/3*NF^4*a^4");
}

#[test]
#[ignore = "pending color.tar.gz su.frm FnFn n=3 case"]
fn su_fnfn_n3() {
    test_initialize();
    let result = parse_lit!(2 * a ^ 3 * nf ^ 2 * NF ^ -1 - 2 * a ^ 3 * nf ^ 2 * NF);

    assert_snapshot!(result.to_bare_ordered_string(), @"-2*NF*a^3*nf^2+2*NF^(-1)*a^3*nf^2");
}

#[test]
#[ignore = "pending color.tar.gz su.frm FnAn n=3 case"]
fn su_fnan_n3() {
    test_initialize();
    let result = parse_lit!(-1𝑖 * a ^ 3 * nf * NF ^ 2 + 1𝑖 * a ^ 3 * nf * NF ^ 4);

    assert_snapshot!(result.to_bare_ordered_string(), @"-1𝑖*NF^2*a^3*nf+1𝑖*NF^4*a^3*nf");
}

#[test]
#[ignore = "pending color.tar.gz su.frm AnAn n=3 case"]
fn su_anan_n3() {
    test_initialize();
    let result = parse_lit!(-2 * a ^ 3 * NF ^ 3 + 2 * a ^ 3 * NF ^ 5);

    assert_snapshot!(result.to_bare_ordered_string(), @"-2*NF^3*a^3+2*NF^5*a^3");
}

#[test]
#[ignore = "pending color.tar.gz su.frm qloop n=3 case"]
fn su_qloop_n3() {
    test_initialize();
    let result = parse_lit!(-a ^ 3 * nf * NF ^ -2 + a ^ 3 * nf * NF ^ 2);

    assert_snapshot!(result.to_bare_ordered_string(), @"-1*NF^(-2)*a^3*nf+NF^2*a^3*nf");
}

#[test]
#[ignore = "pending color.tar.gz su.frm gloop n=3 case"]
fn su_gloop_n3() {
    test_initialize();
    let result = parse_lit!(0);

    assert_snapshot!(result.to_bare_ordered_string(), @"0");
}
