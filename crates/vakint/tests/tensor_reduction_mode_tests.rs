mod test_utils;

use std::{
    hint::black_box,
    time::{Duration, Instant},
};

use spenso::vector_symbol;
use symbolica::atom::{Atom, AtomCore, AtomView};
use test_utils::{TestVakint, get_vakint};
use vakint::{TensorReductionMethod, Vakint, VakintError, VakintSettings, vakint_parse};

#[derive(Clone, Copy)]
struct ProjectedFixture {
    name: &'static str,
    loop_ids: &'static [usize],
    projector_ids: &'static [usize],
}

const PARITY_FIXTURES: &[ProjectedFixture] = &[
    ProjectedFixture {
        name: "one_loop_rank_two",
        loop_ids: &[1, 1],
        projector_ids: &[1, 1],
    },
    ProjectedFixture {
        name: "one_loop_rank_six",
        loop_ids: &[1, 1, 1, 1, 1, 1],
        projector_ids: &[1, 1, 2, 2, 3, 3],
    },
    ProjectedFixture {
        name: "one_loop_rank_ten",
        loop_ids: &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        projector_ids: &[1, 1, 1, 1, 2, 2, 2, 2, 3, 3],
    },
    ProjectedFixture {
        name: "two_loop_rank_four",
        loop_ids: &[1, 1, 2, 2],
        projector_ids: &[1, 1, 2, 2],
    },
    ProjectedFixture {
        name: "two_loop_rank_eight",
        loop_ids: &[1, 1, 1, 1, 2, 2, 2, 2],
        projector_ids: &[1, 1, 1, 1, 2, 2, 3, 3],
    },
    ProjectedFixture {
        name: "two_loop_rank_ten",
        loop_ids: &[1, 1, 2, 2, 2, 2, 2, 2, 2, 2],
        projector_ids: &[1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
    },
    ProjectedFixture {
        name: "three_loop_rank_six",
        loop_ids: &[1, 1, 2, 2, 3, 3],
        projector_ids: &[1, 2, 3, 1, 2, 3],
    },
    ProjectedFixture {
        name: "three_loop_rank_ten",
        loop_ids: &[1, 1, 2, 2, 2, 2, 3, 3, 3, 3],
        projector_ids: &[1, 1, 1, 1, 2, 2, 2, 2, 3, 3],
    },
    ProjectedFixture {
        name: "four_loop_rank_eight",
        loop_ids: &[1, 1, 2, 2, 3, 3, 4, 4],
        projector_ids: &[1, 2, 3, 4, 1, 2, 3, 4],
    },
    ProjectedFixture {
        name: "four_loop_rank_ten",
        loop_ids: &[1, 1, 2, 2, 3, 3, 4, 4, 4, 4],
        projector_ids: &[1, 1, 2, 2, 3, 3, 4, 4, 4, 4],
    },
];

const BENCHMARK_FIXTURES: &[ProjectedFixture] =
    &[PARITY_FIXTURES[9], PARITY_FIXTURES[5], PARITY_FIXTURES[2]];

fn vakint_with(method: TensorReductionMethod) -> TestVakint {
    get_vakint(VakintSettings {
        allow_unknown_integrals: true,
        clean_tmp_dir: true,
        tensor_reduction_method: method,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    })
}

fn feynkit_without_form() -> TestVakint {
    get_vakint(VakintSettings {
        allow_unknown_integrals: true,
        clean_tmp_dir: true,
        form_exe_path: "/definitely/not/a/form/executable".into(),
        tensor_reduction_method: TensorReductionMethod::FeynKit,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    })
}

fn feynkit_expanded_without_form() -> TestVakint {
    get_vakint(VakintSettings {
        allow_unknown_integrals: true,
        clean_tmp_dir: true,
        form_exe_path: "/definitely/not/a/form/executable".into(),
        tensor_reduction_method: TensorReductionMethod::FeynKit,
        use_dot_product_notation: false,
        ..VakintSettings::default()
    })
}

fn collect_dot_dummy_labels(expression: AtomView<'_>, labels: &mut Vec<usize>) {
    if let AtomView::Fun(function) = expression
        && function.get_symbol().get_stripped_name() == "dot_dummy_ind"
        && function.get_nargs() == 1
        && let Some(argument) = function.iter().next()
        && let Ok(label) = usize::try_from(argument)
    {
        labels.push(label);
    }
    for child in expression.children() {
        collect_dot_dummy_labels(child, labels);
    }
}

#[test]
fn tensor_reduction_defaults_to_formless_feynkit() {
    assert_eq!(
        TensorReductionMethod::default(),
        TensorReductionMethod::FeynKit
    );

    let settings = VakintSettings {
        allow_unknown_integrals: true,
        clean_tmp_dir: true,
        form_exe_path: "/definitely/not/a/form/executable".into(),
        use_dot_product_notation: true,
        ..VakintSettings::default()
    };
    assert_eq!(
        settings.tensor_reduction_method,
        TensorReductionMethod::FeynKit
    );

    let vakint = get_vakint(settings);
    let input =
        vakint_parse!("k(1,1)*k(1,2)*p(1,1)*p(1,2)*topo(default_feynkit_reduction)").unwrap();
    let expected =
        vakint_parse!("dot(k(1),k(1))*dot(p(1),p(1))/(4-2*ε)*topo(default_feynkit_reduction)")
            .unwrap();

    let output = vakint.tensor_reduce(input.as_view()).unwrap();
    let difference = (output - &expected).expand().together();
    assert!(
        difference.is_zero(),
        "default FeynKit reduction differs from its expected value: {difference}"
    );
}

fn projected_input(fixture: ProjectedFixture) -> Atom {
    assert_eq!(fixture.loop_ids.len(), fixture.projector_ids.len());
    let numerator = fixture
        .loop_ids
        .iter()
        .zip(fixture.projector_ids)
        .enumerate()
        .flat_map(|(position, (loop_id, projector_id))| {
            let index = position + 1;
            [
                format!("k({loop_id},{index})"),
                format!("p({projector_id},{index})"),
            ]
        })
        .collect::<Vec<_>>()
        .join("*");
    vakint_parse!(&format!("({numerator})*topo({})", fixture.name)).unwrap()
}

fn assert_equivalent_outputs(name: &str, alphaloop: &Atom, feynkit: &Atom) {
    let alphaloop = Vakint::convert_to_dot_notation(alphaloop.as_view());
    let feynkit = Vakint::convert_to_dot_notation(feynkit.as_view());
    let difference = (alphaloop.clone() - &feynkit).expand().together();
    assert!(
        difference.is_zero(),
        "tensor-reduction mismatch for {name}:\nAlphaLoop: {alphaloop}\nFeynKit: {feynkit}\nDifference: {difference}"
    );
}

fn assert_modes_match(alphaloop: &TestVakint, feynkit: &TestVakint, name: &str, input: &Atom) {
    let alphaloop_output = alphaloop.tensor_reduce(input.as_view()).unwrap();
    let feynkit_output = feynkit.tensor_reduce(input.as_view()).unwrap();
    assert_equivalent_outputs(name, &alphaloop_output, &feynkit_output);
}

#[test]
fn projected_high_rank_numerators_match_alphaloop_through_rank_ten() {
    let alphaloop = vakint_with(TensorReductionMethod::AlphaLoop);
    let feynkit = vakint_with(TensorReductionMethod::FeynKit);
    for fixture in PARITY_FIXTURES {
        assert_modes_match(
            &alphaloop,
            &feynkit,
            fixture.name,
            &projected_input(*fixture),
        );
    }
}

#[test]
fn focused_tensor_structures_match_alphaloop() {
    let alphaloop = vakint_with(TensorReductionMethod::AlphaLoop);
    let feynkit = vakint_with(TensorReductionMethod::FeynKit);
    let cases = [
        (
            "odd_rank_vanishes",
            "k(1,1)*k(1,2)*k(1,3)*p(1,1)*p(1,2)*p(1,3)*topo(odd_rank)",
        ),
        ("free_rank_two", "k(1,1)*k(1,2)*topo(free_rank_two)"),
        (
            "free_rank_four",
            "k(1,1)*k(1,2)*k(2,3)*k(2,4)*topo(free_rank_four)",
        ),
        (
            "existing_dot_and_scalar_prefactor",
            "(1+3*x)*dot(k(1),k(2))*k(1,1)*k(2,2)*p(1,1)*p(2,2)*topo(existing_dot)",
        ),
        (
            "sum_of_tensor_monomials",
            "(2*a*k(1,1)*k(1,2)*p(1,1)*p(2,2)-3*b*k(2,3)*k(2,4)*p(1,3)*p(1,4))*topo(tensor_sum)",
        ),
        (
            "compact_mixed_dot_odd_rank",
            "dot(p(1),k(1))*topo(compact_mixed_dot_odd_rank)",
        ),
        (
            "compact_mixed_dots_rank_two",
            "dot(p(1),k(1))*dot(p(2),k(1))*topo(compact_mixed_dots_rank_two)",
        ),
    ];

    for (name, input) in cases {
        let input = vakint_parse!(input).unwrap();
        assert_modes_match(&alphaloop, &feynkit, name, &input);
    }
}

#[test]
fn feynkit_mode_contracts_explicit_metrics_without_form() {
    // AlphaLoop treats metrics as external and its FORM source explicitly
    // marks this contraction as unsafe, so it is not an oracle for this case.
    let feynkit = feynkit_without_form();
    let input = vakint_parse!("g(1,2)*k(1,1)*k(2,2)*topo(feynkit_explicit_metric)").unwrap();
    let expected = vakint_parse!("dot(k(1),k(2))*topo(feynkit_explicit_metric)").unwrap();
    let output = feynkit.tensor_reduce(input.as_view()).unwrap();
    let difference = (expected - &output).expand().together();
    assert!(
        difference.is_zero(),
        "explicit metric contraction differs from its expected value: {difference}"
    );
}

#[test]
fn feynkit_mode_retargets_full_spenso_slots_without_nesting() {
    let feynkit = feynkit_without_form();
    let input = vakint_parse!(
        "k(1,spenso::mink(dim,mu))*k(1,spenso::mink(4,nu))
         *p(7,spenso::mink(dim,mu))*p(8,spenso::mink(4,nu))
         *topo(mixed_full_spenso_slots)"
    )
    .unwrap();
    let expected =
        vakint_parse!("dot(k(1),k(1))*dot(p(7),p(8))/(4-2*ε)*topo(mixed_full_spenso_slots)")
            .unwrap();

    let output = feynkit.tensor_reduce(input.as_view()).unwrap();
    let difference = (output.clone() - expected).expand().together();
    assert!(
        difference.is_zero(),
        "mixed full Spenso slots differ by {difference}: {output}"
    );
    let canonical = output.to_canonical_string();
    assert!(
        !canonical.contains("mink(4+-2*vakint::ε,spenso::")
            && !canonical.contains("mink(vakint::dim,spenso::"),
        "a full Spenso slot was nested inside another slot: {canonical}"
    );
}

#[test]
fn feynkit_mode_retargets_full_metric_slots_without_nesting() {
    let feynkit = feynkit_without_form();
    let input = vakint_parse!(
        "g(spenso::mink(dim,mu),spenso::mink(4,nu))
         *k(1,spenso::mink(dim,mu))*k(2,spenso::mink(4,nu))
         *topo(full_spenso_metric)"
    )
    .unwrap();
    let expected = vakint_parse!("dot(k(1),k(2))*topo(full_spenso_metric)").unwrap();

    let output = feynkit.tensor_reduce(input.as_view()).unwrap();
    let difference = (output.clone() - expected).expand().together();
    assert!(
        difference.is_zero(),
        "full Spenso metric slots differ by {difference}: {output}"
    );
}

#[test]
fn tagged_spatial_components_are_not_minkowski_slots() {
    let _q = vector_symbol!("vakint_tensor_bridge_test::Q_spatial");
    let feynkit = feynkit_without_form();
    let input = vakint_parse!(
        "vakint_tensor_bridge_test::Q_spatial(1,spenso::cind(1))
         *topo(tagged_spatial_component)"
    )
    .unwrap();

    let output = feynkit.tensor_reduce(input.as_view()).unwrap();
    let difference = (output.clone() - &input).expand().together();
    assert!(
        difference.is_zero(),
        "a tagged spatial component was changed by {difference}: {output}"
    );
    let canonical = output.to_canonical_string();
    assert!(
        !canonical.contains("mink("),
        "a Minkowski representation was injected around cind: {canonical}"
    );
}

#[test]
fn compact_tagged_vectors_project_back_to_a_compact_vakint_dot() {
    let _q = vector_symbol!("vakint_tensor_bridge_test::Q_compact");
    let _p = vector_symbol!("vakint_tensor_bridge_test::P_compact");
    let feynkit = feynkit_without_form();
    let input = vakint_parse!(
        "dot(k(1,spenso::mink(4)),vakint_tensor_bridge_test::Q_compact(7,spenso::mink(4)))
         *dot(k(1,spenso::mink(4)),vakint_tensor_bridge_test::P_compact(3,spenso::mink(4)))
         *topo(compact_tagged_vectors)"
    )
    .unwrap();
    let expected = vakint_parse!(
        "dot(k(1),k(1))
         *dot(vakint_tensor_bridge_test::Q_compact(7,spenso::mink(4)),
              vakint_tensor_bridge_test::P_compact(3,spenso::mink(4)))
         /(4-2*ε)*topo(compact_tagged_vectors)"
    )
    .unwrap();

    let output = feynkit.tensor_reduce(input.as_view()).unwrap();
    let difference = (output.clone() - expected).expand().together();
    assert!(
        difference.is_zero(),
        "compact tagged-vector projection differs by {difference}: {output}"
    );
    let canonical = output.to_canonical_string();
    assert!(
        !canonical.contains("feynkit_external_vector")
            && !canonical.contains("feynkit_opaque_index")
            && !canonical.contains("dot_dummy_ind"),
        "temporary FeynKit bridge syntax leaked into the result: {canonical}"
    );
}

#[test]
fn native_spenso_dots_with_loop_momenta_are_tensor_reduced() {
    let feynkit = feynkit_without_form();
    let input = vakint_parse!(
        "spenso::dot(k(1,spenso::mink(4-2*ε)),p(7,spenso::mink(4-2*ε)))
         *spenso::dot(k(1,spenso::mink(4-2*ε)),p(8,spenso::mink(4-2*ε)))
         *topo(native_spenso_dots)"
    )
    .unwrap();
    let expected =
        vakint_parse!("dot(k(1),k(1))*dot(p(7),p(8))/(4-2*ε)*topo(native_spenso_dots)").unwrap();

    let output = feynkit.tensor_reduce(input.as_view()).unwrap();
    let difference = (output.clone() - expected).expand().together();
    assert!(
        difference.is_zero(),
        "native Spenso dots bypassed tensor reduction by {difference}: {output}"
    );
}

#[test]
fn non_minkowski_spenso_dots_are_preserved() {
    let feynkit = feynkit_without_form();
    let input = vakint_parse!(
        "spenso::dot(bridge_test::C(spenso::coad(8)),bridge_test::D(spenso::coad(8)))
         *topo(non_minkowski_spenso_dot)"
    )
    .unwrap();

    let output = feynkit.tensor_reduce(input.as_view()).unwrap();
    let difference = (output.clone() - &input).expand().together();
    assert!(
        difference.is_zero(),
        "a non-Minkowski Spenso dot was changed by {difference}: {output}"
    );
    let canonical = output.to_canonical_string();
    assert!(
        !canonical.contains("dot_dummy_ind"),
        "a non-Minkowski Spenso dot was rewritten as a Lorentz contraction: {canonical}"
    );
}

#[test]
fn malformed_spenso_loop_dots_are_rejected() {
    let feynkit = feynkit_without_form();
    let input = vakint_parse!(
        "spenso::dot(k(1,spenso::mink(4)),bridge_test::C(spenso::coad(8)))
         *topo(malformed_spenso_loop_dot)"
    )
    .unwrap();

    let error = feynkit.tensor_reduce(input.as_view()).unwrap_err();
    assert!(
        matches!(error, VakintError::FeynKitLoopDot(_)),
        "unexpected malformed-loop-dot error: {error}"
    );
}

#[test]
fn gamma_like_opaque_projectors_use_a_paired_boundary_dummy() {
    let feynkit = feynkit_without_form();
    let input = vakint_parse!(
        "k(1,spenso::mink(4,mu))*k(1,spenso::mink(4,nu))
         *spenso::gamma(spenso::bis(4,a),spenso::bis(4,b),spenso::mink(4,mu))
         *p(7,spenso::mink(4,nu))*topo(gamma_like_opaque_projector)"
    )
    .unwrap();
    let expected = vakint_parse!(
        "dot(k(1),k(1))/(4-2*ε)
         *spenso::gamma(spenso::bis(4,a),spenso::bis(4,b),spenso::mink(4,mu))
         *g(spenso::mink(4,mu),dot_dummy_ind(1))*p(7,dot_dummy_ind(1))
         *topo(gamma_like_opaque_projector)"
    )
    .unwrap();

    let feynkit_output = feynkit.tensor_reduce(input.as_view()).unwrap();
    let difference = (feynkit_output.clone() - expected).expand().together();
    assert!(
        difference.is_zero(),
        "gamma-like projector did not use a paired boundary dummy: {difference}; FeynKit: {feynkit_output}"
    );
    let canonical = feynkit_output.to_canonical_string();
    assert!(
        !canonical.contains("feynkit_external_vector")
            && !canonical.contains("feynkit_opaque_index"),
        "temporary FeynKit bridge syntax leaked into gamma-like output: {canonical}"
    );
}

#[test]
fn gamma_like_projector_with_a_tagged_vector_stays_native_spenso() {
    let _q = vector_symbol!("vakint_tensor_bridge_test::Q_gamma_projector");
    let feynkit = feynkit_without_form();
    let input = vakint_parse!(
        "k(1,spenso::mink(4,mu))*k(1,spenso::mink(4,nu))
         *spenso::gamma(spenso::bis(4,a),spenso::bis(4,b),spenso::mink(4,mu))
         *vakint_tensor_bridge_test::Q_gamma_projector(7,spenso::mink(4,nu))
         *topo(gamma_tagged_vector_projector)"
    )
    .unwrap();
    let expected = vakint_parse!(
        "dot(k(1),k(1))/(4-2*ε)
         *spenso::gamma(spenso::bis(4,a),spenso::bis(4,b),spenso::mink(4,mu))
         *spenso::g(spenso::mink(4,mu),spenso::mink(4,dot_dummy_ind(1)))
         *vakint_tensor_bridge_test::Q_gamma_projector(
             7,spenso::mink(4,dot_dummy_ind(1)))
         *topo(gamma_tagged_vector_projector)"
    )
    .unwrap();

    let output = feynkit.tensor_reduce(input.as_view()).unwrap();
    let difference = (output.clone() - expected).expand().together();
    assert!(
        difference.is_zero(),
        "gamma/tagged-vector projector left native Spenso form by {difference}: {output}"
    );
    let canonical = output.to_canonical_string();
    assert!(
        !canonical.contains("feynkit_external_vector")
            && !canonical.contains("feynkit_opaque_index"),
        "temporary FeynKit bridge syntax leaked into tagged-vector output: {canonical}"
    );
}

#[test]
fn opaque_and_free_projector_legs_keep_full_spenso_slots() {
    let feynkit = feynkit_without_form();
    let input = vakint_parse!(
        "k(1,spenso::mink(4,mu))*k(1,nu)
         *spenso::gamma(spenso::bis(4,a),spenso::bis(4,b),spenso::mink(4,mu))
         *topo(opaque_and_free_projector)"
    )
    .unwrap();
    let expected = vakint_parse!(
        "dot(k(1),k(1))/(4-2*ε)
         *spenso::gamma(spenso::bis(4,a),spenso::bis(4,b),spenso::mink(4,mu))
         *g(spenso::mink(4,mu),spenso::mink(4-2*ε,nu))
         *topo(opaque_and_free_projector)"
    )
    .unwrap();

    let output = feynkit.tensor_reduce(input.as_view()).unwrap();
    let difference = (output.clone() - expected).expand().together();
    assert!(
        difference.is_zero(),
        "opaque/free projector lost a full Spenso slot by {difference}: {output}"
    );
}

#[test]
fn opaque_projector_dummies_are_unique_and_noncolliding() {
    let feynkit = feynkit_without_form();
    let input = vakint_parse!(
        "bridge_test::marker(dot_dummy_ind(7))*(
           k(1,spenso::mink(4,mu))*k(1,spenso::mink(4,nu))
            *spenso::gamma(spenso::bis(4,a),spenso::bis(4,b),spenso::mink(4,mu))
            *p(7,spenso::mink(4,nu))
          +k(1,spenso::mink(4,rho))*k(1,spenso::mink(4,sigma))
            *spenso::gamma(spenso::bis(4,c),spenso::bis(4,d),spenso::mink(4,rho))
            *p(8,spenso::mink(4,sigma)))
         *topo(unique_opaque_projector_dummies)"
    )
    .unwrap();

    let output = feynkit.tensor_reduce(input.as_view()).unwrap();
    let mut labels = Vec::new();
    collect_dot_dummy_labels(output.as_view(), &mut labels);
    labels.sort_unstable();
    assert_eq!(labels, vec![7, 7, 8, 8, 9, 9]);
    let canonical = output.to_canonical_string();
    assert!(
        !canonical.contains("feynkit_external_vector")
            && !canonical.contains("feynkit_opaque_index"),
        "temporary FeynKit bridge syntax leaked across additive terms: {canonical}"
    );
}

#[test]
fn expanded_mode_indexes_the_existing_tagged_vector_representation() {
    let _q = vector_symbol!("vakint_tensor_bridge_test::Q_expanded");
    let _p = vector_symbol!("vakint_tensor_bridge_test::P_expanded");
    let feynkit = feynkit_expanded_without_form();
    let input = vakint_parse!(
        "dot(k(1),vakint_tensor_bridge_test::Q_expanded(7,spenso::mink(4)))
         *dot(k(1),vakint_tensor_bridge_test::P_expanded(3,spenso::mink(4)))
         *topo(expanded_tagged_vectors)"
    )
    .unwrap();

    let output = feynkit.tensor_reduce(input.as_view()).unwrap();
    let canonical = output.to_canonical_string();
    assert!(
        canonical.contains("mink(4,vakint::{}::dot_dummy_ind"),
        "expanded tagged vectors did not receive an index in their existing representation: {canonical}"
    );
    assert!(
        !canonical.contains("mink(4),vakint::{}::dot_dummy_ind")
            && !canonical.contains("feynkit_external_vector")
            && !canonical.contains("feynkit_opaque_index"),
        "expanded tagged vectors contain malformed or temporary bridge syntax: {canonical}"
    );
}

#[test]
fn opaque_slot_provenance_is_local_to_each_additive_term() {
    let feynkit = feynkit_without_form();
    let input = vakint_parse!(
        "(k(1,spenso::mink(4,mu))*k(1,spenso::mink(4,nu))
            *bridge_test::A(spenso::mink(4,mu))*bridge_test::B(spenso::mink(4,nu))
          +k(1,spenso::mink(dim,mu))*k(1,spenso::mink(dim,nu))
            *bridge_test::C(spenso::mink(dim,mu))*bridge_test::E(spenso::mink(dim,nu)))
         *topo(opaque_slot_provenance)"
    )
    .unwrap();
    let expected = vakint_parse!(
        "dot(k(1),k(1))/(4-2*ε)
         *(dot(bridge_test::A(spenso::mink(4)),bridge_test::B(spenso::mink(4)))
          +dot(bridge_test::C(spenso::mink(dim)),bridge_test::E(spenso::mink(dim))))
         *topo(opaque_slot_provenance)"
    )
    .unwrap();

    let output = feynkit.tensor_reduce(input.as_view()).unwrap();
    let difference = (output.clone() - expected).expand().together();
    assert!(
        difference.is_zero(),
        "opaque-slot provenance crossed additive terms by {difference}: {output}"
    );
    let canonical = output.to_canonical_string();
    assert!(
        !canonical.contains("g(spenso::mink")
            && !canonical.contains("feynkit_opaque_index")
            && !canonical.contains("dot_dummy_ind"),
        "invalid bridge metric syntax leaked into the result: {canonical}"
    );
}

#[test]
fn untagged_vector_provenance_distinguishes_dimensions_within_one_term() {
    let feynkit = feynkit_without_form();
    let input = vakint_parse!(
        "k(1,spenso::mink(4,mu))*k(1,spenso::mink(dim,nu))
         *bridge_test::V(spenso::mink(4,mu))*bridge_test::V(spenso::mink(dim,nu))
         *topo(untagged_vector_dimensions)"
    )
    .unwrap();
    let expected = vakint_parse!(
        "dot(k(1),k(1))/(4-2*ε)
         *dot(bridge_test::V(spenso::mink(4)),bridge_test::V(spenso::mink(dim)))
         *topo(untagged_vector_dimensions)"
    )
    .unwrap();

    let output = feynkit.tensor_reduce(input.as_view()).unwrap();
    let difference = (output.clone() - expected).expand().together();
    assert!(
        difference.is_zero(),
        "untagged vector dimensions collided by {difference}: {output}"
    );
}

#[test]
fn feynkit_expanded_mode_restores_vakint_vectors_and_metrics() {
    let feynkit = feynkit_expanded_without_form();
    let input = vakint_parse!("(k(1,1)*k(1,2)*p(7,1)*p(8,2)+k(2,3)*k(2,4))*topo(feynkit_expanded)")
        .unwrap();
    let expected = vakint_parse!(
        "(dot(k(1),k(1))*dot(p(7),p(8))+dot(k(2),k(2))*g(3,4))/(4-2*ε)*topo(feynkit_expanded)"
    )
    .unwrap();

    let output = feynkit.tensor_reduce(input.as_view()).unwrap();
    let canonical = output.to_canonical_string();
    assert!(
        !canonical.contains("spenso::"),
        "Spenso bridge syntax leaked into Vakint output: {canonical}"
    );
    assert!(
        !canonical.contains("::dot("),
        "dot products were not expanded as requested: {canonical}"
    );
    assert!(
        canonical.contains("::g("),
        "a residual projector metric was not restored to Vakint syntax: {canonical}"
    );
    assert!(
        canonical.contains("::k(") && canonical.contains("::p("),
        "compact loop or external vectors were not restored to indexed Vakint syntax: {canonical}"
    );

    let redotted = Vakint::convert_to_dot_notation(output.as_view());
    let difference = (expected - &redotted).expand().together();
    assert!(
        difference.is_zero(),
        "expanded FeynKit output differs after restoring dot notation: {difference}"
    );
}

#[test]
fn feynkit_mode_is_formless_and_reaches_ranks_twelve_and_twenty() {
    let feynkit = feynkit_without_form();
    let fixtures = [
        ProjectedFixture {
            name: "feynkit_only_rank_twelve",
            loop_ids: &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            projector_ids: &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        },
        ProjectedFixture {
            name: "feynkit_only_rank_twenty",
            loop_ids: &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            projector_ids: &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        },
    ];

    for fixture in fixtures {
        let output = feynkit
            .tensor_reduce(projected_input(fixture).as_view())
            .unwrap();
        let pair_count = fixture.loop_ids.len() / 2;
        let dimension = vakint_parse!("4-2*ε").unwrap();
        let denominator = (0..pair_count).fold(Atom::one(), |product, offset| {
            product * (dimension.clone() + Atom::num(2 * offset as i64))
        });
        let pairing_multiplicity = (1..=pair_count).map(|pair| 2 * pair - 1).product::<usize>();
        let k_squared = vakint_parse!("dot(k(1),k(1))").unwrap();
        let p_squared = vakint_parse!("dot(p(1),p(1))").unwrap();
        let topology = vakint_parse!(&format!("topo({})", fixture.name)).unwrap();
        let expected = Atom::num(pairing_multiplicity as i64)
            * k_squared.pow(Atom::num(pair_count as i64))
            * p_squared.pow(Atom::num(pair_count as i64))
            / denominator
            * topology;
        let output = Vakint::convert_to_dot_notation(output.as_view());
        let difference = (output - &expected).together();
        assert!(
            difference.is_zero(),
            "{} differs from the closed isotropic rank-{} result: {difference}",
            fixture.name,
            fixture.loop_ids.len()
        );
    }
}

fn timed_reduction(vakint: &TestVakint, input: &Atom) -> (Duration, Atom) {
    let started = Instant::now();
    let output = vakint.tensor_reduce(input.as_view()).unwrap();
    let elapsed = started.elapsed();
    (elapsed, black_box(output))
}

fn median(mut samples: Vec<Duration>) -> Duration {
    samples.sort_unstable();
    samples[samples.len() / 2]
}

fn milliseconds(duration: Duration) -> f64 {
    duration.as_secs_f64() * 1_000.0
}

#[test]
#[ignore = "manual end-to-end runtime report; invokes FORM repeatedly"]
fn tensor_reduction_runtime_report() {
    const REPETITIONS: usize = 3;

    let alphaloop = vakint_with(TensorReductionMethod::AlphaLoop);
    let feynkit = vakint_with(TensorReductionMethod::FeynKit);
    eprintln!(
        "case                         alpha first  feynkit first  alpha median  feynkit median  speedup"
    );

    for fixture in BENCHMARK_FIXTURES {
        let input = projected_input(*fixture);
        let (alphaloop_first, alphaloop_output) = timed_reduction(&alphaloop, &input);
        let (feynkit_first, feynkit_output) = timed_reduction(&feynkit, &input);
        assert_equivalent_outputs(fixture.name, &alphaloop_output, &feynkit_output);

        let mut alphaloop_samples = Vec::with_capacity(REPETITIONS);
        let mut feynkit_samples = Vec::with_capacity(REPETITIONS);
        for repetition in 0..REPETITIONS {
            if repetition.is_multiple_of(2) {
                alphaloop_samples.push(timed_reduction(&alphaloop, &input).0);
                feynkit_samples.push(timed_reduction(&feynkit, &input).0);
            } else {
                feynkit_samples.push(timed_reduction(&feynkit, &input).0);
                alphaloop_samples.push(timed_reduction(&alphaloop, &input).0);
            }
        }
        let alphaloop_median = median(alphaloop_samples);
        let feynkit_median = median(feynkit_samples);
        let speedup = alphaloop_median.as_secs_f64() / feynkit_median.as_secs_f64();
        eprintln!(
            "{:<28} {:>10.3} ms {:>11.3} ms {:>11.3} ms {:>13.3} ms {:>8.2}x",
            fixture.name,
            milliseconds(alphaloop_first),
            milliseconds(feynkit_first),
            milliseconds(alphaloop_median),
            milliseconds(feynkit_median),
            speedup,
        );
    }
}
