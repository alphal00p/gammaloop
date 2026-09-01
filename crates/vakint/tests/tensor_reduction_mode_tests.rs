mod test_utils;

use std::{
    hint::black_box,
    time::{Duration, Instant},
};

use symbolica::atom::{Atom, AtomCore};
use test_utils::{TestVakint, get_vakint};
use vakint::{TensorReductionMethod, Vakint, VakintSettings, vakint_parse};

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
