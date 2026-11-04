// Temporary comparison harness for PR174, retained as review evidence.
// Attach as a test-only child module of utils; remove the tracked copy afterward.
// Baseline bodies are copied from 60eb9f30 and 033910560 respectively.
// Timings cover complete boolean comparisons, including inherited norm checks.
use super::*;
use std::hint;
use std::time::Instant;

#[inline(never)]
fn base_compare<T: FloatLike>(lhs: &Complex<F<T>>, rhs: &Complex<F<T>>, tolerance: &F<T>) -> bool {
    if !lhs.norm().re.approx_eq(&rhs.norm().re, tolerance) {
        return false;
    } else if lhs.norm().is_zero() || rhs.norm().is_zero() {
        return true;
    }
    let two_pi = lhs.re.PI() + lhs.re.PI();
    let arg_lhs = lhs.arg().rem_euclid(&two_pi);
    let arg_rhs = rhs.arg().rem_euclid(&two_pi);
    if !arg_lhs.approx_eq(&arg_rhs, tolerance) {
        return false;
    }
    true
}

fn published_phase_differences<T: FloatLike>(
    lhs: &Complex<F<T>>,
    rhs: &Complex<F<T>>,
) -> (F<T>, F<T>) {
    // Binary scaling retains the represented directions. Dividing each
    // vector by an arbitrary maximum can round distinct ratios to equality.
    // A coarse outward enclosure suffices to choose a bounding binade;
    // no source component is converted through that enclosure.
    let scale = |value: &Complex<F<T>>| {
        let maximum = value.re.abs().max(value.im.abs());
        let (lower, _) = maximum.0.mpfr_enclosure(2);
        let exponent = lower
            .get_exp()
            .expect("norm checks admit finite nonzero directions");
        // Keep roughly half a significand of exponent headroom in each
        // direction: even a subnormal angular gap then has a normal
        // determinant, rather than disappearing in products of fractions.
        let headroom = i32::try_from(maximum.0.get_precision().div_ceil(2))
            .expect("working precision fits a binary scale exponent");
        maximum.from_usize(2).powi(headroom + 1 - exponent)
    };
    let lhs_scale = scale(lhs);
    let rhs_scale = scale(rhs);
    let lhs_re = &lhs.re * &lhs_scale;
    let lhs_im = &lhs.im * &lhs_scale;
    let rhs_re = &rhs.re * &rhs_scale;
    let rhs_im = &rhs.im * &rhs_scale;
    // Compensate the rounded product in the determinant when nearby
    // directions cancel; binary scaling keeps both products bounded.
    let product = &lhs_re * &rhs_im;
    let product_error = lhs_re.mul_add(&rhs_im, &(-&product));
    let cross = lhs_im.mul_add(&rhs_re, &(-product)) - product_error;
    let dot = lhs_re.mul_add(&rhs_re, &(&lhs_im * &rhs_im));
    // Keep the internal atan2 ratio bounded too: the double-double
    // backend can overflow when a near-perpendicular dot is subnormal.
    let difference = if cross.abs() > dot.abs() {
        let angle = cross.PIHALF() - dot.atan2(&cross.abs());
        if cross < cross.zero() { -angle } else { angle }
    } else {
        cross.atan2(&dot)
    };
    let rhs_phase = rhs.arg().rem_euclid(&rhs.re.TAU());
    let relative = if rhs_phase.is_zero() {
        difference.clone()
    } else if dot > dot.zero()
        && !cross.is_zero()
        && difference.0.is_finite()
        && (difference.is_zero() || difference.0.needs_rescaling())
    {
        // In this underflow range atan(t)/t rounds to one. Combine the
        // RHS normalization before rounding the angle: a subnormal (or
        // unrepresentable) difference can still give a finite ratio.
        &cross / (&dot * rhs_phase)
    } else {
        &difference / rhs_phase
    };
    (difference, relative)
}

#[inline(never)]
fn published_compare<T: FloatLike>(
    lhs: &Complex<F<T>>,
    rhs: &Complex<F<T>>,
    tolerance: &F<T>,
) -> bool {
    if !lhs.norm().re.approx_eq(&rhs.norm().re, tolerance) {
        return false;
    } else if lhs.norm().is_zero() || rhs.norm().is_zero() {
        return true;
    }
    let (_, relative_phase) = published_phase_differences(lhs, rhs);
    relative_phase.approx_eq(&relative_phase.zero(), tolerance)
}

#[inline(never)]
fn revised_compare<T: FloatLike>(
    lhs: &Complex<F<T>>,
    rhs: &Complex<F<T>>,
    tolerance: &F<T>,
) -> bool {
    lhs.approx_eq(rhs, tolerance)
}

type Compare<T> = fn(&Complex<F<T>>, &Complex<F<T>>, &F<T>) -> bool;

struct Input<T: FloatLike> {
    lhs: Complex<F<T>>,
    rhs: Complex<F<T>>,
    tolerance: F<T>,
}

fn measure<T: FloatLike>(
    compare: Compare<T>,
    inputs: &[Input<T>],
    iterations: usize,
) -> (f64, usize) {
    let compare = hint::black_box(compare);
    let inputs = hint::black_box(inputs);
    let start = Instant::now();
    let mut checksum = 0;
    for index in 0..iterations {
        let input = hint::black_box(&inputs[index % inputs.len()]);
        checksum += usize::from(hint::black_box(compare(
            hint::black_box(&input.lhs),
            hint::black_box(&input.rhs),
            hint::black_box(&input.tolerance),
        )));
    }
    let elapsed = start.elapsed().as_secs_f64();
    (elapsed * 1e9 / iterations as f64, hint::black_box(checksum))
}

fn benchmark_backend<T: FloatLike>(backend: &str, prototype: F<T>, iterations: usize) {
    let one = prototype.one();
    let two = prototype.from_usize(2);
    let delta = two.powi(-40);
    let tolerance = two.powi(-30);
    let mut ordinary = Vec::new();
    let mut wrap = Vec::new();
    let mut norm_reject = Vec::new();
    // Sixteen inputs prevent a single constant comparison from dominating.
    // Powers of two vary magnitude while remaining valid for every backend.
    for index in 0..16 {
        let scale = two.powi([-200, -50, 0, 50, 200][index % 5]);
        let gap = &delta * prototype.from_usize(index + 1);
        let shifted = &scale * (&one + &gap);
        ordinary.push(Input {
            lhs: Complex::new(shifted.clone(), scale.clone()),
            rhs: Complex::new(scale.clone(), shifted),
            tolerance: tolerance.clone(),
        });
        wrap.push(Input {
            lhs: Complex::new(scale.clone(), -(&scale * &gap)),
            rhs: Complex::new(scale.clone(), prototype.zero()),
            tolerance: tolerance.clone(),
        });
        norm_reject.push(Input {
            lhs: Complex::new(&scale * &two, prototype.zero()),
            rhs: Complex::new(scale, prototype.zero()),
            tolerance: tolerance.clone(),
        });
    }
    let variants: [(&str, Compare<T>); 3] = [
        ("base_60eb9f30", base_compare::<T>),
        ("published_03391056", published_compare::<T>),
        ("revised", revised_compare::<T>),
    ];
    for (case, inputs) in [
        ("ordinary_nearby", ordinary),
        ("wrap_crossing", wrap),
        ("norm_reject", norm_reject),
    ] {
        for (variant, compare) in variants {
            // The base's false rejection on wrap crossing is the intended
            // semantic change. All other workload decisions must agree.
            let expected =
                case != "norm_reject" && !(case == "wrap_crossing" && variant == "base_60eb9f30");
            for input in &inputs {
                assert_eq!(
                    compare(&input.lhs, &input.rhs, &input.tolerance),
                    expected,
                    "unexpected decision: backend={backend}, case={case}, variant={variant}"
                );
            }
            hint::black_box(measure(compare, &inputs, iterations / 20));
        }
        for round in 0..9 {
            // Rotate ordering so one implementation is not always measured
            // first or last under thermal/load drift.
            for offset in 0..variants.len() {
                let (variant, compare) = variants[(round + offset) % variants.len()];
                let (ns_per_call, checksum) = measure(compare, &inputs, iterations);
                let expected_checksum = if case != "norm_reject"
                    && !(case == "wrap_crossing" && variant == "base_60eb9f30")
                {
                    iterations
                } else {
                    0
                };
                assert_eq!(checksum, expected_checksum);
                println!(
                    "phase_bench,{backend},{case},{variant},{round},{iterations},{ns_per_call:.3},{checksum}"
                );
            }
        }
    }
}

#[test]
fn phase_comparison_baseline_controls() {
    for (real, imaginary) in [
        (2.0_f64.powi(500), 2.0_f64.powi(-574)),
        (
            2.0_f64.powi(500) * (1.0 - 2.0_f64.powi(-40)),
            2.0_f64.powi(-575).next_down(),
        ),
    ] {
        let lhs = Complex::new(F(real), F(imaginary.next_up()));
        let rhs = Complex::new(F(real), F(imaginary));
        let tolerance = F(2.0_f64.powi(-56));
        assert!(rhs.norm().re.0.is_normal());
        assert_eq!(lhs.norm().re, rhs.norm().re);
        assert!(base_compare(&lhs, &rhs, &tolerance));
        assert!(published_compare(&lhs, &rhs, &tolerance));
        assert!(!revised_compare(&lhs, &rhs, &tolerance));
    }
    let tiny = Complex::new(F(2.0_f64.powi(-1000)), F(0.0));
    assert_eq!(tiny.norm().re.0, 0.0);
    for compare in [base_compare, published_compare, revised_compare] {
        assert!(compare(&tiny, &tiny, &F(1e-10)));
    }
}

#[test]
#[ignore = "manual PR174 microbenchmark; run alone with --nocapture"]
fn phase_comparison_microbenchmark() {
    println!("record,backend,case,variant,round,iterations,ns_per_call,checksum");
    benchmark_backend("f64", F(0.0_f64), 200_000);
    benchmark_backend::<QuadFloat>("double_double", Default::default(), 5_000);
    benchmark_backend::<SamplingFloat>("mpfr256", Default::default(), 2_000);
}
