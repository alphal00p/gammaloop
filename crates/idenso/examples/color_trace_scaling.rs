//! Colour trace reduction timings, with admission and destruction excluded.
//!
//! cargo run -p idenso --release --example color_trace_scaling -- 8 3
//! Add `--features reference-cases` for a separate diagnostic phase capture.
//! Inputs keep all group invariants symbolic. No tensor sums are expanded.
//! Correctness is checked independently in spenso-hep-lib/color_trace_validation.

use std::{hint::black_box, time::Instant};

use idenso::{
    color::ColorSimplifySettings,
    representations::{ColorAdjoint, ColorFundamental, initialize},
    tensor::{AlgebraContraction, AlgebraSettings, ReductionStatus, SymbolicTensor},
};
use spenso::{network::tags::SPENSO_TAG, structure::representation::RepName};
use symbolica::{atom::Atom, symbol};

fn main() {
    let arguments = std::env::args().collect::<Vec<_>>();
    let maximum: usize = arguments.get(1).map_or(8, |value| value.parse().unwrap());
    let samples: usize = arguments.get(2).map_or(3, |value| value.parse().unwrap());
    initialize();
    let adjoint = ColorAdjoint {}.new_rep(symbol!("color_scaling::Na"));
    let fundamental = ColorFundamental {}.new_rep(symbol!("color_scaling::Nc"));
    let settings = AlgebraSettings {
        gamma: None,
        color: Some(ColorSimplifySettings::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    println!("representation,length,sample,reduction_seconds,rerun_seconds,bytes");
    for is_adjoint in [false, true] {
        let representation = if is_adjoint {
            adjoint.to_symbolic([])
        } else {
            fundamental.to_symbolic([])
        };
        let name = if is_adjoint { "adjoint" } else { "fundamental" };
        for length in 5..=maximum {
            let expression = spenso::trace!(&representation; (0..length).map(|position| {
                let slot = adjoint.to_symbolic([Atom::var(symbol!(&format!("color_scaling::a{position}")))]);
                if is_adjoint {
                    idenso::color_f!(Atom::var(SPENSO_TAG.chain_in), Atom::var(SPENSO_TAG.chain_out), slot)
                } else {
                    idenso::color_t!(slot)
                }
            }));
            for sample in 0..samples {
                let source = SymbolicTensor::infer(expression.clone()).unwrap();
                let start = Instant::now();
                let result = black_box(source.simplify_algebra(&settings).unwrap());
                let elapsed = start.elapsed().as_secs_f64();
                assert_eq!(result.reduction_status(), ReductionStatus::Complete);
                assert_eq!(source.structure(), result.structure());
                let start = Instant::now();
                let repeated = black_box(result.simplify_algebra(&settings).unwrap());
                let rerun = start.elapsed().as_secs_f64();
                assert_eq!(repeated, result);
                assert_eq!(repeated.reduction_status(), ReductionStatus::Complete);
                println!(
                    "{name},{length},{sample},{elapsed:.9},{rerun:.9},{}",
                    result.expression().as_view().get_byte_size()
                );
            }
            #[cfg(feature = "reference-cases")]
            {
                let source = SymbolicTensor::infer(expression.clone()).unwrap();
                let (_, phases) = idenso::reference_cases::timing::measure("colour trace", || {
                    source.simplify_algebra(&settings).unwrap()
                });
                eprintln!("{name},{length},{phases}");
            }
        }
    }
}
