//! Algebra-only ports of the boosted OS regressions. Feynkit has no OS component
//! operation: these certify its tensor algebra, not forest or momentum ownership.

use std::sync::Once;

use idenso::{
    dirac::{AGS, GammaSimplifySettings},
    representations::{Bispinor, initialize},
    tensor::{AlgebraSettings, SymbolicTensor},
};
use spenso::{
    network::{
        ExecutionResult, Sequential, SmallestDegree,
        library::symbolic::{ETS, ExplicitKey, TensorLibrary},
        parsing::{ParseSettings, StrictTensorFilter},
        tags::SPENSO_TAG,
    },
    structure::{
        abstract_index::{AIND_SYMBOLS, AbstractIndex},
        representation::{Minkowski, RepName},
        slot::IsAbstractSlot,
    },
    symbolic_parallelism::{SymbolicParallelism, set_symbolica_rayon_enabled},
    tensors::{
        data::DenseTensor,
        parametric::{MixedTensor, ParamTensor},
    },
};
use spenso_hep_lib::{FUN_LIB, HepNet, hep_lib_atom};
use symbolica::{
    atom::{Atom, AtomCore},
    function, symbol,
};

type Library = TensorLibrary<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>;
type Matrix = [[Atom; 4]; 4];

struct ProjectorEvaluation {
    library: Library,
    spin: [Atom; 4],
    lorentz: [Atom; 2],
    vectors: [Atom; 2],
}

impl ProjectorEvaluation {
    fn new() -> Self {
        static INITIALIZE: Once = Once::new();
        INITIALIZE.call_once(|| {
            initialize();
            set_symbolica_rayon_enabled(SymbolicParallelism::Serial);
        });
        let mut library = hep_lib_atom();
        let mink = Minkowski {}.new_rep(4);
        let bis = Bispinor {}.new_rep(4);
        // Keep all finite contractions exact, including explicit metrics that
        // would otherwise use the mixed library's generic floating-point data.
        for (representation, lorentzian) in [(mink.to_lib(), true), (bis.to_lib(), false)] {
            let key = ExplicitKey::from_iter([representation; 2], ETS.metric, None);
            let entries = (0..16)
                .map(|flat| {
                    Atom::num(if flat / 4 != flat % 4 {
                        0
                    } else if lorentzian && flat != 0 {
                        -1
                    } else {
                        1
                    })
                })
                .collect();
            library.insert_explicit(key.map_canonical(|structure| {
                MixedTensor::Param(ParamTensor::param(
                    DenseTensor::from_storage_data(entries, structure)
                        .unwrap()
                        .into(),
                ))
            }));
        }
        let vectors = ["k", "p"].map(|label| {
            let name = SPENSO_TAG.rank_one_tensor_symbol(&format!("os_algebra_reference::{label}"));
            let entries = (0..4)
                .map(|axis| {
                    Atom::var(symbol!(&format!("os_algebra_reference::{label}{axis}"); Scalar))
                })
                .collect();
            let key = ExplicitKey::from_iter([mink.to_lib()], name, None);
            library.insert_explicit(key.map_canonical(|structure| {
                MixedTensor::Param(ParamTensor::param(
                    DenseTensor::from_storage_data(entries, structure)
                        .unwrap()
                        .into(),
                ))
            }));
            function!(name, mink.to_symbolic([]))
        });
        Self {
            library,
            spin: ["a", "b", "c", "d"].map(|name| {
                bis.slot::<AbstractIndex, _>(symbol!(&format!("os_algebra_reference::{name}")))
                    .to_atom()
            }),
            lorentz: ["mu", "nu"].map(|name| {
                mink.slot::<AbstractIndex, _>(symbol!(&format!("os_algebra_reference::{name}")))
                    .to_atom()
            }),
            vectors,
        }
    }

    fn basis(component: usize, slot: &Atom) -> Atom {
        function!(ETS.delta, function!(AIND_SYMBOLS.cind, component), slot)
    }

    fn projector(&self, velocity: &[Atom; 4], sign: i32, left: usize, right: usize) -> Atom {
        let [row, column] = [&self.spin[left], &self.spin[right]];
        let spatial = Atom::add_many(
            (1..4).map(|axis| &velocity[axis] * Self::basis(axis, &self.lorentz[0])),
        );
        let slash = &velocity[0] * function!(AGS.gamma0, row, column)
            + spatial * function!(AGS.gamma, row, column, &self.lorentz[0]);
        (function!(ETS.metric, row, column) + Atom::num(sign) * slash) / Atom::num(2)
    }

    fn reference(velocity: &[Atom; 4], sign: i32) -> Matrix {
        let [v0, v1, v2, v3] = velocity;
        // Independent Weyl blocks, including the imaginary sigma2 entries.
        let slash = [
            [Atom::Zero, Atom::Zero, v0 - v3, -v1 + Atom::i() * v2],
            [Atom::Zero, Atom::Zero, -v1 - Atom::i() * v2, v0 + v3],
            [v0 + v3, v1 - Atom::i() * v2, Atom::Zero, Atom::Zero],
            [v1 + Atom::i() * v2, v0 - v3, Atom::Zero, Atom::Zero],
        ];
        std::array::from_fn(|row| {
            std::array::from_fn(|column| {
                (Atom::num(i32::from(row == column)) + Atom::num(sign) * &slash[row][column])
                    / Atom::num(2)
            })
        })
    }

    fn component(&self, expression: &Atom, row: usize, column: usize) -> Atom {
        // Select entries with basis vectors, retaining abstract gamma slots
        // until the independent component engine executes the matrices.
        let expression =
            expression * Self::basis(row, &self.spin[0]) * Self::basis(column, &self.spin[2]);
        let settings =
            ParseSettings::default().with_strict_tensor_filter(StrictTensorFilter::ContainsReps);
        let mut network =
            HepNet::<AbstractIndex>::try_from_view(expression.as_view(), &self.library, &settings)
                .unwrap();
        network
            .execute::<Sequential, SmallestDegree, _, _, _>(&self.library, &*FUN_LIB)
            .unwrap();
        match network.result_scalar().unwrap() {
            ExecutionResult::One => Atom::one(),
            ExecutionResult::Zero => Atom::Zero,
            ExecutionResult::Val(value) => value.into_owned(),
        }
    }

    fn simplify(expression: &Atom, gamma: GammaSimplifySettings) -> Atom {
        SymbolicTensor::infer(expression.clone())
            .unwrap()
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(gamma),
                ..Default::default()
            })
            .unwrap()
            .into_expression()
    }

    fn velocities() -> [[Atom; 4]; 3] {
        [
            [Atom::one(), Atom::Zero, Atom::Zero, Atom::Zero],
            [Atom::num((5, 4)), Atom::num((3, 4)), Atom::Zero, Atom::Zero],
            [
                Atom::num((5, 4)),
                Atom::num((1, 4)),
                Atom::num((1, 2)),
                Atom::num((1, 2)),
            ],
        ]
    }
}

#[test]
fn boosted_projector_is_left_multiplication_in_explicit_weyl_components() {
    let evaluation = ProjectorEvaluation::new();
    for velocity in ProjectorEvaluation::velocities() {
        for sign in [1, -1] {
            let projector = evaluation.projector(&velocity, sign, 0, 2);
            let inserted = evaluation.projector(&velocity, sign, 0, 1)
                * ProjectorEvaluation::basis(0, &evaluation.spin[1])
                * ProjectorEvaluation::basis(2, &evaluation.spin[2]);
            let reference = ProjectorEvaluation::reference(&velocity, sign);
            for gamma in [
                GammaSimplifySettings::canonical(),
                GammaSimplifySettings::repeated_pairs(),
            ] {
                let reduced = ProjectorEvaluation::simplify(&projector, gamma);
                let reduced_inserted = ProjectorEvaluation::simplify(&inserted, gamma);
                let mut differs_from_right = false;
                for row in 0..4 {
                    for column in 0..4 {
                        for expression in [&projector, &reduced] {
                            assert_eq!(
                                evaluation.component(expression, row, column),
                                reference[row][column]
                            );
                        }
                        let left = if column == 2 {
                            reference[row][0].clone()
                        } else {
                            Atom::Zero
                        };
                        let right = if row == 0 {
                            reference[2][column].clone()
                        } else {
                            Atom::Zero
                        };
                        for expression in [&inserted, &reduced_inserted] {
                            assert_eq!(evaluation.component(expression, row, column), left);
                        }
                        differs_from_right |= left != right;
                    }
                }
                assert!(differs_from_right, "E02 must distinguish P*E02 from E02*P");
            }
        }
    }
}

#[test]
fn boosted_projectors_obey_exact_weyl_projector_algebra() {
    let evaluation = ProjectorEvaluation::new();
    for velocity in ProjectorEvaluation::velocities() {
        let projectors = [1, -1].map(|sign| evaluation.projector(&velocity, sign, 0, 2));
        for row in 0..4 {
            for column in 0..4 {
                assert_eq!(
                    evaluation.component(&(&projectors[0] + &projectors[1]), row, column),
                    Atom::num(i32::from(row == column))
                );
            }
        }
        for left in [1, -1] {
            for right in [1, -1] {
                let first = evaluation.projector(&velocity, left, 0, 1);
                // Distinct Lorentz dummy slots are essential when both
                // completed projectors contain their own vector contraction.
                let second = evaluation
                    .projector(&velocity, right, 1, 2)
                    .replace(evaluation.lorentz[0].to_pattern())
                    .with(evaluation.lorentz[1].to_pattern());
                let product = first * second;
                let reference = ProjectorEvaluation::reference(&velocity, left);
                for gamma in [
                    GammaSimplifySettings::canonical(),
                    GammaSimplifySettings::repeated_pairs(),
                ] {
                    let reduced = ProjectorEvaluation::simplify(&product, gamma);
                    for (row, entries) in reference.iter().enumerate() {
                        for (column, entry) in entries.iter().enumerate() {
                            let expected = if left == right {
                                entry.clone()
                            } else {
                                Atom::Zero
                            };
                            for expression in [&product, &reduced] {
                                assert_eq!(evaluation.component(expression, row, column), expected);
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn symbolic_spin_reduction_preserves_os_phase_order_and_scalar_factors() {
    let evaluation = ProjectorEvaluation::new();
    let [a, b, c, d] = &evaluation.spin;
    let mu = &evaluation.lorentz[0];
    let mass = Atom::var(symbol!("os_algebra_reference::mass"; Scalar));
    let self_energy = Atom::i()
        * function!(AGS.gamma, a, b, mu)
        * (function!(AGS.gamma, b, d, &evaluation.vectors[0])
            + &mass * function!(ETS.metric, b, d))
        * function!(AGS.gamma, d, c, mu);
    // Use a fresh input row for the insertion, without renaming an existing
    // internal spinor contraction onto another one.
    let fresh = Bispinor {}
        .new_rep(4)
        .slot::<AbstractIndex, _>(symbol!("os_algebra_reference::fresh"))
        .to_atom();
    let inserted = function!(AGS.gamma0, a, &fresh)
        * self_energy.replace(a.to_pattern()).with(fresh.to_pattern());
    let spectator = (Atom::var(symbol!("os_algebra_reference::x"; Scalar))
        + Atom::var(symbol!("os_algebra_reference::y"; Scalar)))
        * (Atom::var(symbol!("os_algebra_reference::z"; Scalar))
            + Atom::var(symbol!("os_algebra_reference::w"; Scalar)))
        .pow(11);
    let dot = spenso::dot!(&evaluation.vectors[0], &evaluation.vectors[1]);
    for original in [self_energy, inserted] {
        for gamma in [
            GammaSimplifySettings::canonical(),
            GammaSimplifySettings::repeated_pairs(),
        ] {
            let reduced = ProjectorEvaluation::simplify(&original, gamma);
            for power in [0, 2, 3] {
                let factor = &spectator * dot.pow(power);
                assert_eq!(
                    ProjectorEvaluation::simplify(&(&factor * &original), gamma),
                    &factor * &reduced,
                    "scalar spectators and dot power {power} must remain opaque"
                );
            }
            for row in 0..4 {
                for column in 0..4 {
                    // The eight vector entries and mass remain independent
                    // symbols. Normalize only this small finite Clifford
                    // coefficient; the literal scalar spectators above are
                    // never distributed, even in a diagnostic copy.
                    let actual = evaluation.component(&original, row, column);
                    let expected = evaluation.component(&reduced, row, column);
                    assert!(
                        (actual - expected).expand().is_zero(),
                        "symbolic phase/order identity, entry ({row},{column})"
                    );
                }
            }
        }
    }
}
