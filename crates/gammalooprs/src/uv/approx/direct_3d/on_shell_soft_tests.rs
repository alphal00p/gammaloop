//! Additional high-precision radial controls of complete local forest sums.
//! Exact operator, hard/soft Taylor and scalar residue certificates live in
//! local_4d::on_shell_forest_tests; these sampled slopes do not replace them.

use super::on_shell_tests::DirectForestFixture;
use crate::utils::{ArbPrec, F};
use color_eyre::Result;
use symbolica::prelude::*;

#[derive(Clone, Copy, Debug)]
enum SoftStratum {
    PhysicalChildGluon,
    ChildOnShellGluon,
    OuterGluon,
    SimultaneousGluons,
    DisjointLeftGluon,
    DisjointRightGluon,
    DisjointGluons,
    ThreeLoopOuterOnShellGluon,
    ThreeLoopMiddleOnShellGluon,
    ThreeLoopInnerOnShellGluon,
    ThreeLoopFlat {
        label: &'static str,
        centers: [[i32; 3]; 3],
        normal_slots: [Option<usize>; 3],
    },
}

impl SoftStratum {
    fn momenta(self, radius: &Atom, ray: usize) -> Vec<[Atom; 3]> {
        let direction = if ray == 0 { [1, 0, 0] } else { [0, 1, 0] };
        let transverse = if ray == 0 { [0, 1, 0] } else { [1, 0, 0] };
        let soft = std::array::from_fn::<_, 3, _>(|index| Atom::num(direction[index]) * radius);
        let hard_parent = [Atom::num(5), Atom::Zero, Atom::Zero];
        let hard_child = [Atom::Zero, Atom::num(9), Atom::Zero];
        match self {
            Self::PhysicalChildGluon => vec![
                hard_parent.clone(),
                std::array::from_fn(|index| &hard_parent[index] + &soft[index]),
            ],
            Self::ChildOnShellGluon | Self::DisjointRightGluon => vec![hard_parent, soft],
            Self::OuterGluon | Self::DisjointLeftGluon => vec![soft, hard_child],
            Self::SimultaneousGluons => vec![
                soft.clone(),
                std::array::from_fn(|index| &soft[index] + Atom::num(transverse[index]) * radius),
            ],
            // Each component owns its own subtraction point. The second
            // soft gluon is measured from zero, without a shift by the first
            // component's loop momentum.
            Self::DisjointGluons => vec![
                soft,
                std::array::from_fn(|index| Atom::num(transverse[index]) * radius),
            ],
            Self::ThreeLoopOuterOnShellGluon => {
                vec![soft, hard_child, [Atom::Zero, Atom::Zero, Atom::num(16)]]
            }
            Self::ThreeLoopMiddleOnShellGluon => {
                vec![hard_parent, soft, [Atom::Zero, Atom::Zero, Atom::num(16)]]
            }
            Self::ThreeLoopInnerOnShellGluon => vec![hard_parent, hard_child, soft],
            Self::ThreeLoopFlat {
                centers,
                normal_slots,
                ..
            } => {
                // Two independent transverse direction sets. Every independent
                // three-vector normal scales uniformly with the same radius.
                let directions = if ray == 0 {
                    [[2, 1, -1], [-1, 3, 2], [1, -2, 3]]
                } else {
                    [[1, 3, 2], [-2, 1, 4], [3, -1, 2]]
                };
                centers
                    .into_iter()
                    .zip(normal_slots)
                    .map(|(center, slot)| {
                        std::array::from_fn(|axis| {
                            let fixed = Atom::num(center[axis]);
                            match slot {
                                Some(slot) => fixed + Atom::num(directions[slot][axis]) * radius,
                                None => fixed,
                            }
                        })
                    })
                    .collect()
            }
        }
    }

    fn radial_measure_power(self) -> i32 {
        match self {
            // Both three-dimensional loop vectors scale together. The radial
            // measure of this six-dimensional stratum is rho^5 d rho.
            Self::SimultaneousGluons | Self::DisjointGluons => 5,
            Self::ThreeLoopFlat { normal_slots, .. } => {
                // Rank r means 3r independent real normal coordinates.
                3 * normal_slots.iter().filter(|slot| slot.is_some()).count() as i32 - 1
            }
            _ => 2,
        }
    }
}

fn check_soft_rays(
    fixture: DirectForestFixture,
    batches: &[([i32; 4], &[SoftStratum])],
) -> Result<()> {
    crate::initialisation::test_initialise()?;
    // Reuse the complete-family exact source certificates and numeric engine.
    // The matched production measurements select explicit direct for the
    // soft-gluon insertion and projected 4D for the other nested fixtures.
    // Small pointwise equivalence tests retain all three representations.
    let route = if matches!(fixture, DirectForestFixture::SoftGluonInQuark) {
        1
    } else {
        2
    };
    let mut evaluator = fixture.full_forest_evaluator(route)?;
    // Numerical external data do not alter the prepared native forest.
    // Reuse its retained families and compiled scalar output probe across
    // the coincident-center and generic-external batches.
    for (external, strata) in batches {
        let external = (*external).map(Atom::num);
        let unit = F(ArbPrec::default()).one();
        let zero = unit.zero();
        for &stratum in *strata {
            for ray in 0..2 {
                let mut weighted = Vec::new();
                for exponent in 3..=7 {
                    let denominator = 1usize << exponent;
                    let radius = Atom::one() / Atom::num(denominator);
                    let values = evaluator.evaluate(&external, &stratum.momenta(&radius, ray))?;
                    let mut norm = zero.clone();
                    for value in values.values() {
                        let component = value.clone().norm().re;
                        assert!(
                            component.is_finite(),
                            "{fixture:?}, {stratum:?}, ray {ray}: finite tensor value"
                        );
                        norm = norm.max(component);
                    }
                    assert!(
                        norm > zero,
                        "{fixture:?}, {stratum:?}, ray {ray}: nonzero full-forest control"
                    );
                    let numeric_radius = &unit / unit.from_usize(denominator);
                    weighted.push(norm * numeric_radius.powi(stratum.radial_measure_power()));
                }
                let flat = match stratum {
                    SoftStratum::ThreeLoopFlat { label, .. } => label,
                    _ => "existing radial stratum",
                };
                crate::debug_tags!(#uv, #local, #soft, #summary;
                    stage = "complete_local_forest_soft_ray",
                    fixture = ?fixture, stratum = ?stratum, flat, route, ray,
                    radial_measure_power = stratum.radial_measure_power(),
                    file.external = ?external, file.weighted = ?weighted,
                    "Complete local CFF radial history at rho=2^-e, e=3..7"
                );
                // On the last two halvings, an uncanceled additional integer pole
                // doubles the expected ratio. The disjoint product is rho^-4,
                // so its rho^5-weighted norm should halve; the other strata only
                // require bounded radial weights. The allowance for subleading
                // terms excludes an additional integer pole in either case.
                // This is a finite-ray control, not a general symbolic IR proof.
                let ratio_denominator = if matches!(stratum, SoftStratum::DisjointGluons) {
                    4
                } else {
                    2
                };
                let maximum_growth = unit.from_usize(3) / unit.from_usize(ratio_denominator);
                for tail in weighted[2..].windows(2) {
                    assert!(
                        tail[1] <= &maximum_growth * &tail[0],
                        "{fixture:?}, {stratum:?}, ray {ray}: radial weight grows toward the soft point: {weighted:?}",
                    );
                }
            }
        }
    }
    Ok(())
}

#[test]
fn complete_os_in_os_forest_has_bounded_radial_soft_controls() -> Result<()> {
    check_soft_rays(
        DirectForestFixture::OsQuarkInQuark,
        &[(
            [5, 0, 0, 0],
            &[
                SoftStratum::PhysicalChildGluon,
                SoftStratum::ChildOnShellGluon,
                SoftStratum::OuterGluon,
                SoftStratum::SimultaneousGluons,
            ],
        )],
    )
}

#[test]
fn complete_os_in_ir_forest_has_bounded_radial_soft_controls() -> Result<()> {
    check_soft_rays(
        DirectForestFixture::OsQuarkInGluon,
        &[(
            [5, 0, 0, 0],
            &[
                SoftStratum::PhysicalChildGluon,
                SoftStratum::ChildOnShellGluon,
            ],
        )],
    )
}

#[test]
fn complete_ir_in_os_forest_has_bounded_radial_soft_controls() -> Result<()> {
    check_soft_rays(
        DirectForestFixture::SoftGluonInQuark,
        &[([5, 0, 0, 0], &[SoftStratum::OuterGluon])],
    )
}

#[test]
fn complete_disjoint_os_forest_has_bounded_radial_soft_controls() -> Result<()> {
    // At the fixed rest-frame OS point each one-loop residue is at most
    // rho^-2. The live bridge has fixed q=(5,0,0,0), m=12, hence no soft
    // pole; the product is at most rho^-4 when both loops scale. Its six-
    // dimensional radial measure is rho^5 d rho. These sampled rays test
    // that bound for the complete (bare + CT)_left * bridge * (bare + CT)_right
    // sum; exact operator and scalar-residue identities are separate tests.
    check_soft_rays(
        DirectForestFixture::DisjointOsQuarks,
        &[(
            [5, 0, 0, 0],
            &[
                SoftStratum::DisjointLeftGluon,
                SoftStratum::DisjointRightGluon,
                SoftStratum::DisjointGluons,
            ],
        )],
    )
}

#[test]
fn complete_three_loop_local_os_forest_has_bounded_uniform_soft_flats() -> Result<()> {
    // Retain the original q=(5,0,0,0) degeneracy control first: the physical
    // outer center coincides with the rest point. The other two carriers
    // stay generic while each of p,k,l separately approaches zero.
    let rest_centers = [
        SoftStratum::ThreeLoopInnerOnShellGluon,
        SoftStratum::ThreeLoopMiddleOnShellGluon,
        SoftStratum::ThreeLoopOuterOnShellGluon,
    ];
    // The factory certifies the physical gluon planes A:p=q, B:k=p, C:l=k
    // and the component-owned rest-point planes D:p=0, E:k=0, F:l=0.
    // A nonzero spatial q separates physical from artificial OS centers.
    // These are the 20 distinct consistent nonempty affine intersections:
    // rank 1: 6, rank 2: 10, rank 3: 4. DE also contains B; EF contains C.
    // This is a uniform finite-ray check of all eight LOCAL forest families.
    // It does not yet include finite integrated prefixes, hierarchical
    // rho/rho^2 approaches, or an all-angle/global integrability claim.
    let q = [1, 2, 3];
    let p = [5, 0, 0];
    let k = [0, 9, 0];
    let l = [0, 0, 16];
    let zero = [0, 0, 0];
    let flats = [
        ("A", [q, k, l], [Some(0), None, None]),
        ("B", [p, p, l], [None, Some(0), None]),
        ("C", [p, k, k], [None, None, Some(0)]),
        ("D", [zero, k, l], [Some(0), None, None]),
        ("E", [p, zero, l], [None, Some(0), None]),
        ("F", [p, k, zero], [None, None, Some(0)]),
        ("AB", [q, q, l], [Some(0), Some(1), None]),
        ("AC", [q, k, k], [Some(0), None, Some(1)]),
        ("BC", [p, p, p], [None, Some(0), Some(1)]),
        ("DE", [zero, zero, l], [Some(0), Some(1), None]),
        ("DF", [zero, k, zero], [Some(0), None, Some(1)]),
        ("EF", [p, zero, zero], [None, Some(0), Some(1)]),
        ("AE", [q, zero, l], [Some(0), Some(1), None]),
        ("AF", [q, k, zero], [Some(0), None, Some(1)]),
        ("BF", [p, p, zero], [None, Some(0), Some(1)]),
        ("CD", [zero, k, k], [Some(0), None, Some(1)]),
        ("ABC", [q, q, q], [Some(0), Some(1), Some(2)]),
        ("DEF", [zero, zero, zero], [Some(0), Some(1), Some(2)]),
        ("ABF", [q, q, zero], [Some(0), Some(1), Some(2)]),
        ("AEF", [q, zero, zero], [Some(0), Some(1), Some(2)]),
    ];
    let strata = flats.map(
        |(label, centers, normal_slots)| SoftStratum::ThreeLoopFlat {
            label,
            centers,
            normal_slots,
        },
    );
    check_soft_rays(
        DirectForestFixture::OsQuarkInQuarkWithSpectator,
        &[
            ([5, 0, 0, 0], rest_centers.as_slice()),
            ([5, 1, 2, 3], strata.as_slice()),
        ],
    )
}
