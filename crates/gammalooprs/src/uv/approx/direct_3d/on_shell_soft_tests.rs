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
}

impl SoftStratum {
    fn momenta(self, radius: &Atom, ray: usize) -> [[Atom; 3]; 2] {
        let direction = if ray == 0 { [1, 0, 0] } else { [0, 1, 0] };
        let transverse = if ray == 0 { [0, 1, 0] } else { [1, 0, 0] };
        let soft = std::array::from_fn::<_, 3, _>(|index| Atom::num(direction[index]) * radius);
        let hard_parent = [Atom::num(5), Atom::Zero, Atom::Zero];
        let hard_child = [Atom::Zero, Atom::num(9), Atom::Zero];
        match self {
            Self::PhysicalChildGluon => [
                hard_parent.clone(),
                std::array::from_fn(|index| &hard_parent[index] + &soft[index]),
            ],
            Self::ChildOnShellGluon | Self::DisjointRightGluon => [hard_parent, soft],
            Self::OuterGluon | Self::DisjointLeftGluon => [soft, hard_child],
            Self::SimultaneousGluons => [
                soft.clone(),
                std::array::from_fn(|index| &soft[index] + Atom::num(transverse[index]) * radius),
            ],
            // Each component owns its own subtraction point. The second
            // soft gluon is measured from zero, without a shift by the first
            // component's loop momentum.
            Self::DisjointGluons => [
                soft,
                std::array::from_fn(|index| Atom::num(transverse[index]) * radius),
            ],
        }
    }

    fn radial_measure_power(self) -> i32 {
        match self {
            // Both three-dimensional loop vectors scale together. The radial
            // measure of this six-dimensional stratum is rho^5 d rho.
            Self::SimultaneousGluons | Self::DisjointGluons => 5,
            _ => 2,
        }
    }
}

fn check_soft_rays(fixture: DirectForestFixture, strata: &[SoftStratum]) -> Result<()> {
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
    let external = [Atom::num(5), Atom::Zero, Atom::Zero, Atom::Zero];
    let unit = F(ArbPrec::default()).one();
    let zero = unit.zero();
    for &stratum in strata {
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
    Ok(())
}

#[test]
fn complete_os_in_os_forest_has_bounded_radial_soft_controls() -> Result<()> {
    check_soft_rays(
        DirectForestFixture::OsQuarkInQuark,
        &[
            SoftStratum::PhysicalChildGluon,
            SoftStratum::ChildOnShellGluon,
            SoftStratum::OuterGluon,
            SoftStratum::SimultaneousGluons,
        ],
    )
}

#[test]
fn complete_os_in_ir_forest_has_bounded_radial_soft_controls() -> Result<()> {
    check_soft_rays(
        DirectForestFixture::OsQuarkInGluon,
        &[
            SoftStratum::PhysicalChildGluon,
            SoftStratum::ChildOnShellGluon,
        ],
    )
}

#[test]
fn complete_ir_in_os_forest_has_bounded_radial_soft_controls() -> Result<()> {
    check_soft_rays(
        DirectForestFixture::SoftGluonInQuark,
        &[SoftStratum::OuterGluon],
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
        &[
            SoftStratum::DisjointLeftGluon,
            SoftStratum::DisjointRightGluon,
            SoftStratum::DisjointGluons,
        ],
    )
}
