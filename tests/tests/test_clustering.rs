use gammalooprs::momentum::{Energy, FourMomentum, ThreeMomentum};
use gammalooprs::observables::{JetAlgorithm, JetClustering};
use gammalooprs::utils::F;
use rand::{Rng, SeedableRng, rngs::SmallRng};
use std::cmp::Ordering;
use std::ffi::c_void;

#[cfg(target_vendor = "apple")]
#[link(name = "c++")]
unsafe extern "C" {}

#[cfg(all(not(target_vendor = "apple"), not(target_env = "msvc")))]
#[link(name = "stdc++")]
unsafe extern "C" {}

#[link(name = "fjcore_test", kind = "static")]
unsafe extern "C" {
    fn fastjet_workspace() -> *mut c_void;
    fn fastjet_free(workspace: *mut c_void);
    fn fastjetppgenkt_(
        workspace: *mut c_void,
        p: *const f64,
        npart: *const i32,
        r: *const f64,
        ptjet_min: *const f64,
        palg: *const f64,
        jets: *mut f64,
        njets: *mut i32,
        whichjets: *mut i32,
    );
}

struct FjcoreWorkspace(*mut c_void);

impl FjcoreWorkspace {
    fn new() -> Self {
        let workspace = unsafe { fastjet_workspace() };
        assert!(!workspace.is_null(), "fjcore workspace allocation failed");
        Self(workspace)
    }
}

impl Drop for FjcoreWorkspace {
    fn drop(&mut self) {
        unsafe { fastjet_free(self.0) };
    }
}

#[derive(Debug, Clone)]
struct FjcoreJet {
    momentum: [f64; 4],
    constituents: Vec<usize>,
}

impl FjcoreJet {
    fn pt(&self) -> f64 {
        (self.momentum[1] * self.momentum[1] + self.momentum[2] * self.momentum[2]).sqrt()
    }

    fn phi(&self) -> f64 {
        let mut phi = self.momentum[2].atan2(self.momentum[1]);
        if phi < 0.0 {
            phi += std::f64::consts::TAU;
        }
        phi
    }

    fn rapidity(&self) -> f64 {
        let px = self.momentum[1];
        let py = self.momentum[2];
        let pz = self.momentum[3];
        let energy = self.momentum[0];
        let pt2 = px * px + py * py;
        if pt2 == 0.0 && energy == pz.abs() {
            let max_rap_here = 1.0e5 + pz.abs();
            if pz >= 0.0 {
                max_rap_here
            } else {
                -max_rap_here
            }
        } else {
            let mass2 = (energy * energy - px * px - py * py - pz * pz).max(0.0);
            let e_plus_abs_pz = energy + pz.abs();
            let mut rapidity = 0.5 * ((pt2 + mass2) / (e_plus_abs_pz * e_plus_abs_pz)).ln();
            if pz > 0.0 {
                rapidity = -rapidity;
            }
            rapidity
        }
    }
}

fn palg(algorithm: JetAlgorithm) -> f64 {
    match algorithm {
        JetAlgorithm::Kt => 1.0,
        JetAlgorithm::CambridgeAachen => 0.0,
        JetAlgorithm::AntiKt => -1.0,
    }
}

fn massless_momentum(pt: f64, rapidity: f64, phi: f64) -> FourMomentum<F<f64>> {
    let px = pt * phi.cos();
    let py = pt * phi.sin();
    let pz = pt * rapidity.sinh();
    let energy = pt * rapidity.cosh();
    FourMomentum {
        temporal: Energy { value: F(energy) },
        spatial: ThreeMomentum {
            px: F(px),
            py: F(py),
            pz: F(pz),
        },
    }
}

fn compare_jets(lhs: &FjcoreJet, rhs: &FjcoreJet) -> Ordering {
    rhs.pt()
        .partial_cmp(&lhs.pt())
        .unwrap_or(Ordering::Equal)
        .then_with(|| {
            lhs.rapidity()
                .partial_cmp(&rhs.rapidity())
                .unwrap_or(Ordering::Equal)
        })
        .then_with(|| lhs.phi().partial_cmp(&rhs.phi()).unwrap_or(Ordering::Equal))
        .then_with(|| lhs.constituents.cmp(&rhs.constituents))
}

fn run_fjcore(
    algorithm: JetAlgorithm,
    r: f64,
    min_jpt: f64,
    momenta: &[FourMomentum<F<f64>>],
) -> Vec<FjcoreJet> {
    let workspace = FjcoreWorkspace::new();
    let npart = i32::try_from(momenta.len()).expect("number of particles should fit into i32");
    let mut input = Vec::with_capacity(momenta.len() * 4);
    for momentum in momenta {
        input.push(momentum.temporal.value.0);
        input.push(momentum.spatial.px.0);
        input.push(momentum.spatial.py.0);
        input.push(momentum.spatial.pz.0);
    }

    let mut jets_out = vec![0.0; momenta.len() * 4];
    let mut which_jets = vec![0_i32; momenta.len()];
    let mut njets = 0_i32;
    let palg = palg(algorithm);

    unsafe {
        fastjetppgenkt_(
            workspace.0,
            input.as_ptr(),
            &npart,
            &r,
            &min_jpt,
            &palg,
            jets_out.as_mut_ptr(),
            &mut njets,
            which_jets.as_mut_ptr(),
        );
    }

    let njets = usize::try_from(njets).expect("number of jets should be non-negative");
    let mut constituents = vec![Vec::new(); njets];
    for (particle_index, jet_index) in which_jets.into_iter().enumerate() {
        if jet_index > 0 {
            constituents[(jet_index - 1) as usize].push(particle_index);
        }
    }

    let mut jets = Vec::with_capacity(njets);
    for (jet_index, jet_constituents) in constituents.into_iter().enumerate() {
        let start = 4 * jet_index;
        jets.push(FjcoreJet {
            momentum: [
                jets_out[start],
                jets_out[start + 1],
                jets_out[start + 2],
                jets_out[start + 3],
            ],
            constituents: jet_constituents,
        });
    }
    jets.sort_unstable_by(compare_jets);
    jets
}

fn assert_close(lhs: f64, rhs: f64, tolerance: f64, label: &str) {
    let scale = lhs.abs().max(rhs.abs()).max(1.0);
    let delta = (lhs - rhs).abs();
    assert!(
        delta <= tolerance * scale,
        "{label} mismatch: lhs={lhs:.16e}, rhs={rhs:.16e}, delta={delta:.16e}, tolerance={tolerance:.16e}"
    );
}

fn compare_against_fjcore(
    algorithm: JetAlgorithm,
    r: f64,
    min_jpt: f64,
    momenta: &[FourMomentum<F<f64>>],
) {
    let rust_clustering = JetClustering::new(algorithm, r, min_jpt, Vec::new());
    let rust_result = rust_clustering.cluster_momenta(momenta);
    let fjcore_result = run_fjcore(algorithm, r, min_jpt, momenta);

    assert_eq!(
        rust_result.jets.len(),
        fjcore_result.len(),
        "jet multiplicity mismatch for {:?}",
        algorithm
    );

    for (rust_jet, fjcore_jet) in rust_result.jets.iter().zip(fjcore_result.iter()) {
        assert_eq!(
            rust_jet.constituent_indices(),
            fjcore_jet.constituents.as_slice(),
            "constituent mismatch for {:?}",
            algorithm
        );
        assert_close(
            rust_jet.momentum.temporal.value.0,
            fjcore_jet.momentum[0],
            1.0e-11,
            "energy",
        );
        assert_close(
            rust_jet.momentum.spatial.px.0,
            fjcore_jet.momentum[1],
            1.0e-11,
            "px",
        );
        assert_close(
            rust_jet.momentum.spatial.py.0,
            fjcore_jet.momentum[2],
            1.0e-11,
            "py",
        );
        assert_close(
            rust_jet.momentum.spatial.pz.0,
            fjcore_jet.momentum[3],
            1.0e-11,
            "pz",
        );
    }
}

#[test]
fn fjcore_regression_fixed_cases() {
    let cases = [
        vec![
            massless_momentum(80.0, 0.2, 0.1),
            massless_momentum(25.0, 0.25, 0.18),
            massless_momentum(35.0, -1.1, 2.4),
        ],
        vec![
            massless_momentum(45.0, -0.7, 5.7),
            massless_momentum(30.0, -0.68, 5.8),
            massless_momentum(12.0, 1.4, 1.3),
            massless_momentum(18.0, 1.35, 1.36),
        ],
        vec![
            massless_momentum(60.0, -2.2, 0.9),
            massless_momentum(14.0, -2.05, 1.0),
            massless_momentum(22.0, 0.4, 3.0),
            massless_momentum(11.0, 0.45, 3.08),
            massless_momentum(8.0, 1.9, 4.7),
        ],
    ];

    for algorithm in [
        JetAlgorithm::Kt,
        JetAlgorithm::CambridgeAachen,
        JetAlgorithm::AntiKt,
    ] {
        for case in &cases {
            compare_against_fjcore(algorithm, 0.6, 5.0, case);
            compare_against_fjcore(algorithm, 0.4, 15.0, case);
        }
    }
}

#[test]
fn fjcore_regression_random_cases() {
    let mut rng = SmallRng::seed_from_u64(0x5eed_cafe_u64);

    for _ in 0..24 {
        let npart = rng.random_range(1..=7);
        let r = rng.random_range(0.35..0.9);
        let min_jpt = rng.random_range(0.0..12.0);
        let mut momenta = Vec::with_capacity(npart);
        for _ in 0..npart {
            let pt = rng.random_range(0.1..120.0);
            let rapidity = rng.random_range(-3.0..3.0);
            let phi = rng.random_range(0.0..std::f64::consts::TAU);
            momenta.push(massless_momentum(pt, rapidity, phi));
        }

        for algorithm in [
            JetAlgorithm::Kt,
            JetAlgorithm::CambridgeAachen,
            JetAlgorithm::AntiKt,
        ] {
            compare_against_fjcore(algorithm, r, min_jpt, &momenta);
        }
    }
}
