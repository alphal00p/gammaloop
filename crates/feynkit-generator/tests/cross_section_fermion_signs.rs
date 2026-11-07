//! Physical spin sums pin cross-section external-fermion signs.
//!
//! Direct Dirac traces for a massive lepton decay and its crossed W decay
//! provide absolute positive rational results, independent of the generator's
//! cycle-count bookkeeping. The complete covariant W/Goldstone final-state
//! sector and its Yukawa relation are kept in the lepton decay.

use std::iter;
use std::sync::Arc;

use feynkit_amplitude::SpinSum;
use feynkit_generator::{GenerationOptions, GenerationType, Process};
use feynkit_graph::expressions::evaluate_overall_factor;
use feynkit_graph::{FeynmanDiagram, symbols};
use feynkit_model::Model;
use idenso::dirac::GammaSimplifySettings;
use idenso::tensor::{AlgebraSettings, ContractSettings, SymbolicTensor};
use spenso::network::library::symbolic::ETS;
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::dimension::Dimension;
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::{function, symbol};

/// Two-body decay of `initial` at rest. Masses are exact integers, so every
/// invariant is rational: `2 p_i.p_j` is stored in `doubled_gram`.
struct Decay {
    initial: &'static str,
    finals: [&'static str; 2],
    masses: [i64; 3],
}

impl Decay {
    fn doubled_gram(&self) -> [[i64; 3]; 3] {
        assert!(self.masses[0] > self.masses[1] + self.masses[2]);
        let [m0, m1, m2] = self.masses.map(|mass| mass * mass);
        [
            [2 * m0, m0 + m1 - m2, m0 + m2 - m1],
            [m0 + m1 - m2, 2 * m1, m0 - m1 - m2],
            [m0 + m2 - m1, m0 - m1 - m2, 2 * m2],
        ]
    }

    /// The sum over forward graphs of the overall factor times the numerator,
    /// with physical spin sums on the sewn initial state, evaluated at the
    /// physical decay kinematics.
    fn spin_summed_forward_numerator(&self, model: &Arc<Model>) -> Atom {
        let generated = Process::new([self.initial], self.finals)
            .generate_cross_section(
                model.clone(),
                &GenerationOptions::default()
                    .with_loop_count(1, 1)
                    .unwrap()
                    .threads(1)
                    .max_vertices(2),
            )
            .unwrap();
        assert!(!generated.diagrams.is_empty());
        let value = generated
            .diagrams
            .iter()
            .map(|diagram| self.forward_numerator(model, diagram))
            .fold(Atom::Zero, |sum, term| sum + term)
            .expand();
        assert!(
            matches!(value.as_view(), AtomView::Num(_)) && !value.is_zero(),
            "{} > {:?} must give a nonzero number, found {value}",
            self.initial,
            self.finals
        );
        value
    }

    fn forward_numerator(&self, model: &Arc<Model>, diagram: &FeynmanDiagram) -> Atom {
        // Assign each edge its physical momentum (0: initial state, 1 and 2:
        // final states) and solve the orientation signs from momentum
        // conservation. A global flip is irrelevant: after the traces every
        // term has an even number of momenta.
        let representatives = Process::new([self.initial], self.finals)
            .covariant_cut_representatives(model, GenerationType::CrossSection)
            .unwrap();
        let final_pdgs = self.finals.map(|name| {
            model
                .particle_by_id(model.particle_id(name).unwrap())
                .unwrap()
                .pdg_code
        });
        let mut edges = Vec::new();
        for (edge, endpoints, data) in diagram.edges() {
            let particle = model.particle_by_id(data.particle).unwrap();
            let antiparticle = model.particle_by_id(particle.antiparticle).unwrap();
            let physical = if data.external.is_some() {
                0
            } else {
                1 + final_pdgs
                    .iter()
                    .position(|pdg| {
                        *pdg == *representatives
                            .get(&particle.pdg_code)
                            .unwrap_or(&particle.pdg_code)
                            || *pdg
                                == *representatives
                                    .get(&antiparticle.pdg_code)
                                    .unwrap_or(&antiparticle.pdg_code)
                    })
                    .unwrap_or_else(|| {
                        panic!(
                            "{} is not a final state of {}",
                            particle.name,
                            diagram.name()
                        )
                    })
            };
            edges.push((edge.0, physical, endpoints));
        }
        let gram = self.doubled_gram();
        let vertices = edges
            .iter()
            .flat_map(|(_, _, endpoints)| [endpoints.source, endpoints.target])
            .flatten()
            .collect::<Vec<_>>();
        let orientations = (0..1_u32 << edges.len())
            .map(|bits| {
                (0..edges.len())
                    .map(|position| if bits >> position & 1 == 1 { -1 } else { 1 })
                    .collect::<Vec<i64>>()
            })
            .find(|signs| {
                vertices.iter().all(|vertex| {
                    // The inflow at each vertex must be orthogonal to the
                    // spanning momenta, hence zero.
                    (0..3).all(|probe| {
                        edges
                            .iter()
                            .zip(signs)
                            .map(|((_, physical, endpoints), sign)| {
                                let inflow = i64::from(endpoints.target == Some(*vertex))
                                    - i64::from(endpoints.source == Some(*vertex));
                                inflow * sign * gram[*physical][probe]
                            })
                            .sum::<i64>()
                            == 0
                    })
                })
            })
            .expect("edge momenta must be consistent with the decay kinematics");

        let mut projector = diagram.projector().expression().clone();
        for (edge, endpoints, data) in diagram.edges() {
            if data.external.is_none() || endpoints.source.is_none() || endpoints.target.is_none() {
                continue;
            }
            let particle = model.particle_by_id(data.particle).unwrap();
            projector = SpinSum::new(particle, model)
                .unwrap()
                .apply(
                    &projector,
                    &symbols::momentum().call(edge.0),
                    edge.0,
                    None,
                    None,
                )
                .unwrap();
        }
        let factor = evaluate_overall_factor(diagram.overall_factor().as_view());
        let integrand = diagram
            .numerator()
            .multiply(diagram.numerator_prefactor())
            .unwrap()
            .multiply(&SymbolicTensor::infer(projector).unwrap())
            .unwrap()
            .with_lorentz_dimension(Dimension::Concrete(4))
            .unwrap();
        let integrand = integrand
            .with_checked_expression(integrand.expression() * &factor)
            .unwrap()
            .contract(ContractSettings::default().without_rank_one_tensors())
            .unwrap()
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::default()),
                epsilon: true,
                ..Default::default()
            })
            .unwrap()
            .contract(ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap()
            .into_expression()
            .expand();

        let mut value = integrand.clone();
        for ((left, left_physical, _), left_sign) in edges.iter().zip(&orientations) {
            for ((right, right_physical, _), right_sign) in edges.iter().zip(&orientations) {
                let dot = function!(
                    SPENSO_TAG.dot,
                    function!(symbols::momentum(), *left, symbol!("sign_test::left___")),
                    function!(symbols::momentum(), *right, symbol!("sign_test::right___"))
                );
                let invariant = Atom::num((
                    left_sign * right_sign * gram[*left_physical][*right_physical],
                    2,
                ));
                value = value.replace(dot.to_pattern()).with(invariant.to_pattern());
                let metric = function!(
                    ETS.metric,
                    function!(symbols::momentum(), *left, symbol!("sign_test::left___")),
                    function!(symbols::momentum(), *right, symbol!("sign_test::right___"))
                );
                value = value
                    .replace(metric.to_pattern())
                    .with(invariant.to_pattern());
            }
        }
        for (name, mass) in iter::once(self.initial).chain(self.finals).zip(self.masses) {
            let particle = model
                .particle_by_id(model.particle_id(name).unwrap())
                .unwrap();
            let symbolic = particle.symbolic_mass(model);
            if !symbolic.is_zero() {
                value = value
                    .replace(symbolic.to_pattern())
                    .with(Atom::num(mass).to_pattern());
            } else {
                assert_eq!(mass, 0, "{name} is massless in the model");
            }
        }
        // Keep the vector/Goldstone Ward relation at each process's masses.
        // Choose e/(sqrt(2)*sw) = 1, hence y_tau = m_tau/M_W.
        let [initial_mass, first_mass, second_mass] = self.masses;
        let (tau_mass, w_mass) = if self.initial.starts_with("ta") {
            (initial_mass, second_mass)
        } else {
            (first_mass, initial_mass)
        };
        for (name, coupling) in [
            ("GC_40", Atom::i()),
            ("GC_96", Atom::num((-tau_mass, w_mass))),
            ("GC_97", Atom::num((tau_mass, w_mass))),
        ] {
            value = value
                .replace(Atom::var(symbol!(&format!("UFO::{name}"))).to_pattern())
                .with(coupling.to_pattern());
        }
        let value = value.expand();
        eprintln!(
            "{} > {:?}: {} particles {:?} factor {factor} => {value}",
            self.initial,
            self.finals,
            diagram.name(),
            diagram
                .edges()
                .map(|(_, _, d)| model.particle_by_id(d.particle).unwrap().name.as_str())
                .collect::<Vec<_>>()
        );
        assert!(
            matches!(value.as_view(), AtomView::Num(_)),
            "{} did not reduce to a number:\n{value}\nfrom\n{integrand}",
            diagram.name()
        );
        value
    }
}

#[test]
fn incoming_fermion_and_crossed_decays_match_positive_physical_spin_sums() {
    let model = Arc::new(
        Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap(),
    );
    for (fermion, boson) in [
        (
            Decay {
                initial: "ta-",
                finals: ["vt", "W-"],
                masses: [173, 0, 80],
            },
            Decay {
                initial: "W-",
                finals: ["ta-", "vt~"],
                masses: [200, 20, 0],
            },
        ),
        (
            Decay {
                initial: "ta+",
                finals: ["vt~", "W+"],
                masses: [173, 0, 80],
            },
            Decay {
                initial: "W+",
                finals: ["ta+", "vt"],
                masses: [200, 20, 0],
            },
        ),
    ] {
        let fermion_value = fermion.spin_summed_forward_numerator(&model);
        let boson_value = boson.spin_summed_forward_numerator(&model);
        // For c = e/(sqrt(2)*sw) = 1, a direct chiral Dirac trace gives
        // (m^2-M^2)(m^2+2M^2)/M^2 for l -> nu W, and
        // (M^2-m^2)(2M^2+m^2)/M^2 for W -> l nu. Neither includes
        // initial-spin averaging. At the two rational mass points above:
        let expected_fermion = Atom::num((1_005_370_641, 6_400));
        let expected_boson = Atom::num(79_596);
        assert_eq!(
            fermion_value, expected_fermion,
            "{} > {:?} must match its positive physical spin sum",
            fermion.initial, fermion.finals
        );
        assert_eq!(
            boson_value, expected_boson,
            "{} > {:?} must match its positive physical spin sum",
            boson.initial, boson.finals
        );
    }
}
