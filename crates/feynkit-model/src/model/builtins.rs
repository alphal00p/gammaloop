use std::collections::{BTreeMap, BTreeSet};

use miniz_oxide::inflate::decompress_to_vec_zlib;

use super::{
    CouplingDefinition as Coupling, LorentzStructureDefinition as LorentzStructure, Model,
    ModelDefinition, Order, ParameterDefinition as Parameter, ParameterNature, ParameterType,
    ParticleDefinition as Particle, PropagatorDefinition as Propagator,
    VertexRuleDefinition as VertexRule,
};
use crate::ComplexValue;

impl Model {
    /// Load a fresh copy of the complete embedded Standard Model (`sm`).
    ///
    /// Includes 43 particles and 153 vertices with their default parameter
    /// values. The compressed data is decoded on demand, without file access
    /// or a Python/UFO runtime.
    pub fn standard_model() -> Self {
        let json = decompress_to_vec_zlib(include_bytes!("../../data/sm.json.zlib"))
            .expect("embedded Standard Model must decompress");
        let definition =
            serde_json::from_slice(&json).expect("embedded Standard Model must contain valid JSON");
        Self::new(definition).expect("embedded Standard Model must be valid")
    }

    /// The Standard Model's QCD sector: six quark flavors, gluons, and gluon
    /// ghosts. Keeps the Standard Model parameters and their default values.
    pub fn qcd() -> Self {
        Self::standard_model_sector("qcd", super::Particle::is_qcd_charged)
    }

    /// The Standard Model's QED sector: photons and all charged fermions,
    /// including quarks. Keeps the Standard Model parameters and their defaults.
    pub fn qed() -> Self {
        Self::standard_model_sector("qed", |particle| {
            particle.pdg_code == 22 || (particle.is_fermion() && !particle.charge.is_zero())
        })
    }

    /// The electroweak Standard Model, including quarks, Higgs, Goldstones,
    /// and electroweak ghosts, with gluons and gluon ghosts removed.
    pub fn electroweak() -> Self {
        Self::standard_model_sector("electroweak", |particle| {
            !particle.is_qcd_charged() || particle.is_fermion()
        })
    }

    /// Strong and electromagnetic interactions of all quarks and charged
    /// leptons, including gluon ghosts, with Standard Model parameters.
    pub fn qcd_qed() -> Self {
        Self::standard_model_sector("qcd_qed", |particle| {
            particle.is_qcd_charged()
                || particle.pdg_code == 22
                || (particle.is_fermion() && !particle.charge.is_zero())
        })
    }

    /// Pure SU(3) Yang–Mills theory: gluons and gluon ghosts, without quarks.
    /// Preserves the Standard Model parameters and their default values.
    pub fn yang_mills() -> Self {
        Self::standard_model_sector("yang_mills", |particle| {
            particle.is_qcd_charged() && !particle.is_fermion()
        })
    }

    fn standard_model_sector(name: &str, keep_particle: impl Fn(&super::Particle) -> bool) -> Self {
        let model = Self::standard_model();
        let selected: BTreeSet<_> = model
            .particles()
            .iter()
            .filter(|p| keep_particle(p))
            .map(|p| p.name.as_str())
            .collect();
        let mut definition = model.definition();
        definition.name = name.to_owned();
        definition
            .particles
            .retain(|p| selected.contains(p.name.as_str()));
        let pdgs: BTreeSet<_> = definition.particles.iter().map(|p| p.pdg_code).collect();
        definition
            .covariant_cut_multiplets
            .retain(|_, members| members.iter().all(|pdg| pdgs.contains(pdg)));
        let particles: BTreeSet<_> = definition.particles.iter().map(|p| &p.name).collect();
        definition
            .propagators
            .retain(|propagator| particles.contains(&propagator.particle));
        definition
            .vertex_rules
            .retain(|vertex| vertex.particles.iter().all(|p| particles.contains(p)));

        let couplings: BTreeSet<_> = definition
            .vertex_rules
            .iter()
            .flat_map(|vertex| vertex.couplings.iter().flatten().flatten())
            .collect();
        definition
            .couplings
            .retain(|coupling| couplings.contains(&coupling.name));
        let lorentz_structures: BTreeSet<_> = definition
            .vertex_rules
            .iter()
            .flat_map(|vertex| &vertex.lorentz_structures)
            .collect();
        definition
            .lorentz_structures
            .retain(|lorentz| lorentz_structures.contains(&lorentz.name));
        let orders: BTreeSet<_> = definition
            .couplings
            .iter()
            .flat_map(|coupling| coupling.orders.keys())
            .collect();
        definition
            .orders
            .retain(|order| orders.contains(&order.name));

        // Keep the complete parameter definitions: expressions may depend on
        // electroweak inputs even after electroweak particles have been removed.
        Self::new(definition).expect("filtered Standard Model sector must be valid")
    }

    /// One real scalar `phi` with L_int = -g phi^3 / 3!.
    /// External parameters `mass` and `g` default to one; the width is zero.
    pub fn phi3() -> Self {
        Self::real_scalar("phi3", 3, "g")
    }

    /// One real scalar `phi` with L_int = -lam phi^4 / 4!.
    /// External parameters `mass` and `lam` default to one; the width is zero.
    pub fn phi4() -> Self {
        Self::real_scalar("phi4", 4, "lam")
    }

    /// One real scalar `phi` with L_int = -g phi^3 / 3! - lam phi^4 / 4!.
    /// External parameters `mass`, `g`, and `lam` default to one; the width is zero.
    pub fn phi_3_4() -> Self {
        let mut definition = Self::phi3().definition();
        let mut quartic = Self::phi4().definition();
        definition.name = "phi_3_4".to_owned();
        let mut coupling_parameter = quartic.parameters.pop().unwrap();
        coupling_parameter.lhacode = Some(vec![3]);
        definition.parameters.push(coupling_parameter);

        quartic.lorentz_structures[0].name = "SCALAR4".to_owned();
        quartic.couplings[0].name = "SCALAR4_COUPLING".to_owned();
        quartic.vertex_rules[0].name = "V_SCALAR4".to_owned();
        quartic.vertex_rules[0].lorentz_structures = vec!["SCALAR4".to_owned()];
        quartic.vertex_rules[0].couplings = vec![vec![Some("SCALAR4_COUPLING".to_owned())]];
        definition
            .lorentz_structures
            .extend(quartic.lorentz_structures);
        definition.couplings.extend(quartic.couplings);
        definition.vertex_rules.extend(quartic.vertex_rules);
        Self::new(definition).expect("built-in cubic and quartic scalar model must be valid")
    }

    fn real_scalar(name: &str, valence: usize, coupling: &str) -> Self {
        let parameters = [("ZERO", 0.0), ("mass", 1.0), (coupling, 1.0)]
            .into_iter()
            .enumerate()
            .map(|(index, (name, value))| Parameter {
                name: name.to_owned(),
                lhablock: (index != 0).then(|| "SCALAR".to_owned()),
                lhacode: (index != 0).then(|| vec![index]),
                nature: if index == 0 {
                    ParameterNature::Internal
                } else {
                    ParameterNature::External
                },
                parameter_type: ParameterType::Real,
                value: Some(ComplexValue::new(value, 0.0)),
                expression: None,
            })
            .collect();
        let definition = ModelDefinition {
            covariant_cut_multiplets: Default::default(),
            name: name.to_owned(),
            restriction: None,
            orders: vec![Order {
                name: "SCALAR".to_owned(),
                expansion_order: 99,
                hierarchy: 1,
            }],
            parameters,
            particles: vec![Particle {
                pdg_code: 9000001,
                name: "phi".to_owned(),
                antiname: "phi".to_owned(),
                spin: 1,
                color: 1,
                mass: "mass".to_owned(),
                width: "ZERO".to_owned(),
                texname: "\\phi".to_owned(),
                antitexname: "\\phi".to_owned(),
                charge: 0.into(),
                ghost_number: 0,
                lepton_number: 0,
                y_charge: None,
                y_charge_right: None,
                propagating: true,
                goldstone: false,
                propagator: Some("phi_prop".to_owned()),
            }],
            propagators: vec![Propagator {
                name: "phi_prop".to_owned(),
                particle: "phi".to_owned(),
                numerator: "1𝑖".to_owned(),
                denominator: "(UFO::P(UFO::idx(1,1)))^2-UFO::mass^2".to_owned(),
            }],
            lorentz_structures: vec![LorentzStructure {
                name: "SCALAR".to_owned(),
                spins: vec![1; valence],
                structure: "1".to_owned(),
            }],
            couplings: vec![Coupling {
                name: "SCALAR_COUPLING".to_owned(),
                expression: format!("-1𝑖*UFO::{coupling}"),
                orders: BTreeMap::from([("SCALAR".to_owned(), 1)]),
                value: Some(ComplexValue::new(0.0, -1.0)),
            }],
            vertex_rules: vec![VertexRule {
                name: "V_SCALAR".to_owned(),
                particles: vec!["phi".to_owned(); valence],
                color_structures: vec!["1".to_owned()],
                lorentz_structures: vec!["SCALAR".to_owned()],
                couplings: vec![vec![Some("SCALAR_COUPLING".to_owned())]],
            }],
            functions: vec![],
            form_factors: vec![],
        };
        Self::new(definition).expect("built-in real scalar model must be valid")
    }

    /// Scalar electrodynamics in Feynman gauge, with `a`, `phi+`, and `phi-`.
    /// `mass` and `e` default to one, and `lam` defaults to zero. The scalar
    /// potential is mass^2 |phi|^2 + lam |phi|^4 / 4; all widths vanish.
    pub fn scalar_qed() -> Self {
        let mut definition = Self::phi4().definition();
        definition.name = "scalar_qed".to_owned();
        let quartic = &mut definition.parameters[2];
        quartic.value = Some(ComplexValue::new(0.0, 0.0));
        let mut charge = quartic.clone();
        charge.name = "e".to_owned();
        charge.lhacode = Some(vec![3]);
        charge.value = Some(ComplexValue::new(1.0, 0.0));
        definition.parameters.push(charge);
        definition.couplings[0].value = Some(ComplexValue::new(0.0, 0.0));
        definition.vertex_rules[0].particles =
            ["phi-", "phi-", "phi+", "phi+"].map(str::to_owned).to_vec();

        let scalar = definition.particles.pop().unwrap();
        let propagator = definition.propagators.pop().unwrap();
        for (name, antiname, sign) in [("phi+", "phi-", 1), ("phi-", "phi+", -1)] {
            let propagator_name = format!("{name}_prop");
            definition.particles.push(Particle {
                name: name.to_owned(),
                antiname: antiname.to_owned(),
                pdg_code: sign * scalar.pdg_code,
                charge: sign.into(),
                texname: format!("\\phi^{{{}}}", if sign > 0 { "+" } else { "-" }),
                antitexname: format!("\\phi^{{{}}}", if sign > 0 { "-" } else { "+" }),
                propagator: Some(propagator_name.clone()),
                ..scalar.clone()
            });
            definition.propagators.push(Propagator {
                name: propagator_name,
                particle: name.to_owned(),
                ..propagator.clone()
            });
        }
        definition.particles.push(Particle {
            pdg_code: 22,
            name: "a".to_owned(),
            antiname: "a".to_owned(),
            spin: 3,
            mass: "ZERO".to_owned(),
            texname: "\\gamma".to_owned(),
            antitexname: "\\gamma".to_owned(),
            propagator: Some("a_prop".to_owned()),
            ..scalar
        });
        definition.propagators.push(Propagator {
            name: "a_prop".to_owned(),
            particle: "a".to_owned(),
            numerator: "-1𝑖*UFO::Metric(UFO::idx(1,1),UFO::idx(1,2))".to_owned(),
            denominator: "(UFO::P(UFO::idx(1,1)))^2".to_owned(),
        });
        definition.orders.push(Order {
            name: "QED".to_owned(),
            expansion_order: 99,
            hierarchy: 1,
        });
        // All momenta are incoming. Match the SM's a G- G+ convention: the
        // cubic rule is -i e (p_- - p_+) and the seagull is 2 i e^2 g(mu,nu).
        for (name, particles, structure, expression, power, value) in [
            (
                "VSS",
                vec!["a", "phi-", "phi+"],
                "-1*UFO::P(UFO::idx(1,1),UFO::idx(1,3))+UFO::P(UFO::idx(1,1),UFO::idx(1,2))",
                "-1𝑖*UFO::e",
                1,
                -1.0,
            ),
            (
                "VVSS",
                vec!["a", "a", "phi-", "phi+"],
                "UFO::Metric(UFO::idx(1,1),UFO::idx(1,2))",
                "2𝑖*UFO::e^2",
                2,
                2.0,
            ),
        ] {
            definition.lorentz_structures.push(LorentzStructure {
                name: name.to_owned(),
                spins: particles
                    .iter()
                    .map(|p| if *p == "a" { 3 } else { 1 })
                    .collect(),
                structure: structure.to_owned(),
            });
            definition.couplings.push(Coupling {
                name: name.to_owned(),
                expression: expression.to_owned(),
                orders: BTreeMap::from([("QED".to_owned(), power)]),
                value: Some(ComplexValue::new(0.0, value)),
            });
            definition.vertex_rules.push(VertexRule {
                name: name.to_owned(),
                particles: particles.into_iter().map(str::to_owned).collect(),
                color_structures: vec!["1".to_owned()],
                lorentz_structures: vec![name.to_owned()],
                couplings: vec![vec![Some(name.to_owned())]],
            });
        }
        Self::new(definition).expect("built-in scalar QED model must be valid")
    }
}

#[cfg(test)]
mod tests {
    use super::Model;
    use symbolica::{atom::Atom, domains::rational::Rational};

    #[test]
    fn standard_model_quantum_numbers_and_antiparticles_are_exact() {
        let model = Model::standard_model();
        for (pdgs, charge, yl, yr, t3) in [
            (
                &[2, 4, 6][..],
                Rational::from((2, 3)),
                Rational::from((1, 3)),
                Some(Rational::from((4, 3))),
                Rational::from((1, 2)),
            ),
            (
                &[1, 3, 5][..],
                Rational::from((-1, 3)),
                Rational::from((1, 3)),
                Some(Rational::from((-2, 3))),
                Rational::from((-1, 2)),
            ),
            (
                &[11, 13, 15][..],
                Rational::from(-1),
                Rational::from(-1),
                Some(Rational::from(-2)),
                Rational::from((-1, 2)),
            ),
            (
                &[12, 14, 16][..],
                Rational::from(0),
                Rational::from(-1),
                None,
                Rational::from((1, 2)),
            ),
        ] {
            for &pdg in pdgs {
                let p = model.particle_by_pdg(pdg).unwrap();
                let anti = model.particle_by_id(p.antiparticle).unwrap();
                assert_eq!(p.charge, charge);
                assert_eq!(p.y_charge, Some(yl.clone()));
                assert_eq!(p.y_charge_right, yr);
                assert_eq!(p.weak_isospin(), Some(t3.clone()));
                assert_eq!(
                    p.weak_isospin_right(),
                    yr.as_ref().map(|_| Rational::zero())
                );
                assert_eq!(anti.charge, -charge.clone());
                assert_eq!(anti.y_charge, yr.as_ref().map(|y| -y.clone()));
                assert_eq!(anti.y_charge_right, Some(-yl.clone()));
                assert_eq!(anti.weak_isospin_right(), Some(-t3.clone()));
                assert_eq!(p.mass, anti.mass);
                assert_eq!(p.width, anti.width);
                assert_eq!(p.spin, anti.spin);
                assert_eq!(p.lepton_number, -anti.lepton_number);
                assert_eq!(
                    p.color,
                    if p.color == 1 {
                        anti.color
                    } else {
                        -anti.color
                    }
                );
            }
        }
        for name in ["H", "G0"] {
            assert!(model.particle(name).unwrap().weak_isospin().is_none());
        }
        for (name, sign) in [("G+", 1), ("G-", -1)] {
            let p = model.particle(name).unwrap();
            assert_eq!(p.y_charge, Some(sign.into()));
            assert_eq!(p.weak_isospin(), Some(Rational::from((sign, 2))));
            assert!(p.is_goldstone());
        }
        assert!(
            model
                .particle("ve")
                .unwrap()
                .symbolic_mass(&model)
                .is_zero()
        );
        assert_eq!(
            model.particle("c").unwrap().symbolic_mass(&model),
            Atom::var(symbolica::symbol!("UFO::MC"))
        );
        for vertex in model.vertex_rules() {
            let particles: Vec<_> = vertex
                .particles
                .iter()
                .map(|&id| model.particle_by_id(id).unwrap())
                .collect();
            assert!(
                particles
                    .iter()
                    .fold(Rational::zero(), |q, p| q + &p.charge)
                    .is_zero(),
                "{} violates charge conservation",
                vertex.name
            );
            assert_eq!(
                particles.iter().map(|p| p.ghost_number).sum::<i64>(),
                0,
                "{}",
                vertex.name
            );
            assert_eq!(
                particles.iter().map(|p| p.lepton_number).sum::<i64>(),
                0,
                "{}",
                vertex.name
            );
        }
    }

    #[test]
    fn quantum_number_json_preserves_exactness_and_unknowns() {
        let model = Model::standard_model();
        let mut json: serde_json::Value = serde_json::from_str(&model.to_json().unwrap()).unwrap();
        let i = json["particles"]
            .as_array()
            .unwrap()
            .iter()
            .position(|p| p["name"] == "c")
            .unwrap();
        assert_eq!(json["particles"][i]["charge"], "2/3");
        for invalid in [
            serde_json::json!(0.6666666666666666),
            serde_json::json!(1e30),
            serde_json::json!("1/0"),
            serde_json::json!("not a charge"),
        ] {
            json["particles"][i]["charge"] = invalid;
            assert!(Model::from_json(&json.to_string()).is_err());
        }
        json["particles"][i]["charge"] =
            serde_json::json!("123456789012345678901234567891/999999999999999999999999999999");
        json["particles"][i]
            .as_object_mut()
            .unwrap()
            .remove("y_charge");
        json["particles"][i]
            .as_object_mut()
            .unwrap()
            .remove("y_charge_right");
        let model = Model::from_json(&json.to_string()).unwrap();
        let p = model.particle("c").unwrap();
        assert_eq!(
            p.charge.to_string(),
            "123456789012345678901234567891/999999999999999999999999999999"
        );
        assert!(p.weak_isospin().is_none());
        assert!(p.weak_isospin_right().is_none());
        let reloaded = Model::from_json(&model.to_json().unwrap()).unwrap();
        assert_eq!(model.fingerprint(), reloaded.fingerprint());
    }

    #[test]
    fn builtin_models_validate_and_round_trip() {
        for (model, name, particles, vertices) in [
            (Model::standard_model(), "sm", 43, 153),
            (Model::qcd(), "qcd", 15, 9),
            (Model::qed(), "qed", 19, 9),
            (Model::electroweak(), "electroweak", 40, 144),
            (Model::qcd_qed(), "qcd_qed", 22, 18),
            (Model::yang_mills(), "yang_mills", 3, 3),
            (Model::phi3(), "phi3", 1, 1),
            (Model::phi4(), "phi4", 1, 1),
            (Model::phi_3_4(), "phi_3_4", 1, 2),
            (Model::scalar_qed(), "scalar_qed", 3, 3),
        ] {
            assert_eq!(model.name(), name);
            assert_eq!(model.particles().len(), particles);
            assert_eq!(model.vertex_rules().len(), vertices);
            let reloaded = Model::from_json(&model.to_json().unwrap()).unwrap();
            assert_eq!(reloaded.fingerprint(), model.fingerprint());
        }
    }

    #[test]
    fn compressed_standard_model_matches_the_source_fixture() {
        let fixture = Model::from_json(include_str!("../../tests/fixtures/sm.json")).unwrap();
        assert_eq!(Model::standard_model().fingerprint(), fixture.fingerprint());
        let disk =
            Model::from_json(include_str!("../../../../assets/models/json/sm/sm.json")).unwrap();
        assert_eq!(disk.fingerprint(), fixture.fingerprint());
        assert!(include_bytes!("../../data/sm.json.zlib").len() < 8 * 1024);
    }
}
