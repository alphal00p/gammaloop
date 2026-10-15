//! Finite physical numerator inputs, never reduction-engine topology dispatch.
use std::collections::BTreeSet;

use feynkit_graph::{DiagramEdge, DiagramVertex, FeynmanDiagram, IntegralFamily};
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use rustred::family::IntegralKey;
use symbolica::{
    atom::{Atom, AtomCore},
    id::Replacement,
    prelude::Integer,
};

use crate::{VakintExpression, diagram_integral::DiagramIntegralOptions, symbols::S, vakint_parse};

pub(super) type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

pub(super) struct Fixture {
    pub name: &'static str,
    pub dot: String,
    pub input: Atom,
    pub pieces: Vec<(i64, Atom)>,
    pub common_numerator: Atom,
    pub base: IntegralKey,
    pub expected: BTreeSet<IntegralKey>,
    pub native_to_common: Vec<Atom>,
    pub scalar_degree: usize,
    pub momentum_degree: usize,
}

fn diagram(banana: bool) -> Result<FeynmanDiagram> {
    let model = Model::phi_3_4();
    let particle = model.particle_id("phi")?;
    let mut builder = FeynmanDiagram::builder(model.clone(), "physical_r2_control");
    let left = builder.add_vertex(DiagramVertex {
        name: "left".into(),
        interaction: None,
        numerator: Atom::one(),
    });
    let right = if banana {
        builder.add_vertex(DiagramVertex {
            name: "right".into(),
            interaction: None,
            numerator: Atom::one(),
        })
    } else {
        left
    };
    let mut edges = Vec::new();
    for _ in 0..if banana { 6 } else { 5 } {
        edges.push(builder.add_edge(left, right, DiagramEdge::new(particle, false))?);
    }
    // Declare the fixture basis explicitly through the native checked service;
    // do not depend on the default spanning-tree selection or graph names.
    let dot = builder
        .build()?
        .with_loop_momentum_edges(&edges[..5])?
        .to_dot()?;
    Ok(FeynmanDiagram::from_dot(model, &dot)?)
}

fn family(diagram: &FeynmanDiagram) -> Result<(IntegralFamily, Atom)> {
    // Explicitly D-dimensional, unlike Kinematics::new()'s four-dimensional
    // default. Scalar projection later substitutes the physical d=4-2 eps.
    let kinematics = Kinematics::in_dimension(&vakint_parse!("d")?)?;
    let original = diagram.propagator_family(&kinematics)?;
    let model = diagram.model();
    let mass = model
        .particle_by_id(model.particle_id("phi")?)?
        .symbolic_mass(model);
    let replacement = [Replacement::new(mass.clone(), Atom::one())];
    let family = IntegralFamily::new(
        original.loop_momenta().to_vec(),
        vec![],
        original
            .denominators()
            .iter()
            .map(|value| value.replace_multiple(&replacement))
            .collect(),
        original.kinematics(),
    )?;
    assert_eq!(diagram.loop_count(), 5);
    diagram.loop_momentum_basis().validate(diagram)?;
    for index in 0..5 {
        let square = family
            .kinematics()
            .scalar_product(&family.loop_momenta()[index], &family.loop_momenta()[index])?;
        assert!(
            (&family.denominators()[index] - (square - Atom::one()))
                .expand()
                .together()
                .is_zero()
        );
    }
    Ok((family, mass))
}

fn input(
    diagram: &FeynmanDiagram,
    family: &IntegralFamily,
    mass: &Atom,
    numerator: &Atom,
    powers: &[i64],
) -> Result<Atom> {
    Ok(Atom::from(VakintExpression::from_diagram(
        diagram,
        family,
        numerator,
        &DiagramIntegralOptions {
            powers: Some(powers.iter().copied().map(Integer::from).collect()),
            parameter_substitutions: vec![(mass.clone(), Atom::one())],
            ..Default::default()
        },
    )?))
}

fn shifted(base: &[i64; 15], slots: &[usize]) -> IntegralKey {
    let mut powers = *base;
    for &slot in slots {
        powers[slot] -= 1;
    }
    IntegralKey::try_new(powers).unwrap()
}

pub(super) fn homogeneous_factor() -> Result<Fixture> {
    let diagram = diagram(false)?;
    let (family, mass) = family(&diagram)?;
    let loops = family.loop_momenta();
    let numerator = family.kinematics().scalar_product(&loops[0], &loops[0])?
        - family.kinematics().scalar_product(&loops[1], &loops[1])?;
    let powers = [2, 3, 1, 2, 3];
    let base = [0, 0, 0, 0, 2, 0, 0, 3, 1, 0, 0, 2, 0, 0, 3];
    let mut first = powers;
    first[0] -= 1;
    let mut second = powers;
    second[1] -= 1;
    Ok(Fixture {
        name: "homogeneous_factor_difference",
        dot: diagram.to_dot()?,
        input: input(&diagram, &family, &mass, &numerator, &powers)?,
        pieces: vec![
            (1, input(&diagram, &family, &mass, &Atom::one(), &first)?),
            (-1, input(&diagram, &family, &mass, &Atom::one(), &second)?),
        ],
        common_numerator: vakint_parse!("dot(k(5),k(5))-dot(k(1)-k(5),k(1)-k(5))")?,
        base: IntegralKey::try_new(base)?,
        expected: [shifted(&base, &[4]), shifted(&base, &[7])].into(),
        native_to_common: factor_coordinates()?,
        scalar_degree: 1,
        momentum_degree: 2,
    })
}

pub(super) fn cross_factor() -> Result<Fixture> {
    let diagram = diagram(false)?;
    let (family, mass) = family(&diagram)?;
    let loops = family.loop_momenta();
    let numerator = family
        .kinematics()
        .scalar_product(&loops[0], &(&loops[1] - &loops[3]))?
        .pow(2);
    let powers = [2, 3, 1, 3, 3];
    let base = [0, 0, 0, 0, 2, 0, 0, 3, 1, 0, 0, 3, 0, 0, 3];
    let slots = [0, 2, 7, 11];
    let expected = (0..4)
        .flat_map(|i| (i..4).map(move |j| shifted(&base, &[slots[i], slots[j]])))
        .collect();
    Ok(Fixture {
        name: "nonzero_cross_factor_contraction",
        dot: diagram.to_dot()?,
        input: input(&diagram, &family, &mass, &numerator, &powers)?,
        pieces: vec![],
        common_numerator: vakint_parse!("dot(k(5),k(1)-k(3))^2")?,
        base: IntegralKey::try_new(base)?,
        expected,
        native_to_common: factor_coordinates()?,
        scalar_degree: 2,
        momentum_degree: 4,
    })
}

pub(super) fn connected() -> Result<Fixture> {
    let diagram = diagram(true)?;
    let (family, mass) = family(&diagram)?;
    let loops = family.loop_momenta();
    let first = &loops[1] + &loops[3] + &loops[0];
    let second = &loops[3] + &loops[0];
    let p1 = family.kinematics().scalar_product(&first, &first)? - Atom::one();
    let p3 = family.kinematics().scalar_product(&second, &second)? - Atom::one();
    let numerator = &p1 - &p3;
    let powers = [1, 1, 1, 1, 1, 5];
    let base = [0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 5];
    Ok(Fixture {
        name: "connected_irreducible_difference",
        dot: diagram.to_dot()?,
        input: input(&diagram, &family, &mass, &numerator, &powers)?,
        pieces: vec![
            (1, input(&diagram, &family, &mass, &p1, &powers)?),
            (-1, input(&diagram, &family, &mass, &p3, &powers)?),
        ],
        common_numerator: vakint_parse!("dot(k(1),k(1))-dot(k(3),k(3))")?,
        base: IntegralKey::try_new(base)?,
        expected: [shifted(&base, &[0]), shifted(&base, &[2])].into(),
        native_to_common: [
            "k(5)",
            "k(1)-k(3)",
            "k(2)-k(3)",
            "k(3)-k(5)",
            "-k(1)-k(2)+k(4)",
        ]
        .iter()
        .map(|text| vakint_parse!(text).map_err(Into::into))
        .collect::<Result<_>>()?,
        scalar_degree: 1,
        momentum_degree: 2,
    })
}

fn factor_coordinates() -> Result<Vec<Atom>> {
    ["k(5)", "k(1)-k(5)", "k(2)-k(3)", "k(3)-k(5)", "k(3)-k(4)"]
        .iter()
        .map(|text| vakint_parse!(text).map_err(Into::into))
        .collect()
}

impl Fixture {
    pub fn verify_native_family_map(&self) -> Result<()> {
        let diagram = FeynmanDiagram::from_dot(Model::phi_3_4(), &self.dot)?;
        let (family, _) = family(&diagram)?;
        let common = super::cancellation_fixtures::common_family(
            family.loop_momenta(),
            family.kinematics(),
            &Atom::one(),
        )?;
        let native_labels = family
            .loop_momenta()
            .iter()
            .enumerate()
            .map(|(i, target)| Replacement::new(S.k.call(i + 1), target.clone()))
            .collect::<Vec<_>>();
        let coordinates = self
            .native_to_common
            .iter()
            .map(|value| value.replace_multiple(&native_labels))
            .collect::<Vec<_>>();
        let witness = family
            .mapping_to(&common, &coordinates)?
            .ok_or("fixture common-family map lacks a native all-row/unit-Jacobian witness")?;
        let expected = self
            .base
            .powers()
            .iter()
            .enumerate()
            .filter_map(|(slot, power)| (*power > 0).then_some(slot))
            .collect::<Vec<_>>();
        assert_eq!(witness.denominator_map(), expected);
        Ok(())
    }

    pub fn transported_input(&self) -> Atom {
        let replacements = self
            .native_to_common
            .iter()
            .enumerate()
            .map(|(i, target)| Replacement::new(S.k.call(i + 1), target.clone()))
            .collect::<Vec<_>>();
        self.input.replace_multiple(&replacements)
    }

    pub fn transported_numerator(&self) -> Result<Atom> {
        let transported = self.transported_input();
        let expression = VakintExpression::try_from(transported.as_view())?;
        assert_eq!(expression.0.len(), 1);
        Ok(expression.0[0].numerator.clone())
    }
}
