//! Generic scalar DOT input for the saved-application example.
//!
//! Graph parsing, routing, mass substitution and physical-slot handling remain
//! in the existing FeynKit/Vakint APIs. This adapter does not select a topology.

use std::{error::Error, path::PathBuf};

use feynkit_graph::{FeynmanDiagram, IntegralFamily};
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use symbolica::{
    atom::{Atom, AtomCore},
    id::Replacement,
    prelude::Integer,
};
use vakint::{VakintExpression, diagram_integral::DiagramIntegralOptions};

pub type Result<T> = std::result::Result<T, Box<dyn Error>>;

#[derive(Debug, PartialEq, Eq)]
pub struct ScalarDotInput {
    pub path: PathBuf,
    pub powers: Vec<i64>,
}

impl ScalarDotInput {
    /// Parse only the explicit DOT mode; leave existing named fixtures intact.
    pub fn from_arguments(arguments: &[PathBuf]) -> Result<Option<Self>> {
        if !arguments
            .iter()
            .skip(1)
            .any(|arg| arg.as_os_str() == "--dot")
        {
            return Ok(None);
        }
        let mut path = None;
        let mut powers = None;
        let mut arguments = arguments.iter().skip(1);
        while let Some(option) = arguments.next() {
            let value = arguments.next().ok_or("missing DOT option value")?;
            match option.to_str() {
                Some("--dot") if path.is_none() => path = Some(value.clone()),
                Some("--powers") if powers.is_none() => {
                    let text = value.to_str().ok_or("powers must be UTF-8 integers")?;
                    powers = Some(
                        text.split(',')
                            .map(|value| value.trim().parse::<i64>())
                            .collect::<std::result::Result<Vec<_>, _>>()?,
                    );
                }
                _ => return Err("DOT mode accepts exactly one --dot PATH and --powers CSV".into()),
            }
        }
        Ok(Some(Self {
            path: path.ok_or("DOT input requires --dot PATH")?,
            powers: powers.ok_or("DOT input requires --powers CSV")?,
        }))
    }

    pub fn prepare(&self) -> Result<PreparedScalarDot> {
        Self::prepare_text(&std::fs::read_to_string(&self.path)?, &self.powers)
    }

    /// The example's phi3/phi4 scalar model supplies the mass symbol. Normalize
    /// its sole mass to one consistently in both the native family and ingress.
    pub fn prepare_text(dot: &str, powers: &[i64]) -> Result<PreparedScalarDot> {
        let model = Model::phi_3_4();
        let diagram = FeynmanDiagram::from_dot(model.clone(), dot)?;
        let original = diagram.propagator_family(&Kinematics::new())?;
        if original.denominators().len() != powers.len() {
            return Err(format!(
                "expected {} powers in native family order, received {}",
                original.denominators().len(),
                powers.len()
            )
            .into());
        }
        let mass = model
            .particle_by_id(model.particle_id("phi")?)?
            .symbolic_mass(&model);
        let substitutions = vec![(mass, Atom::one())];
        let replacements = substitutions
            .iter()
            .map(|(source, target)| Replacement::new(source.clone(), target.clone()))
            .collect::<Vec<_>>();
        let family = IntegralFamily::new(
            original.loop_momenta().to_vec(),
            original.external_momenta().to_vec(),
            original
                .denominators()
                .iter()
                .map(|value| value.replace_multiple(&replacements))
                .collect(),
            original.kinematics(),
        )?;
        let expression = Atom::from(VakintExpression::from_diagram(
            &diagram,
            &family,
            &Atom::one(),
            &DiagramIntegralOptions {
                powers: Some(powers.iter().copied().map(Integer::from).collect()),
                parameter_substitutions: substitutions,
                ..Default::default()
            },
        )?);
        Ok(PreparedScalarDot {
            diagram,
            family,
            expression,
        })
    }
}

pub struct PreparedScalarDot {
    pub diagram: FeynmanDiagram,
    pub family: IntegralFamily,
    pub expression: Atom,
}
