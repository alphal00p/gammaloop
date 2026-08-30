use symbolica::atom::{Atom, AtomCore};
use symbolica::function;

use crate::symbols::S;
use crate::utils::vakint_macros::vk_symbol;
use crate::{ReplacementRules, Topology, VakintSettings, get_integer_from_atom, get_prop_with_id};

use super::artifact::ArtifactFamily;
use super::{RustRedEvaluationError, RustRedEvaluationOptions};

pub(super) struct MatchedScalarFamily {
    pub(super) family: ArtifactFamily,
    pub(super) powers: Vec<i64>,
    pub(super) loop_count: usize,
    pub(super) mass_squared: Atom,
}

impl ArtifactFamily {
    pub(super) fn from_topology(topology: &Topology) -> Result<Self, RustRedEvaluationError> {
        let integral = topology.get_integral();
        let expression = integral.canonical_expression.as_ref().ok_or_else(|| {
            RustRedEvaluationError::UnsupportedMatchedFamily {
                detail: "the matcher did not provide a canonical expression".to_owned(),
            }
        })?;
        match topology {
            Topology::OneLoop(_) if integral.n_loops == 1 && integral.n_props == 1 => {
                require_momentum(expression, 1, function!(S.k, Atom::num(1)))?;
                Ok(Self::UnitMassVacuumK1)
            }
            Topology::TwoLoop(_) if integral.n_loops == 2 && integral.n_props == 3 => {
                require_momentum(expression, 1, function!(S.k, Atom::num(1)))?;
                require_momentum(expression, 2, function!(S.k, Atom::num(2)))?;
                if get_prop_with_id(expression.as_view(), 3).is_some() {
                    require_momentum(
                        expression,
                        3,
                        function!(S.k, Atom::num(1)) + function!(S.k, Atom::num(2)),
                    )?;
                }
                Ok(Self::UnitMassVacuumK3)
            }
            _ => Err(RustRedEvaluationError::UnsupportedMatchedFamily {
                detail: format!(
                    "matcher class has {} loops and {} parent propagators",
                    integral.n_loops, integral.n_props
                ),
            }),
        }
    }
}

impl MatchedScalarFamily {
    pub(super) fn try_new(
        integral_specs: &ReplacementRules,
    ) -> Result<Self, RustRedEvaluationError> {
        Self::try_from_topology(&integral_specs.canonical_topology)
    }

    pub(super) fn try_from_topology(topology: &Topology) -> Result<Self, RustRedEvaluationError> {
        let family = ArtifactFamily::from_topology(topology)?;
        let integral = topology.get_integral();
        let expression = integral
            .canonical_expression
            .as_ref()
            .expect("artifact admission requires a canonical expression");

        let mut powers = Vec::with_capacity(integral.n_props);
        let mut common_mass_squared = None;
        for propagator in 1..=integral.n_props {
            let Some(properties) = get_prop_with_id(expression.as_view(), propagator) else {
                powers.push(0);
                continue;
            };
            let power = properties
                .get(&vk_symbol!("pow_"))
                .and_then(|power| get_integer_from_atom(power.as_view()))
                .ok_or_else(|| RustRedEvaluationError::InvalidPower {
                    propagator,
                    power: properties
                        .get(&vk_symbol!("pow_"))
                        .map_or_else(|| "<missing>".to_owned(), Atom::to_canonical_string),
                })?;
            powers.push(power);

            let propagator_mass = properties
                .get(&vk_symbol!("mUVsq_"))
                .expect("Vakint's canonical propagator matcher captures the mass");
            if let Some(expected) = &common_mass_squared {
                if propagator_mass != expected {
                    return Err(RustRedEvaluationError::InvalidMatchedFamily {
                        detail: "physical propagators do not share one exact mass squared"
                            .to_owned(),
                    });
                }
            } else {
                common_mass_squared = Some(propagator_mass.clone());
            }
        }
        let mass_squared =
            common_mass_squared.ok_or_else(|| RustRedEvaluationError::InvalidMatchedFamily {
                detail: "the matched family has no physical common mass squared".to_owned(),
            })?;
        if mass_squared.as_view().is_zero() {
            return Err(RustRedEvaluationError::InvalidMatchedFamily {
                detail: "the common mass squared is exactly zero".to_owned(),
            });
        }
        family.validate_root_powers(&powers)?;

        Ok(Self {
            family,
            powers,
            loop_count: integral.n_loops,
            mass_squared,
        })
    }

    pub(super) fn evaluate(
        self,
        settings: &VakintSettings,
        numerator: symbolica::atom::AtomView,
        options: &RustRedEvaluationOptions,
    ) -> Result<Atom, RustRedEvaluationError> {
        super::materialize::evaluate(self, settings, numerator, options)
    }
}

fn require_momentum(
    expression: &Atom,
    propagator: usize,
    expected: Atom,
) -> Result<(), RustRedEvaluationError> {
    let properties = get_prop_with_id(expression.as_view(), propagator).ok_or_else(|| {
        RustRedEvaluationError::InvalidMatchedFamily {
            detail: format!("canonical propagator {propagator} is absent"),
        }
    })?;
    let actual = properties
        .get(&vk_symbol!("q_"))
        .expect("Vakint's canonical propagator matcher captures momentum");
    if actual != &expected {
        return Err(RustRedEvaluationError::InvalidMatchedFamily {
            detail: format!(
                "canonical propagator {propagator} has momentum {}, expected {}",
                actual.to_canonical_string(),
                expected.to_canonical_string()
            ),
        });
    }
    Ok(())
}
