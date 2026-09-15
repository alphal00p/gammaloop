use std::sync::Arc;

use symbolica::atom::{Atom, AtomCore};
use symbolica::function;
use symbolica::id::Replacement;

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
    pub(super) parent_coordinates: Option<Arc<[Atom]>>,
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
            Topology::ThreeLoop(_) if integral.n_loops == 3 && integral.n_props == 6 => {
                let (parent, coordinates) = integral.parent_routing.as_ref().ok_or_else(|| {
                    RustRedEvaluationError::InvalidMatchedFamily {
                        detail: "the matched three-loop family has no retained parent routing"
                            .to_owned(),
                    }
                })?;
                if coordinates.len() != 3 {
                    return Err(RustRedEvaluationError::InvalidMatchedFamily {
                        detail: "the retained parent routing does not have three coordinates"
                            .to_owned(),
                    });
                }
                let loops = (1..=3)
                    .map(|axis| function!(S.k, Atom::num(axis)))
                    .collect::<Vec<_>>();
                let expected = [
                    loops[0].clone(),
                    loops[1].clone(),
                    loops[2].clone(),
                    &loops[2] - &loops[0],
                    &loops[0] - &loops[1],
                    &loops[1] - &loops[2],
                ];
                let routing = loops
                    .iter()
                    .zip(coordinates.iter())
                    .map(|(source, target)| {
                        Replacement::new(source.to_pattern(), target.to_pattern())
                    })
                    .collect::<Vec<_>>();
                for (axis, momentum) in expected.iter().enumerate() {
                    // Authenticate the defining parent slots, then the actual
                    // contracted routing already supplied by the matcher.
                    // Neither step invokes another graph match or basis solve.
                    require_momentum(parent, axis + 1, momentum.clone())?;
                    if let Some(properties) = get_prop_with_id(expression.as_view(), axis + 1) {
                        let actual = properties
                            .get(&vk_symbol!("q_"))
                            .expect("Vakint's canonical propagator matcher captures momentum");
                        if !(actual.replace_multiple(&routing) - momentum)
                            .expand()
                            .is_zero()
                        {
                            return Err(RustRedEvaluationError::InvalidMatchedFamily {
                                detail: format!(
                                    "canonical propagator {} does not transport to its retained K6 parent slot",
                                    axis + 1
                                ),
                            });
                        }
                    }
                }
                Ok(Self::UnitMassVacuumK6)
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
            parent_coordinates: (family == ArtifactFamily::UnitMassVacuumK6).then(|| {
                integral
                    .parent_routing
                    .as_ref()
                    .expect("K6 admission authenticates parent routing")
                    .1
                    .clone()
            }),
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vakint;
    use crate::utils::vakint_macros::vk_parse;

    #[test]
    fn k6_matching_preserves_physical_slots_mass_and_retained_routing() {
        let vakint = Vakint::new().unwrap();
        let fixtures = [
            ("topo(I3L(13/10,2,1,1,1,1,1))", [2, 1, 1, 1, 1, 1]),
            ("topo(I3L_pinch_6(13/10,2,1,1,1,1,0))", [2, 1, 1, 1, 1, 0]),
            ("topo(I3L_pinch_3_6(13/10,2,1,0,1,1,0))", [2, 1, 0, 1, 1, 0]),
            ("topo(I3L_pinch_1_6(13/10,0,2,1,1,1,0))", [0, 2, 1, 1, 1, 0]),
            (
                "topo(I3L_pinch_1_3_6(13/10,0,2,0,1,1,0))",
                [0, 2, 0, 1, 1, 0],
            ),
        ];
        for (input, powers) in fixtures {
            let mut matched = vakint
                .topologies
                .match_topologies_to_user_input(vk_parse!(input).unwrap().as_view(), false)
                .unwrap()
                .unwrap();
            matched.apply_replacement_rules().unwrap();
            // Engine admission consumes structural matcher metadata, not its name.
            matched.canonical_topology.get_integral_mut().name = "not-a-dispatch-key".into();
            let admitted = MatchedScalarFamily::try_new(&matched).unwrap();
            assert_eq!(admitted.family, ArtifactFamily::UnitMassVacuumK6);
            assert_eq!(admitted.powers, powers);
            assert_eq!(admitted.loop_count, 3);
            assert_eq!(admitted.mass_squared, vk_parse!("13/10").unwrap());
            assert!(Arc::ptr_eq(
                admitted.parent_coordinates.as_ref().unwrap(),
                &matched
                    .canonical_topology
                    .get_integral()
                    .parent_routing
                    .as_ref()
                    .unwrap()
                    .1
            ));
        }
    }

    #[test]
    fn k6_matching_rejects_missing_or_mismatched_parent_evidence() {
        let vakint = Vakint::new().unwrap();
        let mut matched = vakint
            .topologies
            .match_topologies_to_user_input(
                vk_parse!("topo(I3L_pinch_3_6(13/10,2,1,0,1,1,0))")
                    .unwrap()
                    .as_view(),
                false,
            )
            .unwrap()
            .unwrap();
        matched.apply_replacement_rules().unwrap();
        let original = matched.canonical_topology;
        assert_eq!(
            ArtifactFamily::from_topology(&original).unwrap(),
            ArtifactFamily::UnitMassVacuumK6
        );

        let mut absent = original.clone();
        absent.get_integral_mut().parent_routing = None;
        let mut short = original.clone();
        short.get_integral_mut().parent_routing.as_mut().unwrap().1 =
            Arc::from([vk_parse!("k(1)").unwrap()]);
        let mut wrong_map = original.clone();
        let coordinates = &mut wrong_map
            .get_integral_mut()
            .parent_routing
            .as_mut()
            .unwrap()
            .1;
        let mut changed = coordinates.to_vec();
        changed[0] = &changed[0] + vk_parse!("k(1)").unwrap();
        *coordinates = changed.into();
        let mut wrong_parent = original.clone();
        let parent = &mut wrong_parent
            .get_integral_mut()
            .parent_routing
            .as_mut()
            .unwrap()
            .0;
        // Slot 3 is absent in this contraction, but its defining parent binding
        // is still authenticated; a check of present slots alone would miss it.
        assert!(
            get_prop_with_id(
                original
                    .get_integral()
                    .canonical_expression
                    .as_ref()
                    .unwrap()
                    .as_view(),
                3
            )
            .is_none()
        );
        *parent = parent
            .replace(vk_parse!("prop(3,edge(left_,right_),q_,mass_,power_)").unwrap())
            .with(vk_parse!("prop(3,edge(left_,right_),2*q_,mass_,power_)").unwrap());
        assert_ne!(
            parent,
            &original.get_integral().parent_routing.as_ref().unwrap().0
        );
        let mut wrong_actual = original.clone();
        let expression = wrong_actual
            .get_integral_mut()
            .canonical_expression
            .as_mut()
            .unwrap();
        *expression = expression
            .replace(vk_parse!("k(1)").unwrap())
            .with(vk_parse!("2*k(1)").unwrap());
        for invalid in [absent, short, wrong_map, wrong_parent, wrong_actual] {
            assert!(matches!(
                ArtifactFamily::from_topology(&invalid),
                Err(RustRedEvaluationError::InvalidMatchedFamily { .. })
            ));
        }
    }
}
