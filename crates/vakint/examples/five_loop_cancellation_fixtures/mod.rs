//! Independent positive-line cancellation controls. No saved artifact is needed
//! to prove the numerator/power transport or the structural lowering equality.
use super::five_loop_data as data;

use feynkit_graph::{
    DiagramEdge, DiagramVertex, EdgeId, FeynmanDiagram, IntegralFamily as NativeFamily,
};
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use rustred::algebra::CoefficientContext;
use rustred::family::{AffineDenominator, IntegralFamily, IntegralKey};
use rustred::scalar_numerator::FamilyScalarNumeratorService;
use std::collections::BTreeSet;
use std::error::Error;
use symbolica::atom::{Atom, AtomCore};
use symbolica::id::Replacement;
use symbolica::prelude::Integer;
// Reused from the public example and from the crate's native unit tests.
use super::{DiagramIntegralOptions, Vakint, VakintExpression, vakint_parse, vakint_symbol};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const BASE: [i64; 15] = [0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 5];

pub struct CancellationFixture {
    pub label: String,
    pub cancelled: Atom,
    pub pinched: Atom,
    pub expected_key: IntegralKey,
    pub mass_squared: Atom,
    pub positive_cancellations: Vec<(usize, i64)>,
    /// A separately constructed and reparsed native effective-vertex DOT.
    #[allow(dead_code)] // Used by the executable DOT controls, not every unit-test consumer.
    pub pinched_dot: String,
}

fn parent_endpoints(slot: usize) -> (String, String) {
    let parent = &data::PARENTS[0];
    let signed = (slot + 1) as i64;
    (
        format!(
            "v{}",
            parent
                .vertices
                .iter()
                .position(|vertex| vertex.contains(&signed))
                .unwrap()
        ),
        format!(
            "v{}",
            parent
                .vertices
                .iter()
                .position(|vertex| vertex.contains(&-signed))
                .unwrap()
        ),
    )
}
pub fn full_parent_dot() -> String {
    let edges = data::PARENTS[0]
        .physical_slots
        .iter()
        .map(|&slot| {
            let (source, target) = parent_endpoints(slot);
            format!("{source} -> {target} [particle=\"phi\"];\n")
        })
        .collect::<String>();
    format!("digraph cancellation_parent {{ {edges} }}")
}
fn ordered_physical_slots(diagram: &FeynmanDiagram) -> Result<(Vec<EdgeId>, Vec<usize>, Vec<i64>)> {
    let mut edges = diagram.edges().collect::<Vec<_>>();
    edges.sort_by_key(|(id, _, _)| *id);
    let mut slots = Vec::new();
    let mut signs = Vec::new();
    let mut ids = Vec::new();
    for (id, endpoints, _) in edges {
        let source = &diagram
            .vertex(endpoints.source)
            .ok_or("missing source vertex")?
            .name;
        let target = &diagram
            .vertex(endpoints.target)
            .ok_or("missing target vertex")?
            .name;
        let (slot, sign) = data::PARENTS[0]
            .physical_slots
            .iter()
            .find_map(|&slot| {
                let (left, right) = parent_endpoints(slot);
                if source == &left && target == &right {
                    Some((slot, 1))
                } else if source == &right && target == &left {
                    Some((slot, -1))
                } else {
                    None
                }
            })
            .ok_or("native edge lacks an authoritative parent-slot witness")?;
        ids.push(id);
        slots.push(slot);
        signs.push(sign);
    }
    if slots.iter().collect::<BTreeSet<_>>().len() != 12 {
        return Err("parent physical slots are not bijective".into());
    }
    Ok((ids, slots, signs))
}
fn substitute_mass(family: &NativeFamily, mass: &Atom, value: i64) -> Result<NativeFamily> {
    let replacement = [Replacement::new(mass.clone(), Atom::num(value))];
    Ok(NativeFamily::new(
        family.loop_momenta().to_vec(),
        vec![],
        family
            .denominators()
            .iter()
            .map(|d| d.replace_multiple(&replacement))
            .collect(),
        family.kinematics(),
    )?)
}
pub(super) fn common_family(
    loops: &[Atom],
    kinematics: &Kinematics,
    mass_squared: &Atom,
) -> Result<NativeFamily> {
    let denominators = data::MOMENTA
        .iter()
        .map(|row| {
            let momentum = row
                .iter()
                .zip(loops)
                .map(|(&coefficient, momentum)| Atom::num(coefficient) * momentum)
                .sum::<Atom>();
            kinematics
                .scalar_product(&momentum, &momentum)
                .map(|square| square - mass_squared)
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    Ok(NativeFamily::new(
        loops.to_vec(),
        vec![],
        denominators,
        kinematics,
    )?)
}
pub(super) fn unit_core_family() -> Result<IntegralFamily> {
    let c = CoefficientContext::try_new(["d"])?;
    let denominators = data::MOMENTA
        .iter()
        .map(|row| {
            AffineDenominator::new(
                c.integer(-1),
                (0..5)
                    .flat_map(|a| (a..5).map(move |b| row[a] * row[b] * if a == b { 1 } else { 2 }))
                    .map(|value| c.integer(value))
                    .collect(),
            )
        })
        .collect();
    Ok(IntegralFamily::new(
        "cancellation_common_basis",
        (0..5).map(|i| format!("q{i}")).collect(),
        vec![],
        c.clone(),
        c.parameter("d").unwrap(),
        denominators,
        vec![],
        vec![c.zero(); 15],
    )?)
}

/// Independently construct the final five-loop pinched graph. Native stable DOT
/// supports the resulting degree-ten effective vertex with no invented model
/// interaction. The graph has five nonzero massive self-loop denominators.
fn pinched_diagram(model: &Model) -> Result<(FeynmanDiagram, String)> {
    let particle = model.particle_id("phi")?;
    let mut builder = FeynmanDiagram::builder(model.clone(), "independently_pinched_five_tadpoles");
    let vertex = builder.add_vertex(DiagramVertex {
        name: "pinched".into(),
        interaction: None,
        numerator: Atom::one(),
    });
    for _ in 0..5 {
        builder.add_edge(vertex, vertex, DiagramEdge::new(particle, false))?;
    }
    let dot = builder.build()?.to_dot()?;
    let diagram = FeynmanDiagram::from_dot(model.clone(), &dot)?;
    if diagram.loop_count() != 5 || diagram.edges().count() != 5 {
        return Err("separate pinched DOT changed loop or line count".into());
    }
    Ok((diagram, dot))
}

/// Build six exact controls. In every pair, the left has one or more genuinely
/// positive physical lines canceled by its expanded numerator. The right is a
/// separately parsed pinched DOT; neither side is produced by reducing the other.
pub fn fixtures(_initialized_vakint: &Vakint) -> Result<Vec<CancellationFixture>> {
    if data::PARENTS[0].sector != 32745 {
        return Err("unexpected authoritative parent ordering".into());
    }
    let model = Model::phi_3_4();
    let mass = model
        .particle_by_id(model.particle_id("phi")?)?
        .symbolic_mass(&model);
    let diagram = FeynmanDiagram::from_dot(model.clone(), &full_parent_dot())?;
    let physical = diagram.propagator_family(&Kinematics::new())?;
    let (edge_ids, physical_slots, edge_signs) = ordered_physical_slots(&diagram)?;
    let completed = physical.complete(&[])?;
    if completed.denominators().len() != 15 || diagram.loop_count() != 5 {
        return Err("parent completion has unexpected exact dimensions".into());
    }
    let (pinched, pinched_dot) = pinched_diagram(&model)?;
    let pinched_physical = pinched.propagator_family(&Kinematics::new())?;
    let loop_labels = (1..=5)
        .map(|axis| vakint_parse!(format!("k({axis})")))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let core = unit_core_family()?;
    let lowerer = FamilyScalarNumeratorService::try_new(
        &core,
        vakint_symbol!("dot"),
        loop_labels.clone(),
        Default::default(),
    )?;
    type CancellationCase<'a> = (&'a str, &'a [(usize, i64)], i64);
    let cases: [CancellationCase<'_>; 6] = [
        ("single_positive_line", &[(0, 1)], 1),
        ("dotted_positive_line", &[(3, 2)], 1),
        ("two_positive_lines", &[(0, 1), (5, 1)], 1),
        ("mixed_dotted_lines", &[(1, 2), (6, 3)], 1),
        ("single_line_nonunit_mass", &[(0, 2)], 2),
        ("multiple_lines_nonunit_mass", &[(1, 2), (6, 3)], 2),
    ];
    let mut result = Vec::new();
    for (label, cancellations, mass_value) in cases {
        let mass_squared = Atom::num(mass_value * mass_value);
        let native = substitute_mass(&completed, &mass, mass_value)?;
        let common = common_family(native.loop_momenta(), native.kinematics(), &mass_squared)?;
        // Lift the exact stored chord routing into the authoritative common
        // basis. Unit chord signatures are checked rather than assumed.
        let mut images = vec![None; 5];
        for chord in &diagram.loop_momentum_basis().loop_edges {
            let position = edge_ids
                .iter()
                .position(|id| id == chord)
                .ok_or("unknown native chord")?;
            let coefficients = diagram.loop_momentum_basis().edge_signatures[chord]
                .loops
                .integer_coefficients();
            let nonzero = coefficients
                .iter()
                .enumerate()
                .filter(|(_, c)| **c != 0)
                .collect::<Vec<_>>();
            let [(axis, coefficient)] = nonzero.as_slice() else {
                return Err("native chord is not one independent loop coordinate".into());
            };
            if **coefficient != 1 && **coefficient != -1 {
                return Err("native chord is not unit normalized".into());
            }
            let momentum = data::MOMENTA[physical_slots[position]]
                .iter()
                .zip(common.loop_momenta())
                .map(|(&value, q)| Atom::num(value * edge_signs[position]) * q)
                .sum::<Atom>();
            if images[*axis]
                .replace(Atom::num(**coefficient) * momentum)
                .is_some()
            {
                return Err("native chord axes are repeated".into());
            }
        }
        let images = images
            .into_iter()
            .collect::<Option<Vec<_>>>()
            .ok_or("native chords omit a loop axis")?;
        let original = substitute_mass(&physical, &mass, mass_value)?;
        let mapping = original
            .mapping_to(&common, &images)?
            .ok_or("stored parent routing failed native all-row/unit-Jacobian verification")?;
        if mapping.denominator_map() != physical_slots {
            return Err("native parent map changed the explicit physical-slot ordering".into());
        }

        let mut reordered = vec![Atom::Zero; 15];
        for (source, &slot) in physical_slots.iter().enumerate() {
            reordered[slot] = native.denominators()[source].clone();
        }
        let holes = (0..15)
            .filter(|slot| !physical_slots.contains(slot))
            .collect::<Vec<_>>();
        for (slot, denominator) in holes.iter().zip(&native.denominators()[12..]) {
            reordered[*slot] = denominator.clone();
        }
        let family = NativeFamily::new(
            native.loop_momenta().to_vec(),
            vec![],
            reordered,
            native.kinematics(),
        )?;
        let mut powers = BASE;
        let mut numerator = Atom::one();
        for &(slot, power) in cancellations {
            if BASE[slot] != 0 || !physical_slots.contains(&slot) || power <= 0 {
                return Err("cancellation does not add an absent positive physical line".into());
            }
            powers[slot] = power;
            numerator *= family.denominators()[slot].clone().pow(Atom::num(power));
        }
        let numerator = numerator.expand();
        let cancelled = Atom::from(VakintExpression::from_diagram(
            &diagram,
            &family,
            &numerator,
            &DiagramIntegralOptions {
                powers: Some(powers.into_iter().map(Integer::from).collect()),
                physical_slots: Some(physical_slots.clone()),
                parameter_substitutions: vec![(mass.clone(), Atom::num(mass_value))],
                ..Default::default()
            },
        )?);

        // Structural cancellation in the EXACT same common family as the saved
        // base key, before any entry admission or reduction call.
        let products = (0..5)
            .flat_map(|left| (left..5).map(move |right| (left, right)))
            .map(|(left, right)| {
                Ok(Replacement::new(
                    common.kinematics().scalar_product(
                        &common.loop_momenta()[left],
                        &common.loop_momenta()[right],
                    )?,
                    vakint_symbol!("dot").call((&loop_labels[left], &loop_labels[right])),
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        let common_numerator = mapping
            .apply(&numerator)
            .replace_multiple(&products)
            .expand();
        // Independently replay ingress's simultaneous numerator/denominator
        // transport, including its non-prefix physical_slots arrangement.
        let native_to_vakint = common
            .loop_momenta()
            .iter()
            .zip(&loop_labels)
            .map(|(source, target)| Replacement::new(source.clone(), target.clone()))
            .collect::<Vec<_>>();
        let routing = loop_labels
            .iter()
            .zip(&images)
            .map(|(source, target)| {
                Replacement::new(source.clone(), target.replace_multiple(&native_to_vakint))
            })
            .collect::<Vec<_>>();
        let split = VakintExpression::split_integrals(cancelled.as_view())?;
        if split.len() != 1
            || !(split[0].numerator.replace_multiple(&routing) - &common_numerator)
                .expand()
                .together()
                .is_zero()
        {
            return Err(
                format!("{label}: ingress changed the simultaneously routed numerator").into(),
            );
        }
        let pattern = vakint_parse!("prop(id_,edge(a_,b_),q_,mass_,power_)")?.to_pattern();
        for matched in split[0].integral.pattern_match(&pattern, None, None) {
            let id = usize::try_from(matched[&vakint_symbol!("id_")].as_view())?;
            let position = edge_ids
                .iter()
                .position(|edge| edge.0 == id)
                .ok_or("ingress changed a physical edge identity")?;
            let slot = physical_slots[position];
            let q = matched[&vakint_symbol!("q_")].replace_multiple(&routing);
            let denominator =
                vakint_symbol!("dot").call((&q, &q)) - &matched[&vakint_symbol!("mass_")];
            let expected = common.denominators()[slot].replace_multiple(&products);
            if !(denominator - expected).expand().together().is_zero()
                || matched[&vakint_symbol!("power_")] != Atom::num(powers[slot])
            {
                return Err(
                    format!("{label}: ingress changed a physical denominator/power").into(),
                );
            }
        }
        let lowered = lowerer.lower_weighted(
            &common_numerator,
            &IntegralKey::try_new(powers)?,
            &mass_squared,
        )?;
        if lowered.terms().len() != 1
            || lowered.terms()[0].integral().powers() != BASE
            || lowered.terms()[0].coefficient() != &core.coefficient_context().one()
            || lowered.terms()[0].scalar_spectator() != &Atom::one()
        {
            return Err(format!(
                "{label}: native weighted cancellation did not yield the sole admitted base key"
            )
            .into());
        }

        let active = (0..15).filter(|&slot| BASE[slot] > 0).collect::<Vec<_>>();
        let base_family = NativeFamily::new(
            common.loop_momenta().to_vec(),
            vec![],
            active
                .iter()
                .map(|&slot| common.denominators()[slot].clone())
                .collect(),
            common.kinematics(),
        )?;
        let pinched_family = substitute_mass(&pinched_physical, &mass, mass_value)?;
        let pinched_map = base_family
            .find_mapping(&pinched_family, 1_000_000)?
            .ok_or("separately pinched DOT lacks an exact unit-measure base map")?;
        let mut pinched_powers = vec![Integer::from(0); 5];
        for (source, &target) in pinched_map.denominator_map().iter().enumerate() {
            pinched_powers[target] = Integer::from(BASE[active[source]]);
        }
        let pinched_input = Atom::from(VakintExpression::from_diagram(
            &pinched,
            &pinched_family,
            &Atom::one(),
            &DiagramIntegralOptions {
                powers: Some(pinched_powers),
                parameter_substitutions: vec![(mass.clone(), Atom::num(mass_value))],
                ..Default::default()
            },
        )?);
        result.push(CancellationFixture {
            label: label.into(),
            cancelled,
            pinched: pinched_input,
            expected_key: IntegralKey::try_new(BASE)?,
            mass_squared,
            positive_cancellations: cancellations.to_vec(),
            pinched_dot: pinched_dot.clone(),
        });
    }
    Ok(result)
}

/// Additional no-package ingress control: the numerator and every physical row
/// are replayed under native DOT roundtrip routing, never an inferred rematch.
pub fn verify_fixture_shapes(fixtures: &[CancellationFixture]) -> Result<()> {
    for fixture in fixtures {
        let cancelled = VakintExpression::split_integrals(fixture.cancelled.as_view())?;
        let pinched = VakintExpression::split_integrals(fixture.pinched.as_view())?;
        if cancelled.len() != 1 || pinched.len() != 1 {
            return Err("cancellation fixture unexpectedly split its topology".into());
        }
        let pattern = vakint_parse!("prop(id_,edge(a_,b_),q_,mass_,power_)")?.to_pattern();
        let count = |expression: &Atom| expression.pattern_match(&pattern, None, None).count();
        if count(&cancelled[0].integral) != 5 + fixture.positive_cancellations.len()
            || count(&pinched[0].integral) != 5
        {
            return Err("fixture did not retain genuinely positive canceled lines".into());
        }
        if cancelled[0].numerator.is_one() || !pinched[0].numerator.is_one() {
            return Err("cancellation and pinched numerators have incorrect roles".into());
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn expanded_positive_propagators_cancel_to_the_admitted_base_before_reduction() {
        let vakint = Vakint::new().unwrap();
        let fixtures = fixtures(&vakint).unwrap();
        verify_fixture_shapes(&fixtures).unwrap();
        assert_eq!(fixtures.len(), 6);
        assert!(fixtures.iter().any(|f| f.positive_cancellations.len() > 1));
        assert!(
            fixtures
                .iter()
                .any(|f| f.positive_cancellations.iter().any(|(_, power)| *power > 1))
        );
        assert!(fixtures.iter().any(|f| f.mass_squared == Atom::num(4)));
        for fixture in fixtures {
            assert_eq!(fixture.expected_key.powers(), BASE);
        }
    }
}
