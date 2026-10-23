//! Exact checks of the integral represented by each external-data master DOT.

use std::{collections::BTreeMap, error::Error, fs, path::Path, sync::Arc, time::Instant};

use feynkit_graph::{FeynmanDiagram, IntegralFamily};
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use linnet::parser::DotGraph;
use serde_json::{Value, json};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::integer::Integer,
    id::Replacement,
};
use vakint::{VakintExpression, diagram_integral::DiagramIntegralOptions, symbols::S};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

fn string<'a>(value: &'a Value, field: &str) -> Result<&'a str> {
    value[field]
        .as_str()
        .ok_or_else(|| format!("missing string field {field}").into())
}

fn check(condition: bool, message: &str) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(message.into())
    }
}

fn equal(left: &Atom, right: &Atom, message: &str) -> Result<()> {
    check((left - right).expand().is_zero(), message)
}

fn parse(expression: &str) -> Result<Atom> {
    Ok(Atom::parse(
        expression,
        "external_data_acceptance",
        Default::default(),
    )?)
}

/// Prepare a previously validated master graph through native Vakint ingress.
pub(super) fn input(root: &Path, master: &Value) -> Result<Atom> {
    input_powers(root, master, None)
}

/// Keep a verified graph and numerator while changing its positive line powers.
pub(super) fn input_powers(
    root: &Path,
    master: &Value,
    family_powers: Option<&[i64]>,
) -> Result<Atom> {
    let model = Model::from_json(&fs::read_to_string(root.join("model.json"))?)?;
    let dot = fs::read_to_string(root.join(string(master, "directory")?).join("graph.dot"))?;
    let raw: DotGraph = DotGraph::from_string(&dot).map_err(|error| error.to_string())?;
    let diagram = FeynmanDiagram::from_dot(model, &dot)?;
    let kin = Kinematics::in_dimension(&parse("gammalooprs::dim")?)?;
    let family = diagram.propagator_family(&kin)?;
    let powers = raw
        .iter_edges()
        .map(|(_, id, data)| {
            let attrs = &data.data.statements;
            let power = match family_powers {
                Some(powers) => *powers
                    .get(
                        attrs
                            .get("common_basis_slot")
                            .ok_or("missing common-basis slot")?
                            .parse::<usize>()?,
                    )
                    .ok_or("common-basis slot outside power vector")?,
                None => attrs
                    .get("power")
                    .ok_or("missing edge power")?
                    .parse::<i64>()?,
            };
            check(power > 0, "master graph has a nonpositive physical power")?;
            Ok((id.0, Integer::from(power)))
        })
        .collect::<Result<BTreeMap<_, _>>>()?;
    let options = DiagramIntegralOptions {
        powers: Some(powers.into_values().collect()),
        ..Default::default()
    };
    Ok(Atom::from(VakintExpression::from_diagram(
        &diagram,
        &family,
        diagram.numerator(),
        &options,
    )?))
}

/// An uncut two-vertex parallel-edge vacuum graph has one cycle-matroid block.
/// Requiring unit powers and no numerator identifies the scalar banana control.
pub(super) fn is_scalar_banana(root: &Path, master: &Value) -> Result<bool> {
    let powers: Vec<i64> = serde_json::from_value(master["powers"].clone())?;
    if powers.iter().any(|power| !matches!(power, 0 | 1)) {
        return Ok(false);
    }
    let model = Model::from_json(&fs::read_to_string(root.join("model.json"))?)?;
    let source = fs::read_to_string(root.join(string(master, "directory")?).join("graph.dot"))?;
    let diagram = FeynmanDiagram::from_dot(model, &source)?;
    Ok(diagram.vertices().count() == 2
        && diagram
            .edges()
            .all(|(_, endpoints, _)| endpoints.source != endpoints.target)
        && diagram.numerator().is_one())
}

/// `root` is the loop-pack directory; all paths in its manifest are relative.
pub(super) fn validate(root: &Path, manifest: &Value) -> Result<Value> {
    let masters = manifest["masters"].as_array().ok_or("missing masters")?;
    let families = manifest["families"].as_array().ok_or("missing families")?;
    let model = Arc::new(Model::from_json(&fs::read_to_string(
        root.join("model.json"),
    )?)?);
    let mut results = Vec::with_capacity(masters.len());
    for master in masters {
        let name = string(master, "name")?;
        let family = families
            .iter()
            .find(|family| family["sha256"] == master["family_sha256"])
            .ok_or("master refers to an unknown family")?;
        let result = validate_one(root, master, family, Arc::clone(&model))
            .map_err(|error| format!("{name}: {error}"))?;
        results.push(result);
    }
    Ok(json!({
        "masters_checked": results.len(),
        "dotted_masters": results.iter().filter(|r| r["dotted"] == true).count(),
        "numerator_masters": results.iter().filter(|r| r["numerator_slots"].as_u64().unwrap_or(0) > 0).count(),
        "exact_native_mapping": true,
        "exact_complete_numerator": true,
        "native_dot_roundtrip": true,
        "native_vakint_ingress": true,
        "results": results,
    }))
}

fn validate_one(root: &Path, master: &Value, family: &Value, model: Arc<Model>) -> Result<Value> {
    let started = Instant::now();
    let name = string(master, "name")?;
    eprintln!("DOT {name}: native parse and validation");
    let dot = fs::read_to_string(root.join(string(master, "directory")?).join("graph.dot"))?;
    let raw: DotGraph = DotGraph::from_string(&dot).map_err(|error| error.to_string())?;
    let attrs = &raw.global_data.statements;
    let powers: Vec<i32> = serde_json::from_value(master["powers"].clone())?;
    let stored_powers: Vec<i32> =
        serde_json::from_str(attrs.get("family_powers").ok_or("missing family_powers")?)?;
    check(
        powers == stored_powers,
        "DOT signed powers differ from manifest",
    )?;
    let transform: Vec<Vec<i64>> = serde_json::from_str(
        attrs
            .get("common_basis_to_graph_loops")
            .ok_or("missing common-basis routing")?,
    )?;
    let diagram = FeynmanDiagram::from_dot(model.clone(), &dot)?;
    diagram.validate()?;
    check(diagram.cuts().is_empty(), "master graph contains cuts")?;
    for factor in [
        diagram.overall_factor(),
        diagram.numerator_prefactor(),
        diagram.projector(),
    ] {
        check(factor.is_one(), "unexpected graph prefactor/projector")?;
    }
    // The generated DOT carries symbolic d in its native Spenso scalar products.
    let kin = Kinematics::in_dimension(&parse("gammalooprs::dim")?)?;
    let target = diagram.propagator_family(&kin)?;
    let loop_count = target.loop_momenta().len();
    check(
        master["loop_count"].as_u64() == Some(loop_count as u64)
            && diagram.loop_count() == loop_count,
        "graph changed the integration loop count",
    )?;
    check(
        transform.len() == loop_count && transform.iter().all(|row| row.len() == loop_count),
        "routing matrix has the wrong dimension",
    )?;
    check(
        family["coordinate_index_base"] == 0
            && family["external_momenta"]
                .as_array()
                .is_some_and(Vec::is_empty),
        "expected a zero-based vacuum family",
    )?;
    let shifts = family["power_shifts"]
        .as_array()
        .ok_or("missing power shifts")?;
    check(shifts.len() == powers.len(), "power-shift arity mismatch")?;
    for shift in shifts {
        check(
            parse(shift.as_str().ok_or("invalid power shift")?)?.is_zero(),
            "noninteger power shifts are not represented by this master graph",
        )?;
    }
    let mut edges = diagram.edges().collect::<Vec<_>>();
    edges.sort_by_key(|(id, _, _)| *id);
    let first_edge = edges.first().ok_or("master graph has no propagators")?;
    let mass = model
        .particle_by_id(first_edge.2.particle)?
        .symbolic_mass_parameter(&model)?;
    for (_, _, edge) in &edges {
        equal(
            &model
                .particle_by_id(edge.particle)?
                .symbolic_mass_parameter(&model)?,
            &mass,
            "master graph does not have a common mass",
        )?;
    }
    let mass_squared = mass.pow(2);
    let source_loops = family["loop_momenta"]
        .as_array()
        .ok_or("missing source loop momenta")?
        .iter()
        .map(|name| parse(name.as_str().ok_or("invalid loop name")?))
        .collect::<Result<Vec<_>>>()?;
    check(
        source_loops.len() == loop_count,
        "source loop count mismatch",
    )?;
    let source_kin = kin.clone().with_momenta(source_loops.iter().cloned())?;
    let coordinates = family["coordinates"]
        .as_array()
        .ok_or("missing coordinates")?;
    let products = coordinates
        .iter()
        .map(|coordinate| {
            check(coordinate["kind"] == "loop_loop", "non-vacuum coordinate")?;
            let left = coordinate["left"]
                .as_u64()
                .ok_or("missing left coordinate")? as usize;
            let right = coordinate["right"]
                .as_u64()
                .ok_or("missing right coordinate")? as usize;
            Ok(source_kin.scalar_product(
                source_loops
                    .get(left)
                    .ok_or("left coordinate out of bounds")?,
                source_loops
                    .get(right)
                    .ok_or("right coordinate out of bounds")?,
            )?)
        })
        .collect::<Result<Vec<_>>>()?;
    let denominators = family["denominators"]
        .as_array()
        .ok_or("missing denominators")?
        .iter()
        .map(|denominator| {
            let coefficients = denominator["coefficients"]
                .as_array()
                .ok_or("missing coefficients")?;
            check(
                coefficients.len() == products.len(),
                "denominator coordinate count mismatch",
            )?;
            let mut value = parse(string(denominator, "constant")?)? * &mass_squared;
            for (coefficient, product) in coefficients.iter().zip(&products) {
                value += parse(coefficient.as_str().ok_or("invalid coefficient")?)? * product;
            }
            Ok(value)
        })
        .collect::<Result<Vec<_>>>()?;
    check(
        denominators.len() == powers.len(),
        "denominator/power arity mismatch",
    )?;
    let mut positive = Vec::new();
    let mut positive_powers = Vec::new();
    let mut positive_slots = Vec::new();
    for (slot, (denominator, power)) in denominators.iter().zip(&powers).enumerate() {
        if *power > 0 {
            positive.push(denominator.clone());
            positive_powers.push(*power);
            positive_slots.push(slot);
        }
    }
    eprintln!(
        "DOT {name}: exact denominator and Jacobian witness ({:.3}s)",
        started.elapsed().as_secs_f64()
    );
    let source = IntegralFamily::new(source_loops, vec![], positive, &source_kin)?;
    let images = transform
        .iter()
        .map(|row| {
            row.iter()
                .zip(target.loop_momenta())
                .map(|(coefficient, momentum)| Atom::num(*coefficient) * momentum)
                .sum::<Atom>()
        })
        .collect::<Vec<_>>();
    // This native witness checks the Jacobian and every actual massive denominator.
    let witness = source
        .mapping_to(&target, &images)?
        .ok_or("invalid denominator/Jacobian witness")?;
    // The native momentum substitution is an algebra homomorphism. Applying it
    // to every complete denominator factor and retaining each signed power proves
    // the complete numerator identity without expanding products of quadratics.
    let transported = denominators
        .iter()
        .zip(&powers)
        .filter(|(_, power)| **power < 0)
        .map(|(denominator, power)| witness.apply(denominator).pow(-i64::from(*power)))
        .product::<Atom>();
    eprintln!(
        "DOT {name}: complete factored numerator ({:.3}s)",
        started.elapsed().as_secs_f64()
    );
    equal(
        &transported,
        diagram.numerator(),
        "complete numerator differs after exact routing",
    )?;
    let edge_attrs = raw
        .iter_edges()
        .map(|(_, id, data)| (id.0, data.data.statements.clone()))
        .collect::<BTreeMap<_, _>>();
    check(
        edge_attrs.len() == edges.len(),
        "parsed/native edge count mismatch",
    )?;
    let edge_powers = edges
        .iter()
        .map(|(id, _, _)| {
            Ok(edge_attrs[&id.0]
                .get("power")
                .ok_or("missing edge power")?
                .parse::<i32>()?)
        })
        .collect::<Result<Vec<_>>>()?;
    check(
        witness.map_powers(&positive_powers)? == edge_powers,
        "edge powers differ after routing",
    )?;
    check(
        edges
            .iter()
            .map(|(_, _, edge)| edge.propagator_power)
            .eq(edge_powers.iter().copied()),
        "native DOT parser did not preserve physical propagator powers",
    )?;
    for (source_index, target_index) in witness.denominator_map().iter().enumerate() {
        let id = edges[*target_index].0;
        let slot: usize = edge_attrs[&id.0]
            .get("common_basis_slot")
            .ok_or("missing edge source slot")?
            .parse()?;
        check(
            slot == positive_slots[source_index],
            "edge common-basis slot differs from exact witness",
        )?;
    }
    let basis = diagram.loop_momentum_basis();
    let mut lmb_slots = BTreeMap::new();
    for (id, attrs) in &edge_attrs {
        if let Some(slot) = attrs.get("lmb_id") {
            check(
                lmb_slots.insert(slot.parse::<usize>()?, *id).is_none(),
                "duplicate lmb_id",
            )?;
        }
    }
    check(
        lmb_slots.len() == loop_count
            && basis
                .loop_edges
                .iter()
                .enumerate()
                .all(|(slot, id)| lmb_slots.get(&slot) == Some(&id.0)),
        "lmb_id annotations differ from the native stored loop basis",
    )?;
    let declared_basis: Value = serde_json::from_str(
        attrs
            .get("loop_momentum_basis")
            .ok_or("missing stored routing")?,
    )?;
    for (id, signature) in &basis.edge_signatures {
        let declared: Vec<isize> = serde_json::from_value(
            declared_basis["edge_signatures"][id.0.to_string()]["loops"].clone(),
        )?;
        check(
            signature.loops.integer_coefficients() == declared,
            "serialized/native momentum routing differs",
        )?;
    }
    let options = DiagramIntegralOptions {
        powers: Some(edge_powers.iter().copied().map(Integer::from).collect()),
        ..Default::default()
    };
    eprintln!(
        "DOT {name}: native Vakint ingress ({:.3}s)",
        started.elapsed().as_secs_f64()
    );
    let expression =
        VakintExpression::from_diagram(&diagram, &target, diagram.numerator(), &options)?;
    let mapped = (0..loop_count).map(|i| S.k.call(i + 1)).collect::<Vec<_>>();
    let mut replacements = Vec::new();
    for (i, left) in target.loop_momenta().iter().enumerate() {
        for (j, right) in target.loop_momenta().iter().enumerate().skip(i) {
            replacements.push(Replacement::new(
                target.kinematics().scalar_product(left, right)?,
                S.dot.call((&mapped[i], &mapped[j])),
            ));
        }
    }
    let expected_numerator = transported.replace_multiple(&replacements);
    let mut expected_topology = Atom::one();
    for ((id, endpoints, _), power) in edges.iter().zip(&edge_powers) {
        let momentum = basis.edge_signatures[id]
            .loops
            .apply(&mapped)?
            .unwrap_or_default();
        expected_topology *= S.prop.call((
            id.0,
            S.edge.call((endpoints.source.0, endpoints.target.0)),
            momentum,
            &mass_squared,
            *power,
        ));
    }
    let expected = expected_numerator * S.topo.call(expected_topology);
    equal(
        &Atom::from(expression.clone()),
        &expected,
        "Vakint ingress changed numerator, mass, powers or routing",
    )?;
    // Check both explicit caller powers and the native graph-power defaults.
    eprintln!(
        "DOT {name}: native roundtrip ({:.3}s)",
        started.elapsed().as_secs_f64()
    );
    let roundtrip = FeynmanDiagram::from_dot(model, &diagram.to_dot()?)?;
    let roundtrip_family = roundtrip.propagator_family(&kin)?;
    let roundtrip_expression = VakintExpression::from_diagram(
        &roundtrip,
        &roundtrip_family,
        roundtrip.numerator(),
        &options,
    )?;
    equal(
        &Atom::from(roundtrip_expression),
        &expected,
        "native DOT round-trip changed the integral",
    )?;
    let default_expression = VakintExpression::from_diagram(
        &roundtrip,
        &roundtrip_family,
        roundtrip.numerator(),
        &DiagramIntegralOptions::default(),
    )?;
    equal(
        &Atom::from(default_expression),
        &expected,
        "native DOT power defaults changed the integral",
    )?;
    Ok(json!({
        "name": string(master, "name")?,
        "loop_count": loop_count,
        "propagators": edges.len(),
        "dotted": edge_powers.iter().any(|power| *power > 1),
        "numerator_slots": powers.iter().filter(|power| **power < 0).count(),
        "numerator_rank": powers.iter().filter(|power| **power < 0).map(|power| -i64::from(*power)).sum::<i64>(),
        "lmb_zero_present": lmb_slots.contains_key(&0),
        "exact_mapping_and_ingress": true,
        "native_graph_power_defaults": true,
        "seconds": started.elapsed().as_secs_f64(),
    }))
}
