//! Homogeneous quadratic Taylor coefficient, including every propagator term.

use super::*;
use super::{projector, vacuum::Vacuum};
use idenso::CookSettings;
use spenso::structure::{
    abstract_index::AbstractIndex,
    representation::{LibrarySlot, Minkowski},
    slot::ParseableAind,
};
use symbolica::atom::Symbol;

pub(super) struct Contribution {
    pub powers: Vec<i32>,
    pub numerator: Atom,
}

pub(super) struct Projected {
    pub vacuum: Vacuum,
    pub contributions: Vec<Contribution>,
    pub conditions: BTreeSet<Atom>,
}

fn powers(monomial: &Atom, inverses: &[Atom]) -> Result<Vec<i32>> {
    let factors = match monomial.as_view() {
        AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
        term => vec![term],
    };
    let mut powers = vec![0i32; inverses.len()];
    for factor in factors {
        if factor.is_one() {
            continue;
        }
        let (base, exponent) = if let AtomView::Pow(power) = factor {
            let (base, exponent) = power.get_base_exp();
            (base, i32::try_from(i64::try_from(exponent)?)?)
        } else {
            (factor, 1)
        };
        let slot = inverses
            .iter()
            .position(|inverse| inverse.as_view() == base)
            .ok_or("unexpected non-propagator in native coefficient-list monomial")?;
        if exponent < 0 {
            return Err("Taylor collection returned a negative inverse-propagator power".into());
        }
        powers[slot] = powers[slot]
            .checked_add(exponent)
            .ok_or("propagator power overflow")?;
    }
    Ok(powers)
}

fn coefficient(expression: &Atom, parameter: Symbol, power: i64) -> Result<Atom> {
    Ok(expression
        .series(parameter, 0, 2)?
        .coefficient(power.into())
        .ok_or("native Taylor coefficient is outside the computed remainder")?)
}

type GluonPorts = Vec<Vec<LibrarySlot<AbstractIndex>>>;

fn amputated_source(diagram: &FeynmanDiagram, d: &Atom) -> Result<(Atom, GluonPorts)> {
    let raw = lorentz_dimension(&(diagram.numerator() * diagram.numerator_prefactor()), d);
    // Generated FeynKit indices are function payloads; the local Idenso index
    // type consumes symbols. Keep the existing reversible encoding boundary.
    let cooking = CookSettings::reversible();
    let numerator = cooking
        .try_cook_indices(raw.as_view())
        .map_err(|error| format!("failed to encode generated tensor indices: {error:?}"))?;
    let ports = numerator
        .list_dangling::<AbstractIndex>()?
        .iter()
        .map(|port| LibrarySlot::<AbstractIndex>::try_from(port.as_view()))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let mut legs = BTreeMap::new();
    for (id, ends, _) in diagram.edges() {
        let source = diagram.vertex(ends.source).expect("valid source");
        let target = diagram.vertex(ends.target).expect("valid target");
        let (external, index_head) = if let Some(external) = &source.external {
            (external, symbol!("FeynKit::SinkIndex"))
        } else if let Some(external) = &target.external {
            (external, symbol!("FeynKit::SourceIndex"))
        } else {
            continue;
        };
        let tensor_index = index_head.call((id.0, 1));
        let mut slots = ports
            .iter()
            .filter(|slot| cooking.uncook(slot.aind.to_atom().as_view()) == tensor_index)
            .copied()
            .collect::<Vec<_>>();
        slots.sort();
        if slots.len() != 2 || legs.insert(external.index, slots).is_some() {
            return Err(
                "native amputated source must retain two gluon Lorentz/color ports per leg".into(),
            );
        }
    }
    let legs = legs.into_values().collect::<Vec<_>>();
    if legs.len() != 2 || ports.len() != 4 {
        return Err("native amputated source has unassigned or missing external ports".into());
    }
    Ok((numerator, legs))
}

pub(super) fn prepare(
    generated: feynkit_generator::GenerationResult,
    d: &Atom,
    evidence: &Path,
) -> Result<Projected> {
    projector::controls(d)?;
    if !generated.report.completed || generated.diagrams.len() != 1 {
        return Err("native amputated source requires one completely generated diagram".into());
    }
    let diagram = &generated.diagrams[0];
    let (numerator, legs) = amputated_source(diagram, d)?;
    let mut fermions = 0;
    let mut gluons = 0;
    for (_, ends, edge) in diagram.edges() {
        if is_external_edge(diagram, ends) {
            if !edge.numerator.is_one() {
                return Err("source contains external propagators or wavefunctions".into());
            }
            continue;
        }
        let particle = diagram.model().particle_by_id(edge.particle)?;
        if !particle.symbolic_mass_parameter(diagram.model())?.is_zero() {
            return Err("original QCD numerator is not massless".into());
        }
        match particle.pdg_code.abs() {
            1 => fermions += 1,
            21 => gluons += 1,
            _ => return Err("unexpected internal species".into()),
        }
    }
    if (fermions, gluons) != (10, 4) {
        return Err(
            "native source does not have ten massless fermions and four gluon rungs".into(),
        );
    }
    let vacuum = Vacuum::new(diagram, d)?;
    let q = spenso::vector_symbol!("qcd_ladder::q");
    let t = symbol!("qcd_ladder::taylor_scale");
    let t_atom = Atom::var(t);
    let independent = diagram
        .loop_momentum_basis()
        .external_edges
        .iter()
        .enumerate()
        .filter(|(_, edge)| {
            !diagram
                .loop_momentum_basis()
                .dependent_externals
                .contains(edge)
        })
        .collect::<Vec<_>>();
    let [(external_axis, _)] = independent.as_slice() else {
        return Err("two-point source must have exactly one independent external momentum".into());
    };
    let raw = numerator;
    assert_lorentz_dimension(&raw, d)?;
    let routed = diagram.loop_momentum_basis().route_expression(&raw);
    let args = Atom::var(symbol!("qcd_ladder::tensor_args___"));
    let scaled = routed
        .replace(
            symbols::external_momentum()
                .call((*external_axis, args.clone()))
                .to_pattern(),
        )
        .with((&t_atom * q.call(args)).to_pattern());
    if scaled.contains_symbol(symbols::external_momentum()) {
        return Err(
            "source retained an undeclared external momentum after simultaneous Taylor routing"
                .into(),
        );
    }

    let mut color = Atom::one();
    let mut metric = Atom::one();
    let mut longitudinal = Atom::one();
    for slot in &legs[0] {
        let other = legs[1]
            .iter()
            .find(|other| other.rep_name() == slot.rep_name())
            .ok_or("external representations differ")?;
        if slot.rep_name() == (Minkowski {}).into() {
            let left = lorentz_dimension(&slot.to_atom(), d);
            let right = lorentz_dimension(&other.to_atom(), d);
            metric *= ETS.metric(&left, &right);
            longitudinal *= q.call(left) * q.call(right);
        } else {
            color *= ETS.metric(slot.to_atom(), other.to_atom());
        }
    }
    let mut metric_coefficients = Vec::new();
    let mut longitudinal_coefficients = Vec::new();
    for power in 0..=2 {
        let numerator = coefficient(&scaled, t, power)?;
        eprintln!("qcd_ladder phase=native_gamma_color taylor_order={power}");
        metric_coefficients.push(projector::scalar(&numerator * &color * &metric)?);
        longitudinal_coefficients.push(projector::scalar(numerator * &color * &longitudinal)?);
    }

    let inverses = (0..12)
        .map(|slot| Atom::var(symbol!(&format!("qcd_ladder::inverse_{slot}"))))
        .collect::<Vec<_>>();
    let loops = (0..5)
        .map(|axis| symbols::loop_momentum().call(axis))
        .collect::<Vec<_>>();
    let kin =
        Kinematics::in_dimension(d)?.with_momenta(loops.iter().cloned().chain([Atom::var(q)]))?;
    let q_square = kin.scalar_product(&Atom::var(q), &Atom::var(q))?;
    let rows = loop_rows(diagram);
    let mut denominator = Atom::one();
    for (id, ends, _) in diagram.edges() {
        if is_external_edge(diagram, ends) {
            continue;
        }
        let inverse = &inverses[vacuum.original_slot(id)?];
        let external = diagram.loop_momentum_basis().edge_signatures[&id]
            .external
            .integer_coefficients()[*external_axis];
        let linear =
            Atom::num(2 * external as i64) * kin.scalar_product(&rows[&id], &Atom::var(q))?;
        let quadratic = Atom::num((external * external) as i64) * &q_square;
        // Exactly 1/(D_i+t*b_i+t²*c_i), with w_i=1/D_i. Native series
        // keeps positive powers of w_i, so later mass insertion acts ONLY
        // on D_i and never turns a slash momentum into a massive numerator.
        denominator *= inverse
            / (Atom::one()
                + &t_atom * linear * inverse
                + t_atom.clone().pow(2) * quadratic * inverse);
    }
    let denominator = denominator.series(t, 0, 2)?;
    let mut metric_part = Atom::Zero;
    let mut longitudinal_part = Atom::Zero;
    for power in 0..=2 {
        let factor = denominator
            .coefficient((2 - power as i64).into())
            .ok_or("missing native denominator Taylor coefficient")?;
        metric_part += &metric_coefficients[power] * &factor;
        longitudinal_part += &longitudinal_coefficients[power] * factor;
    }
    let mut groups = BTreeMap::<Vec<i32>, (Atom, Atom)>::new();
    for (monomial, value) in metric_part.coefficient_list::<u16>(&inverses) {
        groups
            .entry(powers(&monomial, &inverses)?)
            .or_insert_with(|| (Atom::Zero, Atom::Zero))
            .0 += value;
    }
    for (monomial, value) in longitudinal_part.coefficient_list::<u16>(&inverses) {
        groups
            .entry(powers(&monomial, &inverses)?)
            .or_insert_with(|| (Atom::Zero, Atom::Zero))
            .1 += value;
    }
    let overall = evaluate_overall_factor(diagram.overall_factor().as_view());
    let mut conditions = BTreeSet::new();
    let mut contributions = Vec::new();
    for (index, (powers, (metric, longitudinal))) in groups.into_iter().enumerate() {
        if powers.len() != 12 || powers.iter().sum::<i32>() < 14 || powers.iter().sum::<i32>() > 16
        {
            return Err(
                "quadratic Taylor coefficient has invalid raised-propagator bookkeeping".into(),
            );
        }
        eprintln!(
            "qcd_ladder phase=native_q_angular group={index} A={}",
            powers.iter().sum::<i32>()
        );
        let numerator =
            projector::quadratic_transverse(&metric, &longitudinal, q, d, 8, &mut conditions)?
                * &overall;
        if numerator.is_zero() {
            continue;
        }
        assert_lorentz_dimension(&numerator, d)?;
        let rank = powers.iter().sum::<i32>() - 10;
        let scale = Atom::var(symbol!("qcd_ladder::loop_scale"));
        let scaled = numerator.replace_map(|view, _, output| {
            if let AtomView::Fun(fun) = view
                && fun.get_symbol() == SPENSO_TAG.dot
            {
                **output = view.to_owned() * scale.clone().pow(2);
            }
        });
        if !(scaled - &numerator * scale.pow(i64::from(2 * rank)))
            .expand()
            .together()
            .is_zero()
        {
            return Err(
                "projected term fails homogeneous kinetic-coefficient power counting A-R=10".into(),
            );
        }
        let numerator = numerator.collect_factors();
        std::fs::write(
            evidence.join(format!("projected-numerator-{index:03}.symbolica")),
            numerator.to_canonical_string(),
        )?;
        contributions.push(Contribution { powers, numerator });
    }
    if contributions.is_empty() {
        return Err("the complete projected source vanished unexpectedly".into());
    }
    Ok(Projected {
        vacuum,
        contributions,
        conditions,
    })
}

#[cfg(test)]
#[test]
fn generated_qcd_external_ports_are_reversibly_encoded() -> Result<()> {
    let generated = dressed_diagram()?;
    let diagram = &generated.diagrams[0];
    let d = dimension()?;
    let (numerator, legs) = amputated_source(diagram, &d)?;
    assert_eq!(legs.len(), 2);
    assert!(legs.iter().all(|slots| slots.len() == 2));
    assert_eq!(
        CookSettings::reversible().uncook(numerator.as_view()),
        lorentz_dimension(&(diagram.numerator() * diagram.numerator_prefactor()), &d)
    );
    Ok(())
}
