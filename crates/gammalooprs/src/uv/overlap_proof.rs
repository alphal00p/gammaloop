//! Sufficient, factor-preserving certificates for the physical soft pole part.
//! Numerator tensors never enter scalar rational-polynomial conversion.

use std::collections::BTreeMap;

use color_eyre::eyre::{Result, ensure, eyre};
use itertools::Itertools;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    function,
    id::Replacement,
    symbol,
};

use super::PhysicalOverlapInput;
use crate::{
    cff::expression::OrientationID,
    debug_tags,
    numerator::{
        soft_energy_certificate::{EnergyZeroOutcome, SoftEnergyAlgebra},
        symbolica_ext::NumeratorAtomExt,
    },
    utils::{GS, W_},
    uv::{
        approx::direct_3d::DirectResidueBranches, overlap_control::ExactControlPoint, uv_graph::UVE,
    },
};

#[path = "overlap_numerators.rs"]
mod numerators;
use numerators::NumeratorJetCache;

const LAST_FORBIDDEN_POWER: i64 = -2;
const MAX_GROUPS: usize = 32768;
const MAX_COEFFICIENT_BYTES: usize = 64 * 1024 * 1024;

struct PhysicalSoftMap {
    routing: Vec<Replacement>,
    deformation: Vec<Replacement>,
    energies: Vec<Replacement>,
    lambda: Symbol,
}

impl PhysicalSoftMap {
    fn new(input: &PhysicalOverlapInput) -> Self {
        let graph = &input.graph;
        let lambda = symbol!("gammalooprs::overlap_certificate::lambda");
        let full = graph.full_filter();
        // Route external edges too: the outgoing photon is dependent, not a
        // second independent external four-momentum. UV-only wrappers skip it.
        let mut routing = Vec::new();
        for (_, edge, _) in graph.iter_edges_of(&full) {
            for head in [GS.emr_vec, GS.emr_mom] {
                let signature = graph
                    .loop_momentum_basis
                    .loop_atom(edge, head, &[W_.x___], true)
                    + graph
                        .loop_momentum_basis
                        .ext_atom(edge, head, &[W_.x___], true);
                routing.push(Replacement::new(
                    function!(head, usize::from(edge), W_.x___).to_pattern(),
                    signature,
                ));
            }
        }
        let energies = graph
            .iter_edges_of(&full)
            .filter(|(pair, _, _)| pair.is_paired())
            .map(|(_, edge, data)| {
                Replacement::new(
                    GS.ose(edge).to_pattern(),
                    GS.ose_full(edge, 0.into(), data.data.mass_atom(), None),
                )
            })
            .collect();
        // Canonical basis [6,5,8]: Q6=λq, Q5=r+(1−λ)q, Q8=b.
        // This keeps Q3/Q4 fixed and deforms Q9=b−λq, exactly as the
        // physical central-edge soft profile. The parameter is not Taylor t.
        let q = function!(GS.emr_vec, 6, W_.a___);
        let r = function!(GS.emr_vec, 5, W_.a___);
        let mut deformation = vec![
            Replacement::new(q.to_pattern(), &q * lambda),
            Replacement::new(r.to_pattern(), &r + (Atom::one() - lambda) * &q),
        ];
        for component in 1..=3 {
            let q = function!(GS.emr_mom, 6, GS.cind(component));
            let r = function!(GS.emr_mom, 5, GS.cind(component));
            deformation.push(Replacement::new(q.to_pattern(), &q * lambda));
            deformation.push(Replacement::new(
                r.to_pattern(),
                &r + (Atom::one() - lambda) * &q,
            ));
        }
        Self {
            routing,
            deformation,
            energies,
            lambda,
        }
    }

    fn apply(&self, atom: &Atom, algebra: &mut SoftEnergyAlgebra) -> Result<Atom> {
        let mut physical = atom
            .unwrap_function(GS.tree_denom_wrapper)
            .replace_multiple(&self.energies)
            .replace_multiple(&self.routing)
            .replace_multiple(&self.deformation);
        physical = physical
            .replace(function!(GS.emr_vec, W_.a_, GS.cind(0)))
            .with(0);
        for component in 1..=3 {
            physical = physical
                .replace(function!(GS.emr_vec, W_.a_, GS.cind(component)))
                .with(function!(GS.emr_mom, W_.a_, GS.cind(component)));
        }
        let mut unfinished = None;
        physical.visitor(&mut |view| {
            if let AtomView::Fun(call) = view {
                if call.get_symbol() == GS.m_uv_vacuum
                    || call.get_symbol() == *crate::uv::approx::direct_3d::LOCAL_3D_MASS_SCOPE
                {
                    unfinished = Some(view.to_owned());
                }
                if [GS.emr_mom, GS.emr_vec].contains(&call.get_symbol()) {
                    let edge = i64::try_from(call.get(0)).ok();
                    let routed = edge.is_some_and(|edge| [0, 5, 6, 8].contains(&edge));
                    let spatial = call.get_symbol() == GS.emr_vec
                        || edge == Some(0)
                        || (1..=3).any(|component| call.get(1) == GS.cind(component).as_view());
                    if !routed || !spatial {
                        unfinished = Some(view.to_owned());
                    }
                }
            }
            unfinished.is_none()
        });
        ensure!(
            unfinished.is_none(),
            "noncanonical final physical leaf: {unfinished:?}"
        );
        algebra.normalize_energies(&physical, self.lambda)
    }
}

/// Split a polynomial in one selected family without distributing products of
/// sums of that family. The same traversal handles known Laurent monomials.
fn split_linear(
    atom: &Atom,
    contains: &impl Fn(AtomView<'_>) -> bool,
    literal: &impl Fn(AtomView<'_>) -> bool,
) -> Result<BTreeMap<Atom, Atom>> {
    fn visit(
        atom: AtomView<'_>,
        scalar: Atom,
        contains: &impl Fn(AtomView<'_>) -> bool,
        literal: &impl Fn(AtomView<'_>) -> bool,
        groups: &mut BTreeMap<Atom, Vec<Atom>>,
    ) -> Result<()> {
        if !contains(atom) {
            groups.entry(Atom::one()).or_default().push(scalar * atom);
        } else if literal(atom) {
            groups.entry(atom.to_owned()).or_default().push(scalar);
        } else {
            match atom {
                AtomView::Add(sum) => {
                    for term in sum {
                        visit(term, scalar.clone(), contains, literal, groups)?;
                    }
                }
                AtomView::Mul(product) => {
                    let (selected, spectators): (Vec<_>, Vec<_>) =
                        product.iter().partition(|part| contains(*part));
                    let scalar = scalar * Atom::mul_many(spectators);
                    if selected.iter().all(|part| literal(*part)) {
                        groups
                            .entry(Atom::mul_many(selected))
                            .or_default()
                            .push(scalar);
                    } else {
                        ensure!(
                            selected.len() == 1,
                            "unproven: refusing to distribute products of numerator sums"
                        );
                        visit(selected[0], scalar, contains, literal, groups)?;
                    }
                }
                _ => return Err(eyre!("unproven: unsupported selected-family shape")),
            }
        }
        ensure!(
            groups.len() <= MAX_GROUPS,
            "unproven: scalar grouping budget"
        );
        Ok(())
    }
    ensure!(
        atom.as_view().get_byte_size() <= MAX_COEFFICIENT_BYTES,
        "unproven: Laurent coefficient byte budget"
    );
    let mut groups = BTreeMap::new();
    visit(atom.as_view(), Atom::one(), contains, literal, &mut groups)?;
    Ok(groups
        .into_iter()
        .map(|(key, terms)| (key, Atom::add_many(terms)))
        .collect())
}

fn power_of(atom: AtomView<'_>, lambda: Symbol) -> Option<i64> {
    if atom.is_one() {
        return Some(0);
    }
    if atom == Atom::var(lambda).as_view() {
        return Some(1);
    }
    if let AtomView::Pow(power) = atom
        && power.get_base() == Atom::var(lambda).as_view()
    {
        return i64::try_from(power.get_exp()).ok();
    }
    None
}

fn laurent_coefficients(atom: &Atom, lambda: Symbol) -> Result<BTreeMap<i64, Atom>> {
    split_linear(atom, &|part| part.contains_symbol(lambda), &|part| {
        power_of(part, lambda).is_some()
    })?
    .into_iter()
    .map(|(power, coefficient)| {
        let power =
            power_of(power.as_view(), lambda).ok_or_else(|| eyre!("non-Laurent soft monomial"))?;
        ensure!(
            !coefficient.contains_symbol(lambda),
            "soft parameter hidden in coefficient"
        );
        Ok((power, coefficient))
    })
    .collect()
}

fn numerator_groups(atom: &Atom) -> Result<BTreeMap<Atom, Atom>> {
    let family = symbol!("gammalooprs::uv::numerator_family");
    let call = |part: AtomView<'_>| matches!(part, AtomView::Fun(f) if f.get_symbol() == family);
    split_linear(atom, &|part| part.contains_symbol(family), &|part| {
        call(part)
            || matches!(part, AtomView::Pow(power)
            if call(power.get_base()) && i64::try_from(power.get_exp()).is_ok_and(|p| (0..=64).contains(&p)))
    })
}

pub(super) fn check(input: &PhysicalOverlapInput) -> Result<()> {
    let complete = input.sum_except(None)?;
    debug_tags!(#uv, #overlap_poles;
        stage = "overlap_poles_signed_input", orientations = input.orientation_count,
        forest_nodes = input.nodes.len(), retained_definitions = complete.numerators().len(),
        "Checked the complete production definition store before physical soft expansion");
    drop(complete);
    let map = PhysicalSoftMap::new(input);
    let mut algebra = SoftEnergyAlgebra::default();
    let mut tensors = NumeratorJetCache::default();
    let family = symbol!("gammalooprs::uv::numerator_family");
    // Each row remains attached to its exact forest node. Removing the control
    // uses these same computed coefficients, not a second generation path.
    let mut rows = BTreeMap::<(i64, Atom), Vec<(usize, Atom)>>::new();
    let mut numerator_bodies = BTreeMap::<Atom, Atom>::new();
    for node in &input.nodes {
        debug_tags!(#uv, #overlap_poles;
            stage = "overlap_poles_node_series_start", node_index = node.node_index,
            node_key = %node.node_key, forest_sign = node.forest_sign,
            "Construct finite physical soft jets of one complete signed forest node");
        let terms = input.node_integrands(node.node_index)?;
        debug_tags!(#uv, #overlap_poles;
            stage = "overlap_poles_node_input", node_index = node.node_index,
            file.roots = ?terms.iter().collect_vec(),
            file.numerators = ?terms.numerators(),
            "Replayable physical-map input with separately retained numerator bodies");
        let terms = terms
            .fallible_map(|atom| map.apply(atom, &mut algebra))?
            .map_numerators(|atom| map.apply(atom, &mut algebra))?;
        let branches = DirectResidueBranches::production(OrientationID::from(0), terms)?;
        let jets = branches
            .series_preserving_numerators_exact(
                map.lambda,
                Atom::Zero.as_view(),
                LAST_FORBIDDEN_POWER,
                DirectResidueBranches::numerator_scope().1,
            )?
            .materialize(false)?;
        let replacements = jets
            .numerators()
            .iter()
            .map(|entry| entry.replacement())
            .collect_vec();
        let mut count = 0usize;
        for (_, root) in jets.iter() {
            for (power, coefficient) in laurent_coefficients(root, map.lambda)? {
                if coefficient.is_zero() {
                    continue;
                }
                ensure!(
                    power <= LAST_FORBIDDEN_POWER,
                    "unexpected untruncated physical jet"
                );
                for (monomial, scalar) in numerator_groups(&coefficient)? {
                    let body = monomial.replace_multiple(&replacements);
                    ensure!(
                        !body.contains_symbol(family) && !body.contains_symbol(map.lambda),
                        "unresolved physical numerator coefficient"
                    );
                    let canonical = tensors.canonicalize(&body)?;
                    if canonical.is_zero {
                        continue;
                    }
                    numerator_bodies
                        .entry(canonical.key.clone())
                        .or_insert(body / &canonical.scalar_prefactor);
                    rows.entry((power, canonical.key))
                        .or_default()
                        .push((node.node_index, scalar * canonical.scalar_prefactor));
                    count += 1;
                }
            }
        }
        debug_tags!(#uv, #overlap_poles;
            stage = "overlap_poles_node_series_done", node_index = node.node_index,
            scalar_groups = count, cumulative_tensor_groups = rows.len(),
            "Retained tensor coefficients canonicalized without graph-numerator expansion");
    }
    ensure!(
        !rows.is_empty(),
        "vacuous overlap certificate: no singular contributions"
    );
    let mut certified = 0usize;
    let mut unresolved = 0usize;
    let mut control_unresolved = 0usize;
    let mut zero_rows = BTreeMap::<i64, usize>::new();
    let mut control_rows = BTreeMap::<i64, Vec<(Atom, Atom)>>::new();
    for ((power, numerator), contributions) in rows {
        let sum = Atom::add_many(contributions.iter().map(|(_, coefficient)| coefficient));
        let omitted = Atom::add_many(
            contributions
                .iter()
                .filter(|(node, _)| *node != input.negative_control_node)
                .map(|(_, coefficient)| coefficient),
        );
        control_rows
            .entry(power)
            .or_default()
            .push((omitted.clone(), numerator_bodies[&numerator].clone()));
        let proof = algebra.scalar_zero(&sum)?;
        let full_zero = matches!(proof, EnergyZeroOutcome::Zero(_));
        match &proof {
            EnergyZeroOutcome::Zero(witness) => {
                certified += 1;
                *zero_rows.entry(power).or_default() += 1;
                debug_tags!(#uv, #overlap_poles;
                    stage = "overlap_poles_coefficient_certificate", power,
                    file.numerator_key = %numerator,
                    file.contributions = ?contributions,
                    file.witness = ?witness,
                    "Exact zero coefficient of one unchanged physical numerator tensor");
            }
            _ => {
                unresolved += 1;
                debug_tags!(#uv, #overlap_poles;
                    stage = "overlap_poles_unproven", power,
                    file.numerator_key = %numerator,
                    file.numerator_body = %numerator_bodies[&numerator],
                    file.contributions = ?contributions, file.outcome = ?proof,
                    "Unresolved tensor/energy relation; not evidence of a physical pole");
            }
        }
        if full_zero && !matches!(algebra.scalar_zero(&omitted)?, EnergyZeroOutcome::Zero(_)) {
            control_unresolved += 1;
        }
    }
    // Include every tensor group at a given Laurent power before deciding
    // nonzero. An unresolved individual group is never a control witness.
    let mut negative_control_witness = false;
    for (point_index, point) in ExactControlPoint::points().into_iter().enumerate() {
        let at_point = (|| -> Result<BTreeMap<i64, Atom>> {
            let mut values = BTreeMap::new();
            for (power, groups) in &control_rows {
                let mut terms = BTreeMap::<i64, Atom>::new();
                for (scalar, body) in groups {
                    let scalar = point.apply(scalar)?;
                    let body = point.apply(body)?;
                    let tensor = tensors.canonicalize(&body)?;
                    ensure!(
                        tensor.is_zero || tensor.key.is_one(),
                        "control numerator remains symbolic or has open tensor slots"
                    );
                    let value = scalar * tensor.scalar_prefactor;
                    let (coefficient, pi_power) = ExactControlPoint::rational_pi_monomial(&value)?;
                    *terms.entry(pi_power).or_insert(Atom::Zero) += coefficient;
                }
                terms.retain(|_, coefficient| !coefficient.is_zero());
                // Keep the known nonzero measure factor π^n explicit. Only a
                // single surviving power is accepted; π is never assigned 1.
                ensure!(
                    terms.len() <= 1,
                    "control sum has unresolved mixed powers of π"
                );
                let value = terms
                    .into_iter()
                    .next()
                    .map_or(Atom::Zero, |(power, coefficient)| {
                        coefficient * Atom::var(Symbol::PI).pow(power)
                    });
                values.insert(*power, value);
            }
            Ok(values)
        })();
        match at_point {
            Ok(values) => {
                if let Some((power, value)) = values.iter().find(|(_, value)| !value.is_zero()) {
                    negative_control_witness = true;
                    debug_tags!(#uv, #overlap_poles;
                        stage = "overlap_poles_negative_control", point_index, power,
                        negative_control_node = input.negative_control_node,
                        file.value = %value, file.all_coefficients = ?values,
                        "Omitting the selected overlap gives an exactly nonzero physical pole coefficient");
                    break;
                }
                debug_tags!(#uv, #overlap_poles;
                    stage = "overlap_poles_control_zero_point", point_index,
                    "This exact point does not witness the omitted-overlap pole");
            }
            Err(error) => {
                debug_tags!(#uv, #overlap_poles;
                    stage = "overlap_poles_control_unsupported", point_index, file.error = %error,
                    "Control point is singular or unsupported; no nonzero claim");
            }
        }
    }
    debug_tags!(#uv, #overlap_poles;
        stage = "overlap_poles_summary", orientations = input.orientation_count,
        forest_nodes = input.nodes.len(), last_forbidden_power = LAST_FORBIDDEN_POWER,
        certified, unresolved, exact_zero_groups_by_power = ?zero_rows,
        negative_control_node = input.negative_control_node, control_unresolved,
        negative_control_witness, generic_soft_pole_certificate = unresolved == 0,
        "Physical soft certificate coverage; a control remainder alone is not a physical nonzero witness");
    ensure!(
        unresolved == 0,
        "overlap pole certificate incomplete: {unresolved} unresolved tensor/energy groups"
    );
    ensure!(
        negative_control_witness,
        "negative-control exact physical witness unavailable ({control_unresolved} unresolved control groups)"
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::parse_lit;

    #[test]
    fn physical_laurent_collection_preserves_tiny_terms() {
        let t = symbol!("overlap_collection_test::t");
        let atom = (parse_lit!(1 / 10 ^ 30) + parse_lit!(a)) / Atom::var(t).pow(3)
            + parse_lit!(b) / Atom::var(t).pow(2);
        let groups = laurent_coefficients(&atom, t).unwrap();
        assert_eq!(groups[&-3], parse_lit!(a + 1 / 10 ^ 30));
        assert_eq!(groups[&-2], parse_lit!(b));
    }

    #[test]
    fn family_collection_refuses_products_of_numerator_sums() {
        let family = symbol!("gammalooprs::uv::numerator_family");
        let a = family.call_args([0]);
        let b = family.call_args([1]);
        let c = family.call_args([2]);
        assert!(numerator_groups(&((&a + &b) * (&b + &c))).is_err());
        let groups = numerator_groups(&(parse_lit!(x) * (&a + &b))).unwrap();
        assert_eq!(groups[&a], parse_lit!(x));
        assert_eq!(groups[&b], parse_lit!(x));
    }
}
