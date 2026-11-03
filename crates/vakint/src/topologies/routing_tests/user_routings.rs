//! Every matched witness must rewrite the user's loop momenta so that each
//! input line momentum becomes (up to sign) the canonical momentum of the line
//! it was identified with. Only then is the numerator substitution the change
//! of integration variables that maps the input integrand onto the canonical
//! one. Every registered topology (one to five loops) is written in routings a
//! user may choose: lines reversed with negated momenta, and loop momenta that
//! are combinations of line momenta instead of line momenta themselves.

use std::collections::BTreeSet;

use super::*;
use crate::{Vakint, get_integer_from_atom, get_prop_with_id};
use symbolica::id::Replacement;

/// `k(i) -> k(i)+k(i+1)` for `i < L`, `k(L) -> -k(L)` (determinant `-1`):
/// apart from `-k(L)`, no loop momentum is a line momentum any more.
fn unimodular_basis_change(loops: usize) -> Vec<Replacement> {
    (1..=loops)
        .map(|axis| {
            let image = if axis < loops {
                S.k.call(axis) + S.k.call(axis + 1)
            } else {
                -S.k.call(axis)
            };
            Replacement::new(S.k.call(axis).to_pattern(), image.to_pattern())
        })
        .collect()
}

/// The canonical expression with unit powers and mass `muvsq`, its loop momenta
/// optionally rewritten and its odd-numbered lines optionally reversed
/// (`edge(b,a)` with the negated momentum).
fn user_routing(topology: &Topology, change_basis: bool, reverse_odd_lines: bool) -> Atom {
    let integral = topology.get_integral();
    let mut expression = integral
        .canonical_expression
        .clone()
        .unwrap()
        .replace(vk_parse!("msq(i_)").unwrap().to_pattern())
        .with(vk_parse!("muvsq").unwrap().to_pattern())
        .replace(vk_parse!("pow(i_)").unwrap().to_pattern())
        .with(Atom::num(1).to_pattern());
    if change_basis {
        let basis = unimodular_basis_change(integral.n_loops);
        expression = expression
            .replace(vk_parse!("prop(id_,e_,q_,m_,p_)").unwrap().to_pattern())
            .with_map(move |m| {
                let q = m.get(vk_symbol!("q_")).unwrap().to_atom();
                S.prop.call((
                    m.get(vk_symbol!("id_")).unwrap().to_atom(),
                    m.get(vk_symbol!("e_")).unwrap().to_atom(),
                    q.replace_multiple(&basis).expand(),
                    m.get(vk_symbol!("m_")).unwrap().to_atom(),
                    m.get(vk_symbol!("p_")).unwrap().to_atom(),
                ))
            });
    }
    if reverse_odd_lines {
        expression = expression
            .replace(
                vk_parse!("prop(id_,edge(a_,b_),q_,m_,p_)")
                    .unwrap()
                    .to_pattern(),
            )
            .with_map(|m| {
                let id = m.get(vk_symbol!("id_")).unwrap().to_atom();
                let (a, b) = (
                    m.get(vk_symbol!("a_")).unwrap().to_atom(),
                    m.get(vk_symbol!("b_")).unwrap().to_atom(),
                );
                let q = m.get(vk_symbol!("q_")).unwrap().to_atom();
                let odd = get_integer_from_atom(id.as_view()).unwrap() % 2 == 1;
                S.prop.call((
                    id,
                    if odd {
                        S.edge.call((b, a))
                    } else {
                        S.edge.call((a, b))
                    },
                    if odd { (-q).expand() } else { q },
                    m.get(vk_symbol!("m_")).unwrap().to_atom(),
                    m.get(vk_symbol!("p_")).unwrap().to_atom(),
                ))
            });
    }
    expression
}

/// Match `input` and check that the numerator substitution maps every input
/// line momentum onto plus or minus the canonical momentum of its image line.
fn assert_witness_maps_every_line(topologies: &Topologies, input: &Atom) {
    let rules = topologies
        .match_topologies_to_user_input(input.as_view(), false)
        .unwrap_or_else(|error| panic!("matching failed for {input}: {error}"))
        .unwrap_or_else(|| panic!("no topology matched {input}"));
    let canonical = rules
        .canonical_topology
        .get_integral()
        .canonical_expression
        .clone()
        .unwrap();
    let substitution = rules
        .numerator_substitutions
        .iter()
        .map(|(source, target)| Replacement::new(source.to_pattern(), target.to_pattern()))
        .collect::<Vec<_>>();
    let input_lines = Topologies::count_propagators_in_integral(input.as_view());
    assert_eq!(
        rules.edge_ids_canonical_to_input_map.len(),
        input_lines,
        "not every line of {input} was matched"
    );
    for (canonical_id, input_id) in &rules.edge_ids_canonical_to_input_map {
        let q_input = get_prop_with_id(input.as_view(), *input_id).unwrap()[&vk_symbol!("q_")]
            .replace_multiple(&substitution)
            .expand();
        let q_canonical = get_prop_with_id(canonical.as_view(), *canonical_id).unwrap()
            [&vk_symbol!("q_")]
            .clone();
        assert!(
            (&q_input - &q_canonical).expand().is_zero()
                || (&q_input + &q_canonical).expand().is_zero(),
            "{}: input line {input_id} becomes {q_input}, but it was matched to canonical line \
             {canonical_id} with momentum {q_canonical}\ninput: {input}\nwitness: {rules}",
            rules.canonical_topology.get_integral().name
        );
    }
}

fn assert_all_user_routings(topology: &Topology, topologies: &Topologies) {
    for change_basis in [false, true] {
        for reverse_odd_lines in [false, true] {
            let input = user_routing(topology, change_basis, reverse_odd_lines);
            assert_witness_maps_every_line(topologies, &input);
        }
    }
}

/// Every one- to three-loop topology, the four-loop parents and every eighth
/// automatic FG pinch (all of them pass too, but take minutes).
#[test]
fn one_to_four_loop_witnesses_map_every_user_routing_onto_the_canonical_lines() {
    let vakint = Vakint::new().unwrap();
    let mut loop_counts = BTreeSet::new();
    let mut four_loop_pinches = 0;
    for topology in vakint.topologies.iter() {
        let integral = topology.get_integral();
        if matches!(topology, Topology::Unknown(_)) {
            continue;
        }
        if integral.n_loops == 4 && integral.name.contains("_pinch_") {
            four_loop_pinches += 1;
            if four_loop_pinches % 8 != 1 {
                continue;
            }
        }
        loop_counts.insert(integral.n_loops);
        if ["I4L_H", "I4L_X"].contains(&integral.name.as_str()) {
            // Without a perfect match the legacy matcher enumerates every
            // isomorphism (about 30 s per input): check the hardest routing only.
            assert_witness_maps_every_line(&vakint.topologies, &user_routing(topology, true, true));
        } else {
            assert_all_user_routings(topology, &vakint.topologies);
        }
    }
    assert_eq!(loop_counts, BTreeSet::from([1, 2, 3, 4]));
}

#[test]
fn five_loop_witnesses_map_every_user_routing_onto_the_canonical_lines() {
    let vakint = Vakint::new().unwrap();
    vakint
        .topologies
        .prepare_five_loop_catalogue_for_tests()
        .unwrap();
    // Parents (no deleted slot) and one contracted and one factorized sector.
    let census = &five_loop::data::CENSUS;
    let selected = census
        .iter()
        .enumerate()
        .filter(|(_, witness)| witness.deleted_slots.is_empty())
        .map(|(index, _)| index)
        .chain(
            census
                .iter()
                .position(|witness| !witness.deleted_slots.is_empty()),
        )
        .chain(census.iter().position(|witness| witness.factorized))
        .collect::<BTreeSet<_>>();
    let five_loop = vakint
        .topologies
        .iter()
        .filter(|topology| matches!(topology, Topology::FiveLoop(_)))
        .collect::<Vec<_>>();
    assert_eq!(five_loop.len(), census.len());
    for index in selected {
        assert_all_user_routings(five_loop[index], &vakint.topologies);
    }
}
