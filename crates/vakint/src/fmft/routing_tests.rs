use super::*;
use symbolica::id::Replacement;

#[test]
fn every_registered_four_loop_route_preserves_fmft_propagator_squares() {
    let vakint = Vakint::new().unwrap();
    // These four independent reference momenta reproduce all scalar-product
    // identities in FMFT's sp2den procedure. No FORM execution is necessary.
    let fmft_momenta = [
        "r1", "r2", "r3", "r4", "r1-r4", "r2-r4", "r3-r4", "r1-r2", "r1-r3", "r1-r2-r3",
    ]
    .map(|expression| vk_parse!(expression).unwrap());
    let mut checked_classes = 0;
    let mut checked_edges = 0;
    for (parent, slots) in [("I4L_H", 9), ("I4L_X", 9), ("I4L_BMW", 8), ("I4L_FG", 8)] {
        let map = FMFT::oriented_edge_map(parent).unwrap();
        // FG owns all registered contracted graph classes. Query their actual
        // short forms rather than duplicating the canonical routing table.
        let masks = if parent == "I4L_FG" { 1 << slots } else { 1 };
        for mask in 0..masks {
            let contractions = (0..slots)
                .filter(|slot| mask & (1 << slot) != 0)
                .map(|slot| (slot + 1).to_string())
                .collect::<Vec<_>>();
            let name = if contractions.is_empty() {
                parent.to_owned()
            } else {
                format!("{parent}_pinch_{}", contractions.join("_"))
            };
            let powers = (0..slots)
                .map(|slot| if mask & (1 << slot) == 0 { "1" } else { "0" })
                .collect::<Vec<_>>()
                .join(",");
            let input = vk_parse!(&format!("topo({name}(muvsq,{powers}))")).unwrap();
            let Some(specs) = vakint
                .topologies
                .match_topologies_to_user_input(input.as_view(), false)
                .unwrap()
            else {
                continue;
            };
            let momenta = specs.get_propagator_property_list("q_");
            let substitutions = (1..=4)
                .map(|axis| {
                    let loop_momentum = function!(S.k, Atom::num(axis));
                    let (slot, _) = momenta
                        .iter()
                        .find(|(_, momentum)| *momentum == &loop_momentum)
                        .expect("registered basis momentum must have an edge");
                    let signed_id = map[slot - 1];
                    let image = Atom::num(signed_id.signum())
                        * &fmft_momenta[signed_id.unsigned_abs() as usize - 1];
                    Replacement::new(loop_momentum.to_pattern(), image.to_pattern())
                })
                .collect::<Vec<_>>();
            for (slot, momentum) in momenta {
                let mapped = momentum.replace_multiple(&substitutions);
                let expected = &fmft_momenta[map[slot - 1].unsigned_abs() as usize - 1];
                assert!(
                    (mapped.pow(Atom::num(2)) - expected.pow(Atom::num(2)))
                        .expand()
                        .is_zero(),
                    "{name} edge {slot}: {momentum} was mapped to {mapped}, expected ±{expected}"
                );
                checked_edges += 1;
            }
            checked_classes += 1;
        }
    }
    assert_eq!(checked_classes, 19);
    assert_eq!(checked_edges, 123);
}
