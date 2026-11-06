//! Published five-loop p-integrals from arXiv:2104.08272v1.
//!
//! Master numbers, vertex labels, and Laurent coefficients below come from
//! https://arxiv.org/src/2104.08272v1/anc/MasterExpansions.m, whose SHA256 is
//! 388e17209e9e97eabf958fb87c310ef918dbd680b261fdc99744c40fc20cec0a.
//! Entries have external legs (v[1],e[1]) and (v[2],e[2]); all listed internal
//! propagators have power one. Equation (11) uses the G-scheme at Euclidean
//! Q²=1. The ten propagators give an even Wick-rotation sign.

use std::{collections::BTreeMap, fmt::Write};

use eyre::Result;
use spenso::{
    network::tags::SPENSO_TAG,
    structure::representation::{Minkowski, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    function,
    transcendental::zeta,
};
use vakint::Vakint;

use crate::{
    debug_tags,
    graph::{Graph, parse::IntoFinalizedRuntimeGraph},
    initialisation::test_initialise,
    model::Model,
    utils::{GS, load_generic_model},
};

use super::{RStarPlan, tests::canonical_tensor};

struct PublishedMaster {
    number: usize,
    edges: [(usize, usize); 10],
    marked: usize,
    proper_spinneys: usize,
}

impl PublishedMaster {
    fn graph(&self, model: &Model) -> Result<Graph> {
        let largest_vertex = self.edges.iter().flat_map(|&(a, b)| [a, b]).max().unwrap();
        let mut components: Vec<_> = (0..=largest_vertex).collect();
        let mut slots = BTreeMap::from([(1, 1), (2, 1)]);
        let mut loops = 0;
        let mut dot = format!(
            "digraph master_{} {{\n\
             graph [projector=1]; node [num=1];\n\
             edge [particle=\"scalar_1\" num=1]; e [style=invis];\n\
             e -> v1:0 [id=10 sink=\"{{ufo_order:0}}\"];\n\
             v2:1 -> e [id=11 source=\"{{ufo_order:0}}\"];\n",
            self.number,
        );
        for (id, &(a, b)) in self.edges.iter().enumerate() {
            let source_slot = *slots.entry(a).or_default();
            *slots.get_mut(&a).unwrap() += 1;
            let sink_slot = *slots.entry(b).or_default();
            *slots.get_mut(&b).unwrap() += 1;
            let basis = if components[a] == components[b] {
                let label = format!("lmb_id={loops}");
                loops += 1;
                label
            } else {
                let from = components[b];
                let into = components[a];
                for component in &mut components {
                    if *component == from {
                        *component = into;
                    }
                }
                String::new()
            };
            writeln!(
                dot,
                "v{a} -> v{b} [id={id} {basis} source=\"{{ufo_order:{source_slot}}}\" sink=\"{{ufo_order:{sink_slot}}}\"];"
            )?;
        }
        assert_eq!(loops, 5);
        dot.push('}');
        let mut graph: Graph = dot.into_finalized_runtime_graph(model)?;
        let edges: Vec<_> = graph.iter_edges().map(|(_, edge, _)| edge).collect();
        for edge in edges {
            let data = &mut graph.underlying[edge];
            data.particle = data.particle.clone().override_mass(Some(Atom::Zero));
        }
        Ok(graph)
    }

    fn assert_pole_part(&self, expected: Atom) -> Result<()> {
        let model = load_generic_model("scalars");
        let graph = self.graph(&model)?;
        let mut settings = super::tests::settings();
        settings.vakint.normalization = "GScheme".into();
        let marked = graph
            .iter_edges()
            .filter_map(|(pair, edge, _)| pair.is_paired().then_some(edge))
            .nth(self.marked)
            .unwrap();
        let plan = RStarPlan::new(&graph, &model, marked, &settings)?;
        assert_eq!(plan.taylor_order, 0);
        assert_eq!(plan.ordinary_uv.len(), self.proper_spinneys + 1);
        let result = plan.pole_part(&graph, &model, &Vakint::new()?, &settings)?;
        let external = graph
            .iter_edges()
            .find_map(|(pair, edge, _)| (!pair.is_paired()).then_some(edge))
            .unwrap();
        let representation = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
        let momentum = function!(GS.emr_mom, usize::from(external), representation);
        let squared = function!(SPENSO_TAG.dot, &momentum, &momentum);
        let actual = canonical_tensor(&result)?
            .replace(squared.to_pattern())
            .with(Atom::num(-1))
            .replace(GS.mu_r_sq)
            .with(Atom::one())
            .expand();
        debug_tags!(#uv, #integrated, #rstar;
            stage = "five_loop_published_master_comparison",
            master = self.number,
            log.actual = actual,
            log.expected = expected,
            "Five-loop bare pole part compared with arXiv:2104.08272"
        );
        let difference = (&actual - &expected).expand().together();
        assert!(
            difference.is_zero(),
            "master {}: reconstructed {actual}; published {expected}; difference {difference}",
            self.number,
        );
        Ok(())
    }
}

#[test]
fn five_loop_master_159_primitive_poles_match_publication() -> Result<()> {
    test_initialise()?;
    PublishedMaster {
        number: 159,
        edges: [
            (1, 6),
            (1, 2),
            (3, 6),
            (1, 3),
            (4, 6),
            (3, 4),
            (2, 4),
            (5, 6),
            (3, 5),
            (2, 5),
        ],
        marked: 0,
        proper_spinneys: 0,
    }
    .assert_pole_part(
        Atom::num(36) * zeta().call(3).pow(2) / (Atom::num(5) * Atom::var(GS.dim_epsilon)),
    )
}

#[test]
fn five_loop_master_169_primitive_poles_match_publication() -> Result<()> {
    test_initialise()?;
    PublishedMaster {
        number: 169,
        edges: [
            (1, 10),
            (2, 4),
            (3, 4),
            (3, 10),
            (4, 5),
            (1, 4),
            (2, 5),
            (5, 10),
            (2, 3),
            (1, 3),
        ],
        marked: 0,
        proper_spinneys: 0,
    }
    .assert_pole_part(Atom::num(441) * zeta().call(7) / (Atom::num(40) * Atom::var(GS.dim_epsilon)))
}

#[test]
fn five_loop_master_76_double_pole_matches_publication() -> Result<()> {
    test_initialise()?;
    let epsilon = Atom::var(GS.dim_epsilon);
    let expected = Atom::num(4) * zeta().call(5) / epsilon.pow(2)
        + (Atom::num(2) * Atom::var(Symbol::PI).pow(6) / Atom::num(189)
            - Atom::num(34) * zeta().call(3).pow(2) / Atom::num(5))
            / epsilon;
    PublishedMaster {
        number: 76,
        edges: [
            (1, 9),
            (2, 6),
            (3, 6),
            (1, 3),
            (3, 4),
            (4, 6),
            (1, 4),
            (3, 9),
            (2, 6),
            (2, 9),
        ],
        marked: 0,
        proper_spinneys: 1,
    }
    .assert_pole_part(expected)
}

#[test]
fn five_loop_master_203_triple_pole_matches_publication() -> Result<()> {
    test_initialise()?;
    let epsilon = Atom::var(GS.dim_epsilon);
    let pi_four = Atom::var(Symbol::PI).pow(4);
    let zeta_three = zeta().call(3);
    let expected = Atom::num(6) * &zeta_three / (Atom::num(5) * epsilon.pow(3))
        + (&pi_four / Atom::num(50) + Atom::num(12) * &zeta_three / Atom::num(5)) / epsilon.pow(2)
        + (pi_four / Atom::num(25) + Atom::num(96) * zeta_three / Atom::num(5)
            - Atom::num(23) * zeta().call(5) / Atom::num(5))
            / epsilon;
    PublishedMaster {
        number: 203,
        edges: [
            (1, 3),
            (2, 8),
            (3, 5),
            (3, 8),
            (4, 5),
            (4, 8),
            (1, 4),
            (2, 5),
            (2, 5),
            (1, 3),
        ],
        marked: 1,
        proper_spinneys: 3,
    }
    .assert_pole_part(expected)
}
