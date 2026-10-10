//! Search census of the reference generation workloads: how many states each search stage keeps,
//! and how many of them lead nowhere.
//!
//! [`Generator::generate`] searches in parallel. The parallel frontier interleaves subtrees and
//! starts split subtrees at children it has already filtered, so stage counts taken there
//! overcount. The census therefore runs the same generator sequentially and counts through its
//! filter, which sees every kept search state in depth-first order. Vertices keep their numbers
//! in every extension of a stage-1 state and stage 2 types one edge per step, so the parent of a
//! state is the nearest earlier state on the depth-first path that it extends.
//!
//! The counts are deterministic. `FEYNKIT_ORACLE_RECORD=1` prints them instead of comparing
//! them; heavy workloads run with `just generation-oracle-heavy`.

use std::sync::{Arc, Mutex};

use linnet::half_edge::generate::Stage;

use super::{
    Generator, ResolvedProcess,
    cases::Case,
    enumeration::{GenerationView, Pruner, VertexSignatures},
};
use crate::{GenerationControl, GenerationError};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct SearchCensus {
    /// Kept complete topologies, and those without any diagram.
    topologies: usize,
    wasted_topologies: usize,
    /// Kept partial topologies; dead ones lead to no kept complete topology, useless ones (dead
    /// ones included) to no diagram.
    partial_topologies: usize,
    dead_partial_topologies: usize,
    useless_partial_topologies: usize,
    /// Kept partial field assignments (complete labellings included); dead ones lead to no
    /// complete labelling, useless ones (dead ones included) to no diagram. Assignments that
    /// linnet's forward checking or budgets (coupling orders, fermion loops) drop never reach the
    /// filter and are not counted.
    dressings: usize,
    dead_dressings: usize,
    useless_dressings: usize,
    diagrams: usize,
    /// States the process filters rejected per stage: topology, fields, complete.
    rejected: [usize; 3],
}

/// A kept state on the depth-first path, with what was found below it.
struct Open<K> {
    key: K,
    complete: bool,
    diagram: bool,
}

/// A stage-1 state: per vertex its degree, self-loops and sorted legs, and the edge
/// multiplicities between earlier vertices.
#[derive(PartialEq, Eq)]
struct Topology {
    vertices: Vec<(usize, usize, Vec<String>)>,
    multiplicities: Vec<Vec<usize>>,
    complete: bool,
}

impl Topology {
    fn of(view: &GenerationView<'_>) -> Self {
        let n = view.n_vertices();
        Topology {
            vertices: (0..n)
                .map(|v| {
                    let mut legs: Vec<_> = view.legs(v).map(|leg| format!("{leg:?}")).collect();
                    legs.sort();
                    (view.degree(v), view.self_loops(v), legs)
                })
                .collect(),
            multiplicities: (0..n)
                .map(|v| (0..v).map(|u| view.multiplicity(u, v)).collect())
                .collect(),
            complete: view.is_topology_complete(),
        }
    }

    /// `other` extends `self`: more vertices, or more labelled legs, on the same vertices,
    /// self-loops and edges.
    fn extended_by(&self, other: &Topology) -> bool {
        self != other
            && self.vertices.len() <= other.vertices.len()
            && self.multiplicities[..] == other.multiplicities[..self.vertices.len()]
            && self.vertices.iter().zip(&other.vertices).all(
                |((degree, loops, legs), (other_degree, other_loops, other_legs))| {
                    degree == other_degree
                        && loops == other_loops
                        && legs.iter().all(|leg| other_legs.contains(leg))
                },
            )
    }
}

#[derive(Default)]
struct Tracker {
    census: SearchCensus,
    topologies: Vec<Open<Topology>>,
    /// Typed edges of the stage-2 states on the path of the current complete topology.
    dressings: Vec<Open<usize>>,
}

impl Tracker {
    fn record(&mut self, view: &GenerationView<'_>, kept: bool) {
        let stage = view.stage();
        if !kept {
            let index = match stage {
                Stage::Topology => 0,
                Stage::Fields => 1,
                Stage::Complete => 2,
            };
            self.census.rejected[index] += 1;
            return;
        }
        match stage {
            Stage::Topology => {
                self.close_dressings(0);
                let topology = Topology::of(view);
                while self
                    .topologies
                    .last()
                    .is_some_and(|open| !open.key.extended_by(&topology))
                {
                    self.close_topology();
                }
                let complete = topology.complete;
                if complete {
                    self.census.topologies += 1;
                } else {
                    self.census.partial_topologies += 1;
                }
                self.topologies.push(Open {
                    key: topology,
                    complete,
                    diagram: false,
                });
            }
            Stage::Fields => {
                let typed = view.edges().filter(|edge| edge.stub.is_some()).count();
                self.close_dressings(typed);
                self.census.dressings += 1;
                self.dressings.push(Open {
                    key: typed,
                    complete: typed == view.edges().count(),
                    diagram: false,
                });
            }
            Stage::Complete => {
                self.census.diagrams += 1;
                for open in &mut self.dressings {
                    open.diagram = true;
                }
                if let Some(topology) = self.topologies.last_mut() {
                    topology.diagram = true;
                }
            }
        }
    }

    /// Close the stage-2 states with at least `typed` typed edges.
    fn close_dressings(&mut self, typed: usize) {
        while self.dressings.last().is_some_and(|open| open.key >= typed) {
            let open = self.dressings.pop().unwrap();
            self.census.dead_dressings += usize::from(!open.complete);
            self.census.useless_dressings += usize::from(!open.diagram);
            if let Some(parent) = self.dressings.last_mut() {
                parent.complete |= open.complete;
                parent.diagram |= open.diagram;
            }
        }
    }

    fn close_topology(&mut self) {
        let open = self.topologies.pop().unwrap();
        if open.key.complete {
            self.census.wasted_topologies += usize::from(!open.diagram);
        } else {
            self.census.dead_partial_topologies += usize::from(!open.complete);
            self.census.useless_partial_topologies += usize::from(!open.diagram);
        }
        if let Some(parent) = self.topologies.last_mut() {
            parent.complete |= open.complete;
            parent.diagram |= open.diagram;
        }
    }

    fn finish(mut self) -> SearchCensus {
        self.close_dressings(0);
        while !self.topologies.is_empty() {
            self.close_topology();
        }
        self.census
    }
}

impl SearchCensus {
    /// The census of `case` with one thread, from the generator that [`Generator::generate`]
    /// builds for it, searched sequentially; the diagram and topology counts are checked
    /// against [`Generator::generate`] cancelled after enumeration.
    fn of(case: &Case) -> Result<Self, GenerationError> {
        let model = case.model();
        let options = case.options(&model, 1);
        let generator = Generator::new(model);
        let model = &generator.model;
        let process = case.process();
        process.validate()?;
        process.validate_in(model)?;
        let options = process
            .restrict_options(model, &options)?
            .resolve_selectors(model)?;
        generator.validate_options(case.kind, &options)?;
        process.validate_covariant_cut_filters(model, &options, case.kind)?;
        let process = process
            .clone()
            .with_final_state_alternatives(process.covariant_cut_states(model, case.kind)?)?;
        let resolved = ResolvedProcess::new(model, &process, case.kind, &options)?;
        let signatures = VertexSignatures::new(model, &options)?;
        let externals = resolved.externals(model)?;
        assert!(!externals.vetoed(model, &options));
        let pruner = Pruner::new(model, &options);
        let tracker = Arc::new(Mutex::new(Tracker::default()));
        let recorder = Arc::clone(&tracker);
        let search = signatures.generator(model, &externals, &options, move |view| {
            let kept = pruner.keeps(view);
            recorder.lock().unwrap().record(view, kept);
            kept
        })?;
        let diagrams = search.iter()?.count();
        drop(search);
        let census = Arc::into_inner(tracker)
            .unwrap()
            .into_inner()
            .unwrap()
            .finish();
        assert_eq!(census.diagrams, diagrams);

        let enumeration = case
            .options(model, 1)
            .progress(|progress| match progress.stage {
                "topologies" => GenerationControl::Continue,
                _ => GenerationControl::Cancel,
            });
        let report = generator
            .generate(&case.process(), &enumeration, case.kind)?
            .report;
        assert_eq!(
            (census.topologies, census.diagrams),
            (report.topology_count, report.interaction_assignment_count),
            "{}: the sequential search differs from generation",
            case.name
        );
        Ok(census)
    }
}

fn check(name: &str, expected: SearchCensus) {
    let census = SearchCensus::of(Case::named(name)).unwrap();
    if std::env::var_os("FEYNKIT_ORACLE_RECORD").is_some() {
        println!("{name}: {census:?}");
    } else {
        assert_eq!(census, expected, "{name}");
    }
}

#[test]
fn ee_to_four_leptons_at_two_loops() {
    check(
        "qed_ee4f_L2",
        SearchCensus {
            topologies: 7909,
            wasted_topologies: 5683,
            partial_topologies: 102805,
            dead_partial_topologies: 87076,
            useless_partial_topologies: 97981,
            dressings: 53951,
            dead_dressings: 13583,
            useless_dressings: 14573,
            diagrams: 8460,
            rejected: [0, 0, 0],
        },
    );
}

#[test]
fn standard_model_ee_to_mumu_at_two_loops() {
    check(
        "sm_eemm_L2",
        SearchCensus {
            topologies: 160,
            wasted_topologies: 62,
            partial_topologies: 1391,
            dead_partial_topologies: 1086,
            useless_partial_topologies: 1190,
            dressings: 175615,
            dead_dressings: 1680,
            useless_dressings: 19734,
            diagrams: 69697,
            rejected: [0, 0, 0],
        },
    );
}

#[test]
fn standard_model_uu_to_dd_with_qed_order_two_at_two_loops() {
    check(
        "sm_uudd_qed2_L2",
        SearchCensus {
            topologies: 230,
            wasted_topologies: 77,
            partial_topologies: 1698,
            dead_partial_topologies: 1325,
            useless_partial_topologies: 1423,
            dressings: 14463,
            dead_dressings: 1761,
            useless_dressings: 2743,
            diagrams: 3286,
            rejected: [0, 0, 0],
        },
    );
}

#[test]
fn gg_to_six_gluon_trees() {
    check(
        "ym_g8_L0",
        SearchCensus {
            topologies: 34300,
            wasted_topologies: 0,
            partial_topologies: 65417,
            dead_partial_topologies: 19,
            useless_partial_topologies: 19,
            dressings: 140735,
            dead_dressings: 0,
            useless_dressings: 0,
            diagrams: 34300,
            rejected: [0, 0, 0],
        },
    );
}

#[test]
#[ignore = "heavy oracle tier: run with `just generation-oracle-heavy`"]
fn gg_to_gg_at_three_loops() {
    check(
        "qcd_gggg_L3",
        SearchCensus {
            topologies: 14340,
            wasted_topologies: 0,
            partial_topologies: 35379,
            dead_partial_topologies: 6324,
            useless_partial_topologies: 6324,
            dressings: 3853181,
            dead_dressings: 369908,
            useless_dressings: 713475,
            diagrams: 686501,
            rejected: [0, 0, 0],
        },
    );
}

#[test]
#[ignore = "heavy oracle tier: run with `just generation-oracle-heavy`"]
fn gg_to_ggg_at_two_loops() {
    check(
        "qcd_ggggg_L2",
        SearchCensus {
            topologies: 11685,
            wasted_topologies: 0,
            partial_topologies: 24107,
            dead_partial_topologies: 1138,
            useless_partial_topologies: 1138,
            dressings: 1171205,
            dead_dressings: 152880,
            useless_dressings: 184835,
            diagrams: 213390,
            rejected: [0, 0, 0],
        },
    );
}

#[test]
#[ignore = "heavy oracle tier: run with `just generation-oracle-heavy`"]
fn uu_to_gg_at_three_loops() {
    check(
        "qcd_uugg_L3",
        SearchCensus {
            topologies: 7513,
            wasted_topologies: 2239,
            partial_topologies: 61781,
            dead_partial_topologies: 51158,
            useless_partial_topologies: 53959,
            dressings: 642646,
            dead_dressings: 49430,
            useless_dressings: 120974,
            diagrams: 114301,
            rejected: [0, 0, 0],
        },
    );
}

#[test]
#[ignore = "heavy oracle tier: run with `just generation-oracle-heavy`"]
fn photon_self_energy_at_five_loops() {
    check(
        "qed_ase_L5",
        SearchCensus {
            topologies: 519,
            wasted_topologies: 0,
            partial_topologies: 12596,
            dead_partial_topologies: 11655,
            useless_partial_topologies: 11655,
            dressings: 2490381,
            dead_dressings: 480051,
            useless_dressings: 896742,
            diagrams: 322704,
            rejected: [0, 0, 0],
        },
    );
}
