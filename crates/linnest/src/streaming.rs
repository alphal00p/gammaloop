use cgmath::{Vector2, Zero};
use figment::{providers::Serialized, Figment, Profile};
use linnet::{
    half_edge::{
        involution::{Flow, HedgePair},
        layout::force::ForceLayoutSession,
        nodestore::DefaultNodeStore,
        NodeIndex,
    },
    parser::set::DotGraphSet,
};

use crate::{LayoutAlgo, TypstEdge, TypstGraph, TypstHedge, TypstNode, TypstPoint};

type NativeForceSession<'a> =
    ForceLayoutSession<'a, TypstEdge, TypstNode, TypstHedge, DefaultNodeStore<TypstNode>>;

self_cell::self_cell! {
    struct ForceLayoutCell {
        owner: TypstGraph,
        #[covariant]
        dependent: NativeForceSession,
    }
}

/// A native, synchronous stream of force-layout positions.
///
/// The graph and solver state share a lifetime so each batch continues the
/// original cooling and auxiliary-depth schedules. Labels are not measured or
/// relaxed by this geometry preview.
pub struct ForceLayoutStream {
    cell: ForceLayoutCell,
}

#[derive(Debug, Clone)]
pub struct LayoutFrame {
    pub nodes: Vec<TypstPoint>,
    pub edges: Vec<TypstPoint>,
    /// Complete source-to-sink carriers, including vertices, bends and the edge anchor.
    pub paths: Vec<Vec<TypstPoint>>,
    pub iteration: usize,
    pub done: bool,
    pub max_movement: f64,
}

impl ForceLayoutStream {
    pub fn from_dot(dot: &str, settings: &[u8]) -> Result<Self, String> {
        let dots = DotGraphSet::from_string(dot).map_err(|error| error.to_string())?;
        if dots.set.len() != 1 {
            return Err("a layout stream needs exactly one DOT graph".to_owned());
        }
        let dot = dots.into_iter().next().expect("checked one graph");
        let value: ciborium::Value = ciborium::de::from_reader(settings)
            .map_err(|error| format!("invalid layout settings: {error}"))?;
        let figment = Figment::from(Serialized::from(value, Profile::Default));
        let config = figment
            .extract()
            .map_err(|error| format!("invalid layout settings: {error}"))?;
        let owner = TypstGraph::from_dot_with_layout_config(dot, config);
        Self::from_graph(owner)
    }

    /// Continue the existing force solver from graph geometry, including route bends.
    pub fn from_graph(owner: TypstGraph) -> Result<Self, String> {
        owner.validate_layout()?;
        if !matches!(owner.layout_config.layout_algo, LayoutAlgo::Force) {
            return Err("streaming currently supports only the force layout".to_owned());
        }
        if owner.layout_config.layout_nodes.nodes_are_fixed() {
            return Err(
                "streaming does not support layout-nodes=fixed; use explicit DOT pins".to_owned(),
            );
        }
        let cell = ForceLayoutCell::new(owner, |owner| {
            let (state, energy) = owner.layout_energy_state();
            owner.force_session(state, &energy, None)
        });
        Ok(Self { cell })
    }

    pub fn step(&mut self, iterations: usize) -> LayoutFrame {
        self.cell.with_dependent_mut(|owner, session| {
            let snapshot = session.step(iterations);
            let node_point = |index: NodeIndex| {
                let point = snapshot.vertex_points[index]
                    + owner[index].shift.unwrap_or_else(Vector2::zero);
                TypstPoint {
                    x: point.x,
                    y: point.y,
                }
            };
            let paths = owner
                .graph
                .iter_edges()
                .map(|(pair, edge, _)| {
                    let shift = owner[edge].shift.unwrap_or_else(Vector2::zero);
                    let interior = |hedge| {
                        snapshot.route_points[hedge]
                            .iter()
                            .map(|point| {
                                let point = *point + shift;
                                TypstPoint {
                                    x: point.x,
                                    y: point.y,
                                }
                            })
                            .collect::<Vec<_>>()
                    };
                    let anchor = snapshot.edge_points[edge] + shift;
                    let anchor = TypstPoint {
                        x: anchor.x,
                        y: anchor.y,
                    };
                    let (source, sink) = match pair {
                        HedgePair::Paired { source, sink }
                        | HedgePair::Split { source, sink, .. } => (Some(source), Some(sink)),
                        HedgePair::Unpaired {
                            hedge,
                            flow: Flow::Source,
                        } => (Some(hedge), None),
                        HedgePair::Unpaired {
                            hedge,
                            flow: Flow::Sink,
                        } => (None, Some(hedge)),
                    };
                    let mut path = Vec::new();
                    if let Some(hedge) = source {
                        path.push(node_point(owner.node_id(hedge)));
                        path.extend(interior(hedge));
                    }
                    path.push(anchor);
                    if let Some(hedge) = sink {
                        path.extend(interior(hedge).into_iter().rev());
                        path.push(node_point(owner.node_id(hedge)));
                    }
                    path
                })
                .collect();
            LayoutFrame {
                paths,
                nodes: snapshot
                    .vertex_points
                    .into_iter()
                    .map(|(i, point)| {
                        let point = point + owner[i].shift.unwrap_or_else(Vector2::zero);
                        TypstPoint {
                            x: point.x,
                            y: point.y,
                        }
                    })
                    .collect(),
                edges: snapshot
                    .edge_points
                    .into_iter()
                    .map(|(i, point)| {
                        let point = point + owner[i].shift.unwrap_or_else(Vector2::zero);
                        TypstPoint {
                            x: point.x,
                            y: point.y,
                        }
                    })
                    .collect(),
                iteration: snapshot.iteration,
                done: snapshot.done,
                max_movement: snapshot.max_move,
            }
        })
    }

    pub fn node_names(&self) -> Vec<String> {
        let owner = self.cell.borrow_owner();
        (0..owner.n_nodes())
            .map(|i| {
                owner[NodeIndex(i)]
                    .name
                    .clone()
                    .unwrap_or_else(|| format!("n{i}"))
            })
            .collect()
    }

    pub fn endpoints(&self) -> Vec<(Option<usize>, Option<usize>)> {
        let owner = self.cell.borrow_owner();
        owner
            .graph
            .iter_edges()
            .map(|(pair, _, _)| match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    (Some(owner.node_id(source).0), Some(owner.node_id(sink).0))
                }
                HedgePair::Unpaired {
                    hedge,
                    flow: Flow::Source,
                } => (Some(owner.node_id(hedge).0), None),
                HedgePair::Unpaired {
                    hedge,
                    flow: Flow::Sink,
                } => (None, Some(owner.node_id(hedge).0)),
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use linnet::half_edge::involution::EdgeIndex;

    use super::*;

    #[test]
    fn streamed_batches_match_full_layout_and_preserve_constraints() {
        let dot = r#"digraph {
            a [id=0 pos="1,2!"]
            b [id=1 pos="x:@column!,y:0"]
            c [id=2 pos="x:@column!,y:3"]
            ext [style=invis]
            a -> b [id=0 "spring-length"="1.75"]
            b -> c [id=1]
            c -> ext [id=2]
        }"#;
        let settings = BTreeMap::from([
            ("layout-algo", "force"),
            ("steps", "9"),
            ("epochs", "3"),
            ("seed", "7"),
            ("early-tol", "0"),
        ]);
        let mut bytes = Vec::new();
        ciborium::ser::into_writer(&settings, &mut bytes).unwrap();
        let figment = Figment::from(Serialized::from(settings, Profile::Default));
        let parsed = DotGraphSet::from_string(dot).unwrap();
        let mut graph = TypstGraph::from_dot(parsed.into_iter().next().unwrap(), &figment);
        graph.layout();

        for every in [1, 7, 64] {
            let mut stream = ForceLayoutStream::from_dot(dot, &bytes).unwrap();
            assert_eq!(stream.node_names(), ["a", "b", "c"]);
            assert_eq!(
                stream.endpoints(),
                [(Some(0), Some(1)), (Some(1), Some(2)), (Some(2), None)]
            );
            let initial = stream.step(0);
            assert_eq!(initial.iteration, 0);
            assert!(!initial.done);
            let mut frame = initial.clone();
            while !frame.done {
                assert_eq!(frame.nodes[0], TypstPoint { x: 1.0, y: 2.0 });
                assert_eq!(frame.nodes[1].x, frame.nodes[2].x);
                frame = stream.step(every);
            }
            assert_eq!(frame.iteration, 27);
            for (i, point) in frame.nodes.iter().enumerate() {
                let expected = graph[NodeIndex(i)].pos;
                assert_eq!((point.x, point.y), (expected.x, expected.y));
            }
            for (i, point) in frame.edges.iter().enumerate() {
                let expected = graph[EdgeIndex(i)].pos;
                assert_eq!((point.x, point.y), (expected.x, expected.y));
            }
            assert_eq!(stream.step(every).nodes, frame.nodes);
            assert_eq!(initial.iteration, 0);
        }
    }

    #[test]
    fn routed_stream_paths_match_full_layout_in_both_directions() {
        use cgmath::Point2;
        use linnet::half_edge::involution::Hedge;

        let mut graph = TypstGraph::parse(
            r#"digraph {
            ext [style=invis]
            a [id=0 pos="0,0!"]
            b [id=1 pos="4,0!"]
            ext -> a:0 [id=0 pos="-2,0"]
            a:1 -> b:2 [id=1 pos="2,0"]
            b:3 -> ext [id=2 pos="6,0"]
        }"#,
        )
        .unwrap();
        for (hedge, points) in [
            vec![Point2::new(-1.0, 1.0)],
            vec![Point2::new(0.5, 1.0), Point2::new(1.5, 1.0)],
            vec![Point2::new(3.0, -1.0)],
            vec![Point2::new(5.0, 1.0)],
        ]
        .into_iter()
        .enumerate()
        {
            graph.graph[Hedge(hedge)].route_points = points;
        }
        let settings = BTreeMap::from([
            ("layout-algo", "force"),
            ("steps", "10"),
            ("epochs", "2"),
            ("depth-scale", "0"),
            ("early-tol", "0"),
            ("label-steps", "0"),
        ]);
        graph.layout_config = crate::LayoutConfig::from_figment(&Figment::from(Serialized::from(
            settings,
            Profile::Default,
        )));
        let mut expected = graph.clone();
        expected.layout();
        for batch in [1, 7, 100] {
            let mut stream = ForceLayoutStream::from_graph(graph.clone()).unwrap();
            let initial = stream.step(0);
            assert_eq!(
                initial.paths.iter().map(Vec::len).collect::<Vec<_>>(),
                [3, 6, 3]
            );
            assert_eq!(initial.paths[0][0], initial.edges[0]);
            assert_eq!(initial.paths[0][2], initial.nodes[0]);
            assert_eq!(initial.paths[1][0], initial.nodes[0]);
            assert_eq!(initial.paths[1][5], initial.nodes[1]);
            assert_eq!(initial.paths[2][0], initial.nodes[1]);
            assert_eq!(initial.paths[2][2], initial.edges[2]);
            let mut frame = initial.clone();
            while !frame.done {
                frame = stream.step(batch);
            }
            for (hedge, edge, path_index) in [(0, 0, 1), (1, 1, 1), (2, 1, 4), (3, 2, 1)] {
                let point = expected.graph[Hedge(hedge)].route_points[0];
                assert_eq!(
                    frame.paths[edge][path_index],
                    TypstPoint {
                        x: point.x,
                        y: point.y
                    }
                );
            }
            assert_ne!(frame.paths[1], initial.paths[1]);
        }
    }

    #[test]
    fn edge_spring_length_changes_the_streamed_solution() {
        let settings = BTreeMap::from([
            ("layout-algo", "force"),
            ("steps", "30"),
            ("epochs", "1"),
            ("depth-scale", "0"),
            ("early-tol", "0"),
        ]);
        let mut bytes = Vec::new();
        ciborium::ser::into_writer(&settings, &mut bytes).unwrap();
        let mut solutions = Vec::new();
        for scale in [1, 3] {
            let dot =
                format!("digraph {{ a -> b [\"spring-length\"=\"{scale}\"]; b -> c; c -> a }}");
            let mut stream = ForceLayoutStream::from_dot(&dot, &bytes).unwrap();
            solutions.push(stream.step(30).nodes);
        }
        assert_ne!(solutions[0], solutions[1]);
    }

    #[test]
    fn empty_graphs_and_zero_budgets_finish_with_finite_snapshots() {
        for dot in ["digraph {}", "digraph { a -> b }"] {
            for (steps, epochs) in [(0, 2), (2, 0), (2, 2)] {
                let settings = BTreeMap::from([
                    ("layout-algo", "force".to_owned()),
                    ("steps", steps.to_string()),
                    ("epochs", epochs.to_string()),
                ]);
                let mut bytes = Vec::new();
                ciborium::ser::into_writer(&settings, &mut bytes).unwrap();
                let mut stream = ForceLayoutStream::from_dot(dot, &bytes).unwrap();
                let initial = stream.step(0);
                let frame = stream.step(4);
                assert!(frame.done);
                assert!(frame.max_movement.is_finite());
                assert!(frame
                    .nodes
                    .iter()
                    .chain(&frame.edges)
                    .all(|point| { point.x.is_finite() && point.y.is_finite() }));
                if steps * epochs == 0 {
                    assert!(initial.done);
                    assert_eq!(frame.iteration, 0);
                    assert_eq!(initial.nodes, frame.nodes);
                    assert_eq!(initial.edges, frame.edges);
                }
            }
        }
    }

    #[test]
    fn invalid_dot_and_settings_return_errors() {
        let settings = BTreeMap::from([("layout-algo", "force")]);
        let mut bytes = Vec::new();
        ciborium::ser::into_writer(&settings, &mut bytes).unwrap();
        for dot in [
            "this is not DOT",
            "digraph { a -> }",
            "digraph {} digraph {}",
            "digraph { a -> b [\"spring-length\"=\"0\"] }",
            "digraph { a [\"pos-z\"=\"NaN\"] }",
        ] {
            assert!(ForceLayoutStream::from_dot(dot, &bytes).is_err(), "{dot}");
        }
        assert!(ForceLayoutStream::from_dot("digraph {}", &[]).is_err());
        for (key, value) in [
            ("layout-algo", "tree"),
            ("layout-nodes", "fixed"),
            ("depth-scale", "-1"),
            ("flattening-end", "2"),
            ("steps", "invalid"),
        ] {
            let mut settings = BTreeMap::from([("layout-algo", "force")]);
            settings.insert(key, value);
            let mut bytes = Vec::new();
            ciborium::ser::into_writer(&settings, &mut bytes).unwrap();
            assert!(ForceLayoutStream::from_dot("digraph { a -> b }", &bytes).is_err());
        }
    }
}
