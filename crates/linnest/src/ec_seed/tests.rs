use super::embedding::{ConstraintKind, ConstraintTree, Expansion, Graph, OrderedMap};
use super::spqr::{
    Decomposition, SplitDecomposer, SpqrComponent, SpqrDecomposer, SpqrKind, SpqrRequest,
};
use super::*;
use indexmap::IndexMap;
use serde_json::json;
use std::collections::HashMap;

/// Reference seeds and constrained embeddings of the native C++ engine.
const SEEDS: &str = include_str!("../../layout/native/fixtures/seeds.json");
const CONSTRAINTS: &str = include_str!("../../layout/native/fixtures/constraints.json");
/// SPQR decompositions by the native engine's OGDF primitive (`ec-spqr`
/// `decompose` responses) of every block the fixtures decompose, keyed by
/// request.
const NATIVE_SPQR: &str = include_str!("fixtures/native_spqr.json");
/// Seeds of the fixture diagrams with the split decomposer; see
/// `regenerate_split_seeds`.
const SPLIT_SEEDS: &str = include_str!("fixtures/split_seeds.json");

/// Replays the native engine's SPQR decompositions.
struct NativeSpqr(HashMap<String, Value>);

impl NativeSpqr {
    fn load() -> Self {
        Self(serde_json::from_str(NATIVE_SPQR).unwrap())
    }

    fn decomposition(response: &Value) -> Result<Option<Decomposition>> {
        if !response["planar"].as_bool().unwrap() {
            return Ok(None);
        }
        require(
            response["biconnected"].as_bool().unwrap(),
            "SPQR requires a biconnected block",
        )?;
        serde_json::from_value(response.clone())
            .map(Some)
            .map_err(|error| error.to_string())
    }
}

impl SpqrDecomposer for NativeSpqr {
    fn decompose(&self, block: &SpqrRequest<'_>) -> Result<Option<Decomposition>> {
        let request = serde_json::to_string(block).unwrap();
        let response = self
            .0
            .get(&request)
            .unwrap_or_else(|| panic!("no native decomposition of {request}"));
        Self::decomposition(response)
    }
}

#[derive(Deserialize)]
struct SeedFixture {
    name: String,
    diagram: Value,
    expected: Value,
}

fn seed_fixtures() -> Vec<SeedFixture> {
    serde_json::from_str(SEEDS).unwrap()
}

fn seed(diagram: &Value, spqr: &dyn SpqrDecomposer) -> Value {
    let request = SeedRequest {
        diagram: Diagram::try_from(diagram.clone()).unwrap(),
        scale: SeedRequest::default_scale(),
        external_sides: SeedRequest::default_external_sides(),
    };
    serde_json::to_value(request.initialize_with(spqr).unwrap()).unwrap()
}

/// The checks of the native `test_layout`, with bit-identical coordinates.
fn check_seed_fixtures(spqr: &dyn SpqrDecomposer) {
    for fixture in seed_fixtures() {
        let (name, expected) = (&fixture.name, &fixture.expected);
        let result = seed(&fixture.diagram, spqr);
        for key in ["routes", "node_ids", "external_ids", "edge_endpoints"] {
            assert_eq!(result[key], expected[key], "{name}: {key}");
        }
        let positions = result["positions"].as_object().unwrap();
        let expected_positions = expected["positions"].as_object().unwrap();
        assert_eq!(positions.len(), expected_positions.len(), "{name}");
        for (point, coordinate) in expected_positions {
            for axis in 0..2 {
                assert_eq!(
                    positions[point][axis].as_f64().unwrap().to_bits(),
                    coordinate[axis].as_f64().unwrap().to_bits(),
                    "{name}: {point}"
                );
            }
        }
        let crossings = result["crossings"].as_array().unwrap();
        let expected_crossings = expected["crossings"].as_array().unwrap();
        assert_eq!(crossings.len(), expected_crossings.len(), "{name}");
        for (actual, expected) in crossings.iter().zip(expected_crossings) {
            for key in ["id", "edges", "route_points"] {
                assert_eq!(actual[key], expected[key], "{name}: crossing {key}");
            }
        }
        for key in ["contacts", "overlaps", "degenerate"] {
            assert_eq!(
                result["report"]["geometry"][key],
                json!([]),
                "{name}: {key}"
            );
        }
    }
}

#[test]
fn native_spqr_reproduces_reference_seeds() {
    check_seed_fixtures(&NativeSpqr::load());
}

#[derive(Deserialize)]
struct ConstraintFixture {
    request: ConstraintRequest,
    expected: ConstraintExpectation,
}

#[derive(Deserialize)]
struct ConstraintRequest {
    nodes: Set,
    edges: Vec<FixtureEdge>,
    rotation: BTreeMap<Id, Ids>,
    /// Constraints apply in file order.
    constraints: IndexMap<Id, Value>,
}

#[derive(Deserialize)]
struct FixtureEdge {
    id: Id,
    source: Id,
    target: Id,
}

#[derive(Deserialize)]
struct ConstraintExpectation {
    /// Feasibility by exhaustive rotation enumeration.
    planar: bool,
    #[serde(default)]
    rotation: BTreeMap<Id, Ids>,
    #[serde(default)]
    edges: BTreeMap<Id, Ends>,
}

fn constraint_tree(value: &Value) -> ConstraintTree {
    if let Some(leaf) = value.as_str() {
        return ConstraintTree::Leaf(leaf.to_owned());
    }
    let kind = match value["kind"].as_str().unwrap() {
        "group" => ConstraintKind::Group,
        "mirror" => ConstraintKind::Mirror,
        "oriented" => ConstraintKind::Oriented,
        other => panic!("unknown constraint kind {other}"),
    };
    let children = value["children"]
        .as_array()
        .unwrap()
        .iter()
        .map(constraint_tree)
        .collect();
    ConstraintTree::Node { kind, children }
}

/// Constrained embedding of every fixture collapsed to its original vertices,
/// after checking feasibility and physical incidence.
fn constrained_embeddings(spqr: &dyn SpqrDecomposer) -> Vec<(ConstraintFixture, Option<Graph>)> {
    let fixtures: Vec<ConstraintFixture> = serde_json::from_str(CONSTRAINTS).unwrap();
    fixtures
        .into_iter()
        .map(|fixture| {
            let request = &fixture.request;
            let edges = request
                .edges
                .iter()
                .map(|e| (e.id.clone(), [e.source.clone(), e.target.clone()]))
                .collect();
            let mut graph = Graph::new(request.nodes.clone(), edges).unwrap();
            graph.rotation = request.rotation.clone();
            let constraints: OrderedMap<_> = request
                .constraints
                .iter()
                .map(|(n, t)| (n.clone(), constraint_tree(t)))
                .collect();
            let expansion = Expansion::new(&graph, constraints).unwrap();
            let embedded = expansion.graph.try_embed(spqr).unwrap();
            assert_eq!(embedded.is_some(), fixture.expected.planar);
            // Collapsing validates the embedding and the constraints.
            let collapsed = embedded.map(|embedded| expansion.collapse(&embedded).unwrap());
            for (edge, ends) in collapsed.iter().flat_map(|g| g.edges.iter()) {
                assert_eq!(
                    *ends, fixture.expected.edges[edge],
                    "physical incidence of {edge}"
                );
            }
            (fixture, collapsed)
        })
        .collect()
}

#[test]
fn native_spqr_reproduces_reference_constrained_embeddings() {
    for (fixture, collapsed) in constrained_embeddings(&NativeSpqr::load()) {
        if let Some(collapsed) = collapsed {
            assert_eq!(collapsed.rotation, fixture.expected.rotation);
        }
    }
}

#[test]
fn split_decomposer_realizes_feasible_constrained_embeddings() {
    for (fixture, collapsed) in constrained_embeddings(&SplitDecomposer) {
        if let Some(collapsed) = collapsed {
            assert_eq!(collapsed.nodes, fixture.request.nodes);
        }
    }
}

type ComponentShape = (SpqrKind, Ids, Vec<Ids>, BTreeSet<Vec<Ids>>);

/// Decomposition up to identities: every skeleton edge is named by the real
/// edges it represents, and R skeleton rotations are compared up to a common
/// reflection.
fn shape(decomposition: &Decomposition) -> BTreeSet<ComponentShape> {
    let components: BTreeMap<&str, &SpqrComponent> = decomposition
        .components
        .iter()
        .map(|c| (c.id.as_str(), c))
        .collect();
    fn behind(components: &BTreeMap<&str, &SpqrComponent>, component: &str, edge: &str) -> Ids {
        let e = components[component]
            .edges
            .iter()
            .find(|e| e.id == edge)
            .unwrap();
        if let Some(real) = &e.real_edge {
            return vec![real.clone()];
        }
        let twin = e.twin.as_ref().unwrap();
        let mut out: Ids = components[twin.component.as_str()]
            .edges
            .iter()
            .filter(|f| f.id != twin.edge)
            .flat_map(|f| behind(components, &twin.component, &f.id))
            .collect();
        out.sort();
        out
    }
    decomposition
        .components
        .iter()
        .map(|c| {
            let label = |e: &str| behind(&components, &c.id, e);
            let mut edges: Vec<Ids> = c.edges.iter().map(|e| label(&e.id)).collect();
            edges.sort();
            let rotations = |reflect: bool| -> BTreeSet<Vec<Ids>> {
                c.rotation_cw
                    .values()
                    .map(|r| {
                        let mut r: Vec<Ids> = r.iter().map(|e| label(e)).collect();
                        if reflect {
                            r.reverse();
                        }
                        let start = (0..r.len()).min_by_key(|&i| &r[i]).unwrap();
                        r.rotate_left(start);
                        r
                    })
                    .collect()
            };
            let rotation = match c.kind {
                SpqrKind::R => rotations(false).min(rotations(true)),
                SpqrKind::S | SpqrKind::P => BTreeSet::new(),
            };
            let mut nodes = c.nodes.clone();
            nodes.sort();
            (c.kind, nodes, edges, rotation)
        })
        .collect()
}

#[test]
fn split_decomposer_matches_native_spqr_trees() {
    for (request, response) in NativeSpqr::load().0 {
        let request: Value = serde_json::from_str(&request).unwrap();
        let nodes: Set = serde_json::from_value(request["nodes"].clone()).unwrap();
        let edges = request["edges"]
            .as_array()
            .unwrap()
            .iter()
            .map(|e| {
                (
                    json_key(&e["id"]),
                    [json_key(&e["source"]), json_key(&e["target"])],
                )
            })
            .collect();
        let block = Graph::new(nodes, edges).unwrap();
        let native = NativeSpqr::decomposition(&response).unwrap();
        let split = SplitDecomposer
            .decompose(&SpqrRequest::new(&block))
            .unwrap();
        assert_eq!(native.is_some(), split.is_some(), "planarity of {request}");
        if let (Some(native), Some(split)) = (native, split) {
            assert_eq!(shape(&native), shape(&split), "SPQR tree of {request}");
        }
    }
}

#[test]
fn split_decomposer_seeds_are_stable() {
    let expected: BTreeMap<String, Value> = serde_json::from_str(SPLIT_SEEDS).unwrap();
    for fixture in seed_fixtures() {
        let result = seed(&fixture.diagram, &SplitDecomposer);
        for key in ["node_ids", "external_ids", "edge_endpoints"] {
            assert_eq!(
                result[key], fixture.expected[key],
                "{}: {key}",
                fixture.name
            );
        }
        let geometry = &result["report"]["geometry"];
        for key in ["contacts", "overlaps", "degenerate"] {
            assert_eq!(geometry[key], json!([]), "{}: {key}", fixture.name);
        }
        assert_eq!(
            result, expected[&fixture.name],
            "{}: run regenerate_split_seeds after intended changes",
            fixture.name
        );
    }
}

#[test]
#[ignore = "rewrites fixtures/split_seeds.json"]
fn regenerate_split_seeds() {
    let seeds: BTreeMap<String, Value> = seed_fixtures()
        .into_iter()
        .map(|f| (f.name, seed(&f.diagram, &SplitDecomposer)))
        .collect();
    let entries: Vec<String> = seeds
        .iter()
        .map(|(name, seed)| format!("{}:{seed}", json!(name)))
        .collect();
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/src/ec_seed/fixtures/split_seeds.json"
    );
    std::fs::write(path, format!("{{\n{}\n}}\n", entries.join(",\n"))).unwrap();
}

#[test]
fn seed_protocol_reads_edge_pairs_and_defaults() {
    let request: SeedRequest = serde_json::from_value(json!({
        "diagram": {
            "vertices": [0, 1],
            "edges": [
                [{"source": 0, "target": 1}, {}],
                [{"source": 0, "target": 1}, {}],
                [{"source": null, "target": 0}, {"external": {"state": "incoming"}}],
                [{"source": 1, "target": null}, {"external": {"state": "outgoing"}}],
            ],
        },
    }))
    .unwrap();
    assert_eq!((request.scale, request.external_sides), (2.4, true));
    let seed = serde_json::to_value(request.initialize().unwrap()).unwrap();
    assert_eq!(seed["external_ids"], json!({"2": "x:2", "3": "x:3"}));
    assert_eq!(seed["incoming_ids"], json!(["x:2"]));
    assert_eq!(seed["routes"]["0"], json!(["v:0", "b:0:0", "v:1"]));
    assert_eq!(seed["report"]["status"], "ec-planar");
}

#[test]
fn paired_rows_share_one_order_on_both_sides() {
    // A forward box whose opened initial states pair opposite corners: rows
    // can align only by crossing two internal lines.
    let edge = |id: usize, source: Value, target: Value, state: Value, row: Value| json!({"id": id, "source": source, "target": target, "state": state, "row": row});
    let request: SeedRequest = serde_json::from_value(json!({
        "diagram": {
            "vertices": [0, 1, 2, 3],
            "edges": [
                edge(0, json!(1), Value::Null, json!("outgoing"), json!("edge:0")),
                edge(1, json!(0), Value::Null, json!("outgoing"), json!("edge:1")),
                edge(2, json!(1), json!(0), Value::Null, Value::Null),
                edge(3, json!(0), json!(3), Value::Null, Value::Null),
                edge(4, json!(1), json!(2), Value::Null, Value::Null),
                edge(5, json!(3), json!(2), Value::Null, Value::Null),
                edge(6, Value::Null, json!(3), json!("incoming"), json!("edge:0")),
                edge(7, Value::Null, json!(2), json!("incoming"), json!("edge:1")),
            ],
        },
    }))
    .unwrap();
    let seed = serde_json::to_value(request.initialize().unwrap()).unwrap();
    assert_eq!(seed["report"]["status"], "ec-planarized");
    assert_eq!(seed["crossings"].as_array().unwrap().len(), 1);
    // Partners keep the same vertical order on both sides.
    let height = |id: &str| seed["positions"][format!("x:{id}")][1].as_f64().unwrap();
    assert_eq!(height("0") > height("1"), height("6") > height("7"));
    // Without rows the box stays planar, with its legs in opposite orders.
    let mut unpaired = request.clone();
    for edge in &mut unpaired.diagram.edges {
        edge.row = None;
    }
    let seed = serde_json::to_value(unpaired.initialize().unwrap()).unwrap();
    assert_eq!(seed["report"]["status"], "ec-planar");
}

#[test]
fn rows_sharing_a_fan_follow_their_partners() {
    // A crossed Compton forward graph: both outgoing legs fan from one vertex,
    // so only the fan's spread can keep each row level.
    let edge = |id: usize, source: Value, target: Value, state: Value, row: Value| json!({"id": id, "source": source, "target": target, "state": state, "row": row});
    let request: SeedRequest = serde_json::from_value(json!({
        "diagram": {
            "vertices": [0, 1, 2, 3],
            "edges": [
                edge(0, json!(0), Value::Null, json!("outgoing"), json!("edge:0")),
                edge(1, json!(0), Value::Null, json!("outgoing"), json!("edge:1")),
                edge(2, json!(2), json!(0), Value::Null, Value::Null),
                edge(3, json!(1), json!(2), Value::Null, Value::Null),
                edge(4, json!(3), json!(1), Value::Null, Value::Null),
                edge(5, json!(2), json!(3), Value::Null, Value::Null),
                edge(6, Value::Null, json!(3), json!("incoming"), json!("edge:0")),
                edge(7, Value::Null, json!(1), json!("incoming"), json!("edge:1")),
            ],
        },
    }))
    .unwrap();
    let seed = serde_json::to_value(request.initialize().unwrap()).unwrap();
    assert_eq!(seed["report"]["status"], "ec-planar");
    let height = |id: &str| seed["positions"][format!("x:{id}")][1].as_f64().unwrap();
    assert_eq!(height("0") > height("1"), height("6") > height("7"));
}

#[test]
fn seed_protocol_rejects_invalid_requests() {
    let error = |request: Value| {
        serde_json::from_value::<SeedRequest>(request)
            .map_err(|e| e.to_string())
            .and_then(|request| request.initialize().map(|_| ()))
            .unwrap_err()
    };
    let diagram = json!({"vertices": [0], "edges": [{"source": 0, "target": 1}]});
    assert_eq!(
        error(json!({ "diagram": diagram })),
        "Unknown edge endpoint"
    );
    let diagram = json!({"vertices": [0], "edges": [{"source": null, "target": null}]});
    assert_eq!(error(json!({ "diagram": diagram })), "Edge has no endpoint");
    let diagram = json!({"vertices": [0], "edges": []});
    assert_eq!(
        error(json!({"diagram": diagram, "scale": 0.0})),
        "Spacing must be finite and positive"
    );
}

#[test]
fn coordinate_mean_rounds_once() {
    let mean = |values: &[f64]| Magnitude::mean(values.iter().copied(), values.len());
    // Summing in f64 would round 0.1 + 0.2 to 0.30000000000000004 first.
    assert_eq!(mean(&[0.1, 0.2, 0.3]), 0.2);
    assert_eq!(mean(&[1e16, 1.0, -1e16]), 1.0 / 3.0);
    assert_eq!(mean(&[-1.5, -2.5]), -2.0);
    assert_eq!(mean(&[f64::MIN_POSITIVE, -f64::MIN_POSITIVE]).to_bits(), 0);
    assert_eq!(mean(&[]), 0.0);
}
