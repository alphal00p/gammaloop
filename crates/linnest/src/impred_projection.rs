//! Native ownership and constraint projection for constrained ImPrEd layouts.
//! Route samples remain drawing geometry and never alter physical incidence.
use crate::{TreeInitCfg, TypstGraph};
use cgmath::{InnerSpace, Point2};
use linnet::half_edge::involution::{EdgeIndex, Flow, HedgePair};
use serde_json::{Value, json};
#[cfg(any(test, feature = "projection-service"))]
use std::collections::BTreeSet;
use std::{collections::BTreeMap, error::Error};

pub type Result<T> = std::result::Result<T, Box<dyn Error>>;

enum CoordinateProjection {
    InitialSeed,
    Proposal,
}

pub struct ProjectionSession {
    template: Value,
    initial: BTreeMap<String, Point2<f64>>,
    routes: Value,
    nodes: BTreeMap<usize, String>,
    anchors: BTreeMap<usize, String>,
    #[cfg(any(test, feature = "projection-service"))]
    half_routes: BTreeMap<usize, Vec<String>>,
    fixed_axes: BTreeMap<String, [bool; 2]>,
    records: Vec<Value>,
    external_pull_scales: BTreeMap<String, f64>,
}

fn coordinates(value: &Value) -> Result<BTreeMap<String, Point2<f64>>> {
    value
        .as_object()
        .ok_or("positions must be an object")?
        .iter()
        .map(|(key, value)| Ok((key.clone(), point(value)?)))
        .collect()
}

fn id_keys(value: &Value, name: &str) -> Result<BTreeMap<usize, String>> {
    value
        .as_object()
        .ok_or_else(|| format!("{name} must be an object"))?
        .iter()
        .map(|(key, value)| {
            Ok((
                key.parse()?,
                value.as_str().ok_or("invalid position key")?.to_owned(),
            ))
        })
        .collect()
}

fn position_path(kind: &str, index: usize) -> String {
    if kind == "node" {
        format!("/graph/node_store/node_data/{index}")
    } else {
        format!("/graph/edge_store/data/{index}/0")
    }
}

impl ProjectionSession {
    pub fn new(graph: TypstGraph, seed: Value) -> Result<Self> {
        let template = serde_json::to_value(graph)?;
        let initial = coordinates(&seed["positions"])?;
        let nodes = id_keys(&seed["node_ids"], "node_ids")?;
        let anchors = id_keys(&seed["anchor_ids"], "anchor_ids")?;
        let graph: TypstGraph = serde_json::from_value(template.clone())?;
        if nodes.keys().copied().ne(0..graph.n_nodes())
            || anchors.keys().copied().ne(0..graph.n_edges())
        {
            return Err("seed IDs do not match graph topology".into());
        }
        let routes: BTreeMap<usize, Vec<String>> = serde_json::from_value(seed["routes"].clone())?;
        if routes.keys().copied().ne(0..graph.n_edges()) {
            return Err("route IDs do not match graph topology".into());
        }
        let fixed_nodes = template["layout_config"]["layout-nodes"] == "fixed";
        let mut fixed_axes = BTreeMap::new();
        let mut records = Vec::new();
        for (kind, keys) in [("node", &nodes), ("edge", &anchors)] {
            for (&id, key) in keys {
                let item = template
                    .pointer(&position_path(kind, id))
                    .ok_or("missing graph point")?;
                let p = initial.get(key).ok_or("missing seed point")?;
                let stored: Point2<f64> = serde_json::from_value(item["pos"].clone())?;
                if (stored.x - p.x).abs() > 1e-10 || (stored.y - p.y).abs() > 1e-10 {
                    return Err(format!("seed and projected graph disagree at {key}").into());
                }
                let frozen = ["x", "y"].map(|axis| {
                    (kind == "node" && fixed_nodes) || item["constraints"][axis] == "Fixed"
                });
                fixed_axes.insert(key.clone(), frozen);
                records.push(json!({"kind":kind,"id":id,"key":key,
                    "constraints":item["constraints"],"shift":item["shift"],
                    "fixed_axes":frozen,"initial_position":[p.x,p.y]}));
            }
        }
        let mut half_routes = BTreeMap::new();
        for (&edge, route) in &routes {
            if route.len() < 2 || route.iter().any(|key| !initial.contains_key(key)) {
                return Err(format!("invalid route {edge}").into());
            }
            let anchor = &anchors[&edge];
            let slots: Vec<_> = route
                .iter()
                .enumerate()
                .filter(|(_, key)| *key == anchor)
                .map(|(i, _)| i)
                .collect();
            if slots.len() != 1 {
                return Err(format!("edge {edge} needs one distinguished anchor").into());
            }
            let slot = slots[0];
            let (_, pair) = &graph[&EdgeIndex(edge)];
            let halves = match *pair {
                HedgePair::Paired { source, sink } => {
                    if slot == 0
                        || slot == route.len() - 1
                        || route.first() != nodes.get(&graph.node_id(source).0)
                        || route.last() != nodes.get(&graph.node_id(sink).0)
                    {
                        return Err(format!("paired edge {edge} identity mismatch").into());
                    }
                    vec![
                        (source, route[1..slot].to_vec()),
                        (
                            sink,
                            route[slot + 1..route.len() - 1]
                                .iter()
                                .rev()
                                .cloned()
                                .collect(),
                        ),
                    ]
                }
                HedgePair::Unpaired { hedge, flow } => {
                    let source = flow == Flow::Source;
                    if slot != if source { route.len() - 1 } else { 0 }
                        || (if source { route.first() } else { route.last() })
                            != nodes.get(&graph.node_id(hedge).0)
                    {
                        return Err(format!("dangling edge {edge} identity mismatch").into());
                    }
                    let mut keys = route[1..route.len() - 1].to_vec();
                    if !source {
                        keys.reverse();
                    }
                    vec![(hedge, keys)]
                }
                HedgePair::Split { .. } => {
                    return Err("split edges must be separated before ImPrEd projection".into());
                }
            };
            for (hedge, keys) in halves {
                let stored: Vec<Point2<f64>> = serde_json::from_value(
                    template["graph"]["hedge_data"][hedge.0]["route_points"].clone(),
                )?;
                if keys.len() != stored.len()
                    || keys.iter().zip(&stored).any(|(key, p)| {
                        (initial[key].x - p.x).abs() > 1e-10 || (initial[key].y - p.y).abs() > 1e-10
                    })
                {
                    return Err(format!(
                        "seed route and projected graph disagree at half-edge {}",
                        hedge.0
                    )
                    .into());
                }
                let frozen = fixed_axes[anchor];
                for key in &keys {
                    if fixed_axes.insert(key.clone(), frozen).is_some() {
                        return Err("interior route point has multiple owners".into());
                    }
                    let p = initial[key];
                    records.push(json!({"kind":"route","key":key,"edge":edge,"hedge":hedge.0,
                        "fixed_axes":frozen,"initial_position":[p.x,p.y],
                        "constraints":{"x":if frozen[0]{"Fixed"}else{"Free"},"y":if frozen[1]{"Fixed"}else{"Free"}}}));
                }
                half_routes.insert(hedge.0, keys);
            }
        }
        if fixed_axes.len() != initial.len() {
            return Err("seed contains unowned position keys".into());
        }
        let (state, energy) = graph.layout_energy_state();
        let external_pull_scales = anchors
            .iter()
            .filter(|(edge, _)| matches!(graph[&EdgeIndex(**edge)].1, HedgePair::Unpaired { .. }))
            .map(|(&edge, key)| {
                (
                    key.clone(),
                    energy.external_pull_scale(state.external_pull_topology(EdgeIndex(edge))),
                )
            })
            .collect();
        Ok(Self {
            template,
            initial,
            routes: seed["routes"].clone(),
            nodes,
            anchors,
            #[cfg(any(test, feature = "projection-service"))]
            half_routes,
            fixed_axes,
            records,
            external_pull_scales,
        })
    }

    pub fn ready(&self, positions: &Value) -> Value {
        json!({"ok":true,"ready":true,"protocol":"native-projection-jsonl-v3",
            "method":"native constraint projection without force integration",
            "node_ids":self.nodes,"anchor_ids":self.anchors,"routes":self.routes,
            "positions":positions,"constraints":self.records,
            "external_pull_scales":self.external_pull_scales})
    }

    #[cfg(any(test, feature = "projection-service"))]
    pub fn request(&mut self, request: &Value) -> Result<Value> {
        match request.get("op").and_then(Value::as_str) {
            None if request.get("op").is_none() => self.project(request),
            Some("project") => self.project(request),
            Some("update_external_routes") => self.update_external_routes(request),
            Some("result") => self.result(request),
            _ => Err("unknown projection operation".into()),
        }
    }

    #[cfg(any(test, feature = "projection-service"))]
    fn update_external_routes(&mut self, request: &Value) -> Result<Value> {
        let positions = coordinates(&request["positions"])?;
        let removed: BTreeSet<_> = self
            .initial
            .keys()
            .filter(|key| !positions.contains_key(*key))
            .cloned()
            .collect();
        if removed
            .iter()
            .any(|key| self.fixed_axes[key].iter().any(|fixed| *fixed))
        {
            return Err("external updates must retain every pinned route point".into());
        }
        let old_positions: BTreeMap<_, _> = self
            .initial
            .iter()
            .map(|(key, initial)| {
                let p = positions.get(key).unwrap_or(initial);
                (key.clone(), vec![p.x, p.y])
            })
            .collect();
        let projected = self.project(&json!({"positions":old_positions}))?;
        let projected = coordinates(&projected["positions"])?;
        if projected.iter().any(|(key, p)| {
            positions.get(key).is_some_and(|position| {
                (p.x - position.x).abs() > 1e-10 || (p.y - position.y).abs() > 1e-10
            })
        }) {
            return Err("refinement input must already satisfy native constraints".into());
        }
        let routes: BTreeMap<usize, Vec<String>> =
            serde_json::from_value(request["routes"].clone())?;
        let old_routes: BTreeMap<usize, Vec<String>> = serde_json::from_value(self.routes.clone())?;
        if routes.keys().ne(old_routes.keys()) {
            return Err("refinement must retain physical edge IDs".into());
        }
        let graph: TypstGraph = serde_json::from_value(self.template.clone())?;
        let mut inserted = BTreeSet::new();
        let mut initial: BTreeMap<_, _> = self
            .initial
            .iter()
            .filter(|(key, _)| !removed.contains(*key))
            .map(|(key, p)| (key.clone(), *p))
            .collect();
        let mut contracted = BTreeSet::new();
        let mut template = self.template.clone();
        for (&edge, route) in &routes {
            let old = &old_routes[&edge];
            if route == old {
                continue;
            }
            let HedgePair::Unpaired { hedge, flow } = graph[&EdgeIndex(edge)].1 else {
                return Err("only external routes may be updated".into());
            };
            if route.len() < 2
                || route.len() > 5
                || route.first() != old.first()
                || route.last() != old.last()
                || route.iter().any(|key| !positions.contains_key(key))
            {
                return Err(
                    "external updates must retain endpoints and at most three interior points"
                        .into(),
                );
            }
            let retained_old: Vec<_> = old.iter().filter(|key| !removed.contains(*key)).collect();
            let retained_new: Vec<_> = route
                .iter()
                .filter(|key| self.initial.contains_key(*key))
                .collect();
            if retained_old != retained_new {
                return Err("external updates must retain surviving route keys in order".into());
            }
            contracted.extend(
                old[1..old.len() - 1]
                    .iter()
                    .filter(|key| removed.contains(*key))
                    .cloned(),
            );
            for window in route.windows(3) {
                let key = &window[1];
                if self.initial.contains_key(key) {
                    continue;
                }
                if !old
                    .windows(2)
                    .any(|segment| segment[0] == window[0] && segment[1] == window[2])
                {
                    return Err(
                        "insert one midpoint per existing segment, separately from contraction"
                            .into(),
                    );
                }
                if !inserted.insert(key.clone()) {
                    return Err("inserted midpoint must have a new unique position key".into());
                }
                let p = *positions
                    .get(key)
                    .ok_or("missing inserted midpoint position")?;
                let a = positions[&window[0]];
                let b = positions[&window[2]];
                let scale = (b.x - a.x).hypot(b.y - a.y);
                if scale <= 1e-12
                    || (p.x - (a.x + b.x) * 0.5).hypot(p.y - (a.y + b.y) * 0.5)
                        > 1e-10 * scale.max(1.0)
                {
                    return Err(
                        "new external route point must be a nondegenerate segment midpoint".into(),
                    );
                }
                initial.insert(key.clone(), p);
            }
            let mut keys = route[1..route.len() - 1].to_vec();
            if flow != Flow::Source {
                keys.reverse();
            }
            template["graph"]["hedge_data"][hedge.0]["route_points"] =
                json!(keys.iter().map(|key| initial[key]).collect::<Vec<_>>());
        }
        if (inserted.is_empty() && removed.is_empty())
            || contracted != removed
            || positions.keys().ne(initial.keys())
        {
            return Err("external updates must change only owned free route points".into());
        }
        // Retain the original template and every existing fixed-coordinate
        // reference. A newly inserted point gets its own reference at insertion.
        let candidate = Self::new(
            serde_json::from_value(template)?,
            json!({
                "positions":initial.iter().map(|(key,p)|(key,[p.x,p.y])).collect::<BTreeMap<_,_>>(),
                "node_ids":self.nodes,"anchor_ids":self.anchors,"routes":routes,
            }),
        )?;
        let output = candidate.project(request)?;
        let actual = coordinates(&output["positions"])?;
        if positions.iter().any(|(key, p)| {
            (p.x - actual[key].x).abs() > 1e-10 || (p.y - actual[key].y).abs() > 1e-10
        }) {
            return Err("native update changed the supplied route geometry".into());
        }
        let mut ready = candidate.ready(&output["positions"]);
        ready["updated"] = json!(true);
        ready["inserted_points"] = json!(inserted.len());
        ready["removed_points"] = json!(removed.len());
        *self = candidate;
        Ok(ready)
    }

    #[cfg(any(test, feature = "projection-service"))]
    fn project(&self, request: &Value) -> Result<Value> {
        if request
            .get("routes")
            .is_some_and(|routes| routes != &self.routes)
        {
            return Err(
                "routes are immutable in this session; start a new session after adaptation".into(),
            );
        }
        let mut positions = coordinates(&request["positions"])?;
        if positions.keys().ne(self.initial.keys()) {
            return Err("positions must contain exactly the session keys".into());
        }
        for (key, p) in &mut positions {
            let frozen = self.fixed_axes[key];
            if frozen[0] {
                p.x = self.initial[key].x;
            }
            if frozen[1] {
                p.y = self.initial[key].y;
            }
        }
        let mut candidate = self.template.clone();
        for (kind, keys) in [("node", &self.nodes), ("edge", &self.anchors)] {
            for (&id, key) in keys {
                let p = positions[key];
                candidate
                    .pointer_mut(&position_path(kind, id))
                    .ok_or("missing graph point")?["pos"] = json!({"x":p.x,"y":p.y});
            }
        }
        for (&hedge, keys) in &self.half_routes {
            candidate["graph"]["hedge_data"][hedge]["route_points"] =
                json!(keys.iter().map(|key| positions[key]).collect::<Vec<_>>());
        }
        let mut graph: TypstGraph = serde_json::from_value(candidate)?;
        graph.project_coordinates(CoordinateProjection::Proposal)?;
        let output = serde_json::to_value(&graph)?;
        for (kind, keys) in [("node", &self.nodes), ("edge", &self.anchors)] {
            for (&id, key) in keys {
                let item = output
                    .pointer(&position_path(kind, id))
                    .ok_or("missing result point")?;
                positions.insert(key.clone(), serde_json::from_value(item["pos"].clone())?);
            }
        }
        for (&hedge, keys) in &self.half_routes {
            let values: Vec<Point2<f64>> = serde_json::from_value(
                output["graph"]["hedge_data"][hedge]["route_points"].clone(),
            )?;
            if keys.len() != values.len() {
                return Err("native projection changed route shape".into());
            }
            for (key, p) in keys.iter().zip(values) {
                positions.insert(key.clone(), p);
            }
        }
        for (key, p) in &positions {
            if !p.x.is_finite() || !p.y.is_finite() {
                return Err("native projection returned nonfinite geometry".into());
            }
            let frozen = self.fixed_axes[key];
            if (frozen[0] && p.x != self.initial[key].x)
                || (frozen[1] && p.y != self.initial[key].y)
            {
                return Err(format!("native projection moved frozen coordinate at {key}").into());
            }
        }
        Ok(
            json!({"ok":true,"positions":positions.iter().map(|(key,p)|(key,[p.x,p.y])).collect::<BTreeMap<_,_>>()}),
        )
    }
}

fn point(value: &Value) -> Result<Point2<f64>> {
    let xy = value.as_array().ok_or("expected coordinate pair")?;
    if xy.len() != 2 {
        return Err("expected two coordinates".into());
    }
    let x = xy[0].as_f64().ok_or("invalid x")?;
    let y = xy[1].as_f64().ok_or("invalid y")?;
    if !x.is_finite() || !y.is_finite() {
        return Err("non-finite coordinate".into());
    }
    Ok(Point2::new(x, y))
}

impl TypstGraph {
    /// Export physical incidence for the embedding-constrained initializer.
    /// Legs sharing a grouped y coordinate share a `row`, as the two opened
    /// halves of a cross-section initial state do.
    pub fn impred_diagram(&self) -> Result<Value> {
        let mut edges = Vec::new();
        for i in 0..self.n_edges() {
            let row = match self.graph[EdgeIndex(i)].constraints.y {
                crate::Constraint::Grouped(crate::LayoutPointIndex::Edge(edge), _) => {
                    Some(format!("edge:{}", edge.0))
                }
                crate::Constraint::Grouped(crate::LayoutPointIndex::Node(node), _) => {
                    Some(format!("node:{}", node.0))
                }
                _ => None,
            };
            let (source, target, state) = match self[&EdgeIndex(i)].1 {
                HedgePair::Paired { source, sink } => (
                    Some(self.node_id(source).0),
                    Some(self.node_id(sink).0),
                    None,
                ),
                HedgePair::Unpaired {
                    hedge,
                    flow: Flow::Source,
                } => (Some(self.node_id(hedge).0), None, Some("outgoing")),
                HedgePair::Unpaired {
                    hedge,
                    flow: Flow::Sink,
                } => (None, Some(self.node_id(hedge).0), Some("incoming")),
                HedgePair::Split { .. } => {
                    return Err("split edges must be separated before ImPrEd layout".into());
                }
            };
            let mut edge = json!({"id":i,"source":source,"target":target,"state":state});
            if let (Some(row), Some(_)) = (row, state) {
                edge["row"] = json!(row);
            }
            edges.push(edge);
        }
        Ok(json!({"vertices":(0..self.n_nodes()).collect::<Vec<_>>(), "edges":edges}))
    }

    fn project_coordinates(&mut self, stage: CoordinateProjection) -> Result<()> {
        self.validate_layout()?;
        // Group references are user input: reject invalid targets before indexing.
        for constraint in (0..self.n_nodes())
            .map(|i| self[crate::NodeIndex(i)].constraints)
            .chain((0..self.n_edges()).map(|i| self[EdgeIndex(i)].constraints))
        {
            for axis in [constraint.x, constraint.y] {
                if let crate::Constraint::Grouped(reference, _) = axis {
                    let valid = match reference {
                        crate::LayoutPointIndex::Node(node) => node.0 < self.n_nodes(),
                        crate::LayoutPointIndex::Edge(edge) => edge.0 < self.n_edges(),
                    };
                    if !valid {
                        return Err("unknown native grouped-coordinate reference".into());
                    }
                }
            }
        }
        self.layout_initialized = true;
        let (mut nodes, mut edges) =
            self.initial_solver_positions(TreeInitCfg { dx: 1.0, dy: 1.0 });
        // The EC seed already supplies separated coordinates. Initialize a shared
        // axis from its canonical reference: averaging opposite-side seed rows can
        // collapse distinct groups (for example, reversed cross-section cut pairs).
        // Later proposals retain the native start/averaging and sign projection.
        if matches!(stage, CoordinateProjection::Proposal) {
            self.apply_initial_grouped_constraints(&mut nodes, &mut edges);
        }
        self.apply_group_references(
            &mut nodes,
            &mut edges,
            self.layout_config.layout_nodes.nodes_are_fixed(),
        );
        self.apply_layout_constraints(&mut nodes, &mut edges);
        let routes = self.initial_solver_route_points();
        self.update_positions(nodes, edges, routes, None, None);
        Ok(())
    }
}

impl ProjectionSession {
    /// Attach a certified EC carrier to the original styled graph, retaining its
    /// pins, coordinate groups, shifts, identities and user drawing attributes.
    pub fn initialize(mut graph: TypstGraph, mut seed: Value) -> Result<(Self, Value)> {
        let mut positions = coordinates(&seed["positions"])?;
        let nodes = id_keys(&seed["node_ids"], "node_ids")?;
        let mut routes: BTreeMap<usize, Vec<String>> =
            serde_json::from_value(seed["routes"].clone())?;
        if nodes.keys().copied().ne(0..graph.n_nodes())
            || routes.keys().copied().ne(0..graph.n_edges())
        {
            return Err("seed IDs do not match physical graph topology".into());
        }
        let before_crossings = carrier_crossings(&positions, &routes)?;
        let mut anchors = BTreeMap::new();
        let mut halves = BTreeMap::new();
        for (&edge, route) in &mut routes {
            if route.len() < 2 || route.iter().any(|key| !positions.contains_key(key)) {
                return Err(format!("invalid route {edge}").into());
            }
            let pair = graph[&EdgeIndex(edge)].1;
            let slot = match pair {
                HedgePair::Paired { source, sink } => {
                    if route.first() != nodes.get(&graph.node_id(source).0)
                        || route.last() != nodes.get(&graph.node_id(sink).0)
                    {
                        return Err(
                            format!("paired edge {edge} endpoint identities disagree").into()
                        );
                    }
                    if route.len() == 2 {
                        let key = format!("native:anchor:{edge}");
                        if positions.contains_key(&key) {
                            return Err("reserved native anchor identity is in use".into());
                        }
                        positions.insert(
                            key.clone(),
                            positions[&route[0]]
                                + (positions[&route[1]] - positions[&route[0]]) * 0.5,
                        );
                        route.insert(1, key);
                        1
                    } else {
                        let spans: Vec<_> = route
                            .windows(2)
                            .map(|p| (positions[&p[1]] - positions[&p[0]]).magnitude())
                            .collect();
                        let length: f64 = spans.iter().sum();
                        let mut distance = 0.0;
                        let mut nearest = (f64::INFINITY, 1);
                        for (i, span) in spans.iter().enumerate().take(spans.len() - 1) {
                            distance += span;
                            let error = (distance - length * 0.5).abs();
                            if error < nearest.0 {
                                nearest = (error, i + 1);
                            }
                        }
                        nearest.1
                    }
                }
                HedgePair::Unpaired { hedge, flow } => {
                    let source = flow == Flow::Source;
                    if (if source { route.first() } else { route.last() })
                        != nodes.get(&graph.node_id(hedge).0)
                    {
                        return Err(
                            format!("external edge {edge} endpoint identity disagrees").into()
                        );
                    }
                    if source { route.len() - 1 } else { 0 }
                }
                HedgePair::Split { .. } => {
                    return Err("split edges must be separated before ImPrEd layout".into());
                }
            };
            anchors.insert(edge, route[slot].clone());
            match pair {
                HedgePair::Paired { source, sink } => {
                    halves.insert(source.0, route[1..slot].to_vec());
                    halves.insert(
                        sink.0,
                        route[slot + 1..route.len() - 1]
                            .iter()
                            .rev()
                            .cloned()
                            .collect(),
                    );
                }
                HedgePair::Unpaired { hedge, flow } => {
                    let mut keys = route[1..route.len() - 1].to_vec();
                    if flow == Flow::Sink {
                        keys.reverse();
                    }
                    halves.insert(hedge.0, keys);
                }
                HedgePair::Split { .. } => unreachable!(),
            }
        }
        let fixed_nodes = graph.layout_config.layout_nodes.nodes_are_fixed();
        let mut authored_pins = Vec::new();
        for (&id, key) in &nodes {
            let data = &mut graph.graph[crate::NodeIndex(id)];
            let mut p = positions[key];
            let fixed = [
                fixed_nodes || matches!(data.constraints.x, crate::Constraint::Fixed),
                fixed_nodes || matches!(data.constraints.y, crate::Constraint::Fixed),
            ];
            if fixed[0] {
                p.x = data.pos.x;
            }
            if fixed[1] {
                p.y = data.pos.y;
            }
            authored_pins.push((key.clone(), fixed, p));
            positions.insert(key.clone(), p);
            data.pos = p;
        }
        for (&id, key) in &anchors {
            let data = &mut graph.graph[EdgeIndex(id)];
            let mut p = positions[key];
            let fixed = [
                matches!(data.constraints.x, crate::Constraint::Fixed),
                matches!(data.constraints.y, crate::Constraint::Fixed),
            ];
            if fixed[0] {
                p.x = data.pos.x;
            }
            if fixed[1] {
                p.y = data.pos.y;
            }
            authored_pins.push((key.clone(), fixed, p));
            positions.insert(key.clone(), p);
            data.pos = p;
        }
        for (&hedge, keys) in &halves {
            graph.graph[crate::Hedge(hedge)].route_points =
                keys.iter().map(|key| positions[key]).collect();
        }
        graph.project_coordinates(CoordinateProjection::InitialSeed)?;
        for (&id, key) in &nodes {
            positions.insert(key.clone(), graph[crate::NodeIndex(id)].pos);
        }
        for (&id, key) in &anchors {
            positions.insert(key.clone(), graph[EdgeIndex(id)].pos);
        }
        for (&hedge, keys) in &halves {
            for (key, &p) in keys.iter().zip(&graph[crate::Hedge(hedge)].route_points) {
                positions.insert(key.clone(), p);
            }
        }
        for (key, axes, original) in authored_pins {
            let actual = positions[&key];
            if (axes[0] && actual.x != original.x) || (axes[1] && actual.y != original.y) {
                return Err(format!("native seed projection moved authored pin at {key}").into());
            }
        }
        // Anchoring can subdivide a segment, so compare physical edge pairs,
        // including multiplicity, rather than changing segment indices.
        if carrier_crossings(&positions, &routes)? != before_crossings {
            return Err("native seed constraints changed the embedding's crossings".into());
        }
        seed["positions"] = json!(
            positions
                .iter()
                .map(|(key, p)| (key, [p.x, p.y]))
                .collect::<BTreeMap<_, _>>()
        );
        seed["routes"] = json!(routes);
        seed["anchor_ids"] = json!(anchors);
        if let Some(object) = seed.as_object_mut() {
            object.remove("route_cubics");
        }
        let projected_graph = serde_json::to_value(&graph)?;
        let session = Self::new(graph, seed.clone())?;
        let ready = session.ready(&seed["positions"]);
        Ok((
            session,
            json!({"ok":true,"state":seed,"graph":projected_graph,"ready":ready}),
        ))
    }

    /// Return dense display-coordinate arrays for the renderer's shared import.
    #[cfg(any(test, feature = "projection-service"))]
    pub fn result(&self, request: &Value) -> Result<Value> {
        let projected = self.project(request)?;
        let positions = coordinates(&projected["positions"])?;
        let mut graph: TypstGraph = serde_json::from_value(self.template.clone())?;
        for (&id, key) in &self.nodes {
            graph.graph[crate::NodeIndex(id)].pos = positions[key];
        }
        for (&id, key) in &self.anchors {
            graph.graph[EdgeIndex(id)].pos = positions[key];
        }
        for (&id, keys) in &self.half_routes {
            graph.graph[crate::Hedge(id)].route_points =
                keys.iter().map(|key| positions[key]).collect();
        }
        graph.project_coordinates(CoordinateProjection::Proposal)?;
        let pair = |p: Point2<f64>| [p.x, p.y];
        let geometry = json!({
            "nodes":(0..graph.n_nodes()).map(|i| pair(graph[crate::NodeIndex(i)].pos)).collect::<Vec<_>>(),
            "edges":(0..graph.n_edges()).map(|i| pair(graph[EdgeIndex(i)].pos)).collect::<Vec<_>>(),
            "routes":(0..graph.n_hedges()).map(|i| graph[crate::Hedge(i)].route_points.iter().copied().map(pair).collect::<Vec<_>>()).collect::<Vec<_>>(),
        });
        Ok(json!({"ok":true,"geometry":geometry,"graph":graph}))
    }
}

/// Physical crossing witnesses for the seed handoff. Proper crossings are
/// permitted; contacts, coincident spans and degenerate segments are rejected.
fn carrier_crossings(
    positions: &BTreeMap<String, Point2<f64>>,
    routes: &BTreeMap<usize, Vec<String>>,
) -> Result<Vec<(usize, usize)>> {
    let mut segments = Vec::new();
    for (&edge, route) in routes {
        if route.len() < 2 {
            return Err("a carrier needs two endpoints".into());
        }
        for pair in route.windows(2) {
            let a = *positions.get(&pair[0]).ok_or("missing carrier point")?;
            let b = *positions.get(&pair[1]).ok_or("missing carrier point")?;
            if (b - a).magnitude() <= 1e-8 {
                return Err(format!("degenerate carrier segment on edge {edge}").into());
            }
            segments.push((edge, &pair[0], &pair[1], a, b));
        }
    }
    let cross = |a: Point2<f64>, b: Point2<f64>, c: Point2<f64>| {
        (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
    };
    let on = |a: Point2<f64>, b: Point2<f64>, p: Point2<f64>| {
        cross(a, b, p).abs() <= 1e-8
            && p.x >= a.x.min(b.x) - 1e-8
            && p.x <= a.x.max(b.x) + 1e-8
            && p.y >= a.y.min(b.y) - 1e-8
            && p.y <= a.y.max(b.y) + 1e-8
    };
    let mut crossings = Vec::new();
    for (i, &(edge, aid, bid, a, b)) in segments.iter().enumerate() {
        for &(other, cid, did, c, d) in &segments[i + 1..] {
            let v = [
                cross(a, b, c),
                cross(a, b, d),
                cross(c, d, a),
                cross(c, d, b),
            ];
            if v.iter().all(|value| value.abs() <= 1e-8) {
                let interval = |p: Point2<f64>| {
                    if (a.x - b.x).abs() >= (a.y - b.y).abs() {
                        p.x
                    } else {
                        p.y
                    }
                };
                let extent = interval(a)
                    .max(interval(b))
                    .min(interval(c).max(interval(d)))
                    - interval(a)
                        .min(interval(b))
                        .max(interval(c).min(interval(d)));
                if extent > 1e-8 {
                    return Err(format!("overlapping carrier edges {edge} and {other}").into());
                }
            }
            let s = v.map(|value| {
                if value.abs() <= 1e-8 {
                    0.0
                } else {
                    value.signum()
                }
            });
            if s[0] * s[1] < 0.0 && s[2] * s[3] < 0.0 {
                crossings.push((edge.min(other), edge.max(other)));
            } else if aid != cid
                && aid != did
                && bid != cid
                && bid != did
                && (on(a, b, c) || on(a, b, d) || on(c, d, a) || on(c, d, b))
            {
                return Err(format!("contact between carrier edges {edge} and {other}").into());
            }
        }
    }
    crossings.sort_unstable();
    Ok(crossings)
}

/// In-process JSON protocol shared by native layout hosts and the command tool.
#[cfg(any(test, feature = "projection-service"))]
#[derive(Default)]
pub struct ProjectionService {
    session: Option<ProjectionSession>,
}

#[cfg(any(test, feature = "projection-service"))]
impl ProjectionService {
    pub fn handle(&mut self, request: Value) -> std::result::Result<Value, String> {
        self.dispatch(request).map_err(|error| error.to_string())
    }

    fn dispatch(&mut self, request: Value) -> Result<Value> {
        match request["op"].as_str() {
            Some("describe") => {
                let graph = TypstGraph::from_layout_snapshot(&request["graph"])?;
                Ok(json!({"ok":true,"diagram":graph.impred_diagram()?}))
            }
            Some("initialize") => {
                let graph = TypstGraph::from_layout_snapshot(&request["graph"])?;
                let (session, response) =
                    ProjectionSession::initialize(graph, request["state"].clone())?;
                self.session = Some(session);
                Ok(response)
            }
            Some("resume") => {
                let graph = TypstGraph::from_layout_snapshot(&request["graph"])?;
                let state = request["state"].clone();
                let session = ProjectionSession::new(graph, state.clone())?;
                let response = session.ready(&state["positions"]);
                self.session = Some(session);
                Ok(response)
            }
            _ => self
                .session
                .as_mut()
                .ok_or("initialize a projection session first")?
                .request(&request),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn edge_seed() -> Value {
        json!({
            "positions":{"v:0":[0.0,0.0],"v:1":[8.0,0.0]},
            "node_ids":{"0":"v:0","1":"v:1"},
            "routes":{"0":["v:0","v:1"]},
            "external_ids":{},"incoming_ids":[],"outgoing_ids":[],"report":{},
        })
    }

    #[test]
    fn initializer_preserves_authored_pins_and_original_payload() {
        let graph = TypstGraph::parse(
            r#"digraph {
            a [id=0 pos="x:2!" label="keep node"]
            b [id=1 pos="y:-1!"]
            a -> b [id=0 pos="x:5!" label="keep edge" particle="g"]
        }"#,
        )
        .unwrap();
        let (mut session, initialized) =
            ProjectionSession::initialize(graph.clone(), edge_seed()).unwrap();
        assert_eq!(initialized["state"]["positions"]["v:0"], json!([2.0, 0.0]));
        assert_eq!(initialized["state"]["positions"]["v:1"], json!([8.0, -1.0]));
        assert_eq!(
            initialized["state"]["positions"]["native:anchor:0"],
            json!([5.0, 0.0])
        );
        let final_state = session
            .request(&json!({"op":"result","positions":{
                "v:0":[-100.0,1.0],"v:1":[9.0,100.0],"native:anchor:0":[100.0,1.0]
            }}))
            .unwrap();
        assert_eq!(
            final_state["geometry"]["nodes"],
            json!([[2.0, 1.0], [9.0, -1.0]])
        );
        assert_eq!(final_state["geometry"]["edges"], json!([[5.0, 1.0]]));
        let result: TypstGraph = serde_json::from_value(final_state["graph"].clone()).unwrap();
        assert_eq!(
            result[EdgeIndex(0)].statements,
            graph[EdgeIndex(0)].statements
        );
        assert_eq!(
            result[crate::NodeIndex(0)].statements,
            graph[crate::NodeIndex(0)].statements
        );
    }

    #[test]
    fn external_refinement_retains_identity_and_exports_half_routes() {
        let graph = TypstGraph::parse(
            r#"digraph {
            ext [style=invis]
            a [id=0]
            ext -> a:0 [id=0]
            a:1 -> ext [id=1]
        }"#,
        )
        .unwrap();
        let diagram = graph.impred_diagram().unwrap();
        assert_eq!(diagram["edges"][0]["state"], "incoming");
        assert_eq!(diagram["edges"][1]["state"], "outgoing");
        let seed = json!({
            "positions":{"v:0":[0.0,0.0],"x:0":[-3.0,0.0],"x:1":[3.0,0.0]},
            "node_ids":{"0":"v:0"},"routes":{"0":["x:0","v:0"],"1":["v:0","x:1"]},
            "external_ids":{"0":"x:0","1":"x:1"},"incoming_ids":["x:0"],"outgoing_ids":["x:1"],"report":{}
        });
        let (mut session, reply) = ProjectionSession::initialize(graph, seed).unwrap();
        let mut state = reply["state"].clone();
        state["positions"]["b:new"] = json!([1.5, 0.0]);
        state["routes"]["1"] = json!(["v:0", "b:new", "x:1"]);
        let updated = session.request(&json!({"op":"update_external_routes","positions":state["positions"],"routes":state["routes"]})).unwrap();
        assert_eq!(updated["inserted_points"], 1);
        let result = session
            .result(&json!({"positions":state["positions"]}))
            .unwrap();
        assert_eq!(result["geometry"]["routes"], json!([[], [[1.5, 0.0]]]));
        assert_eq!(
            result["geometry"]["edges"],
            json!([[-3.0, 0.0], [3.0, 0.0]])
        );
        state["positions"].as_object_mut().unwrap().remove("b:new");
        state["routes"]["1"] = json!(["v:0", "x:1"]);
        let contracted = session.request(&json!({"op":"update_external_routes","positions":state["positions"],"routes":state["routes"]})).unwrap();
        assert_eq!(contracted["removed_points"], 1);
        assert_eq!(contracted["constraints"], reply["ready"]["constraints"]);
    }

    #[test]
    fn initializer_keeps_reversed_cut_groups_separate() {
        let graph = TypstGraph::parse(
            r#"digraph {
                left [id=0]
                right [id=1]
                ext [style=invis]
                right -> ext [id=0 pos="3,-1" pin="x:@+right,y:@cut-10" "group-start-y"=true]
                right -> ext [id=1 pos="3,1" pin="x:@+right,y:@cut-20" "group-start-y"=true]
                ext -> left [id=2 pos="-3,-1" pin="x:@-left,y:@cut-20" "group-start-y"=true]
                ext -> left [id=3 pos="-3,1" pin="x:@-left,y:@cut-10" "group-start-y"=true]
                left -> right [id=4]
            }"#,
        )
        .unwrap();
        let seed = json!({
            "positions":{"v:0":[-1.0,0.0],"v:1":[1.0,0.0],
                "x:0":[3.0,-1.0],"x:1":[3.0,1.0],
                "x:2":[-3.0,-1.0],"x:3":[-3.0,1.0]},
            "node_ids":{"0":"v:0","1":"v:1"},
            "routes":{"0":["v:1","x:0"],"1":["v:1","x:1"],
                "2":["x:2","v:0"],"3":["x:3","v:0"],"4":["v:0","v:1"]},
            "external_ids":{"0":"x:0","1":"x:1","2":"x:2","3":"x:3"},
            "incoming_ids":["x:2","x:3"],"outgoing_ids":["x:0","x:1"],"report":{}
        });
        let (session, initialized) = ProjectionSession::initialize(graph, seed).unwrap();
        let mut positions = initialized["state"]["positions"].clone();
        assert_eq!(positions["x:0"], json!([3.0, -1.0]));
        assert_eq!(positions["x:1"], json!([3.0, 1.0]));
        assert_eq!(positions["x:2"], json!([-3.0, 1.0]));
        assert_eq!(positions["x:3"], json!([-3.0, -1.0]));
        let mut layout = session.impred_layout().unwrap();
        layout
            .solve(linnet::half_edge::layout::impred::ImpredConfig {
                steps: 0,
                ..Default::default()
            })
            .unwrap();

        // Subsequent proposals still average group-start members, rather than
        // granting their canonical representative a different movement weight.
        for (edge, y) in [(0, -2.0), (1, 2.0), (2, 6.0), (3, -6.0)] {
            positions[format!("x:{edge}")][1] = json!(y);
        }
        let projected = session.project(&json!({"positions":positions})).unwrap();
        assert_eq!(projected["positions"]["x:0"][1], -4.0);
        assert_eq!(projected["positions"]["x:3"][1], -4.0);
        assert_eq!(projected["positions"]["x:1"][1], 4.0);
        assert_eq!(projected["positions"]["x:2"][1], 4.0);
    }

    #[test]
    fn initializer_rejects_pins_that_destroy_the_carrier() {
        let graph = TypstGraph::parse(
            r#"digraph {
            a [id=0 pos="x:8!"]
            b [id=1]
            a -> b [id=0]
        }"#,
        )
        .unwrap();
        let error = ProjectionSession::initialize(graph, edge_seed())
            .err()
            .unwrap();
        assert!(error.to_string().contains("overlapping carrier"));
    }

    #[test]
    fn native_solver_import_preserves_pins_and_physical_payload() {
        let graph = TypstGraph::parse(
            r#"digraph {
                a [id=0 pos="x:2!" label="keep node"]
                b [id=1 pos="y:-1!"]
                a:0 -> b:1 [id=0 pos="x:5!" label="keep edge"]
            }"#,
        )
        .unwrap();
        let (session, _) = ProjectionSession::initialize(graph.clone(), edge_seed()).unwrap();
        let mut layout = session.impred_layout().unwrap();
        let imported = session.apply_impred_layout(&layout).unwrap();
        assert_eq!(imported[crate::NodeIndex(0)].pos, Point2::new(2.0, 0.0));
        assert_eq!(imported[crate::NodeIndex(1)].pos, Point2::new(8.0, -1.0));
        assert_eq!(imported[EdgeIndex(0)].pos, Point2::new(5.0, 0.0));
        // Authored statements survive; the paired label gains its carrier gap.
        let mut statements = graph[EdgeIndex(0)].statements.clone();
        statements.insert("layout-label-gap".into(), "0".into());
        assert_eq!(imported[EdgeIndex(0)].statements, statements);
        assert_eq!(imported.n_nodes(), graph.n_nodes());
        assert_eq!(imported.n_edges(), graph.n_edges());
        layout.positions[layout.node_points[0]][0] += 1.0;
        assert!(
            session
                .apply_impred_layout(&layout)
                .err()
                .unwrap()
                .to_string()
                .contains("authored fixed coordinate")
        );
    }

    #[test]
    fn native_solver_import_keeps_source_and_sink_route_orientation() {
        let graph =
            TypstGraph::parse(r#"digraph { a [id=0] b [id=1] a:0 -> b:1 [id=0] }"#).unwrap();
        let seed = json!({
            "positions":{"v:0":[0.0,0.0],"v:1":[8.0,0.0],
                "b:0":[1.0,1.0],"b:1":[4.0,2.0],"b:2":[7.0,1.0]},
            "node_ids":{"0":"v:0","1":"v:1"},
            "routes":{"0":["v:0","b:0","b:1","b:2","v:1"]},
            "external_ids":{},"incoming_ids":[],"outgoing_ids":[],"report":{}
        });
        let (session, _) = ProjectionSession::initialize(graph, seed).unwrap();
        let layout = session.impred_layout().unwrap();
        let imported = session.apply_impred_layout(&layout).unwrap();
        assert_eq!(imported[EdgeIndex(0)].pos, Point2::new(4.0, 2.0));
        assert_eq!(
            imported.graph[linnet::half_edge::involution::Hedge(0)].route_points,
            vec![Point2::new(1.0, 1.0)]
        );
        assert_eq!(
            imported.graph[linnet::half_edge::involution::Hedge(1)].route_points,
            vec![Point2::new(7.0, 1.0)]
        );
    }

    #[test]
    fn impred_labels_follow_paired_carriers_in_every_label_layout() {
        let mut graph = TypstGraph::parse(
            r#"digraph {
                ext [style=invis]
                a [id=0]
                b [id=1]
                ext -> a:0 [id=0 "layout-label-gap"="0.7"]
                a:1 -> b:2 [id=1 "layout-label-gap"="0.7"]
            }"#,
        )
        .unwrap();
        let seed = json!({
            "positions":{"v:0":[0.0,0.0],"v:1":[4.0,0.0],"x:0":[-3.0,0.0]},
            "node_ids":{"0":"v:0","1":"v:1"},
            "routes":{"0":["x:0","v:0"],"1":["v:0","v:1"]},
            "external_ids":{"0":"x:0"},"incoming_ids":["x:0"],"outgoing_ids":[],"report":{}
        });
        for label_layout in [
            crate::LabelLayout::Normal,
            crate::LabelLayout::DanglingTangent,
            crate::LabelLayout::FixedLength,
            crate::LabelLayout::FixedGap,
        ] {
            graph.layout_config.label_layout = label_layout;
            let (session, _) = ProjectionSession::initialize(graph.clone(), seed.clone()).unwrap();
            let imported = session
                .apply_impred_layout(&session.impred_layout().unwrap())
                .unwrap();
            let gaps: Vec<_> = (0..2)
                .map(|edge| {
                    imported[EdgeIndex(edge)]
                        .statements
                        .get("layout-label-gap")
                        .cloned()
                })
                .collect();
            // Dangling gaps from an earlier force layout never survive ImPrEd.
            assert_eq!(gaps, [None, Some("0".to_owned())]);
        }
    }
}

impl ProjectionSession {
    /// Build the shared native solver's dense coordinate view from the same
    /// immutable references used by projection and adaptive route ownership.
    pub fn impred_layout(&self) -> Result<linnet::half_edge::layout::impred::ImpredLayout> {
        use linnet::half_edge::layout::impred::{
            EdgeLabel, ExternalFlow, ImpredLayout,
            movement::{AxisConstraint, PointRecord},
        };
        let keys: Vec<_> = self.initial.keys().cloned().collect();
        let ids: BTreeMap<_, _> = keys
            .iter()
            .enumerate()
            .map(|(id, key)| (key.clone(), id))
            .collect();
        let native_ids: BTreeMap<_, _> = self
            .records
            .iter()
            .filter_map(|record| {
                record["id"].as_u64().map(|id| {
                    (
                        (record["kind"].as_str().unwrap(), id as usize),
                        record["key"].as_str().unwrap(),
                    )
                })
            })
            .collect();
        let records: BTreeMap<_, _> = self
            .records
            .iter()
            .map(|record| (record["key"].as_str().unwrap(), record))
            .collect();
        let graph: TypstGraph = serde_json::from_value(self.template.clone())?;
        let mut constraints = Vec::with_capacity(keys.len());
        for key in &keys {
            let record = records[key.as_str()];
            let mut axes = [AxisConstraint::Free, AxisConstraint::Free];
            for (axis, name) in ["x", "y"].into_iter().enumerate() {
                let constraint: crate::Constraint =
                    serde_json::from_value(record["constraints"][name].clone())?;
                axes[axis] = match constraint {
                    crate::Constraint::Free => AxisConstraint::Free,
                    crate::Constraint::Fixed => AxisConstraint::Fixed,
                    crate::Constraint::Grouped(reference, direction) => {
                        let identity = match reference {
                            crate::LayoutPointIndex::Node(node) => ("node", node.0),
                            crate::LayoutPointIndex::Edge(edge) => ("edge", edge.0),
                        };
                        let key = native_ids
                            .get(&identity)
                            .ok_or("unknown native grouped-coordinate reference")?;
                        AxisConstraint::Grouped {
                            point: ids[*key],
                            direction,
                        }
                    }
                };
            }
            let shift = if record["shift"].is_null() {
                [0.0, 0.0]
            } else {
                let value: cgmath::Vector2<f64> = serde_json::from_value(record["shift"].clone())?;
                [value.x, value.y]
            };
            let reference = point(&record["initial_position"])?;
            constraints.push(PointRecord {
                constraints: axes,
                shift,
                reference: [reference.x, reference.y],
                fixed_axes: serde_json::from_value(record["fixed_axes"].clone())?,
                is_node: record["kind"] == "node",
            });
        }
        let routes: BTreeMap<usize, Vec<String>> = serde_json::from_value(self.routes.clone())?;
        let external: Vec<_> = (0..graph.n_edges())
            .map(|edge| match graph[&EdgeIndex(edge)].1 {
                HedgePair::Unpaired {
                    flow: Flow::Source, ..
                } => Some(ExternalFlow::Outgoing),
                HedgePair::Unpaired {
                    flow: Flow::Sink, ..
                } => Some(ExternalFlow::Incoming),
                _ => None,
            })
            .collect();
        let pull_scales = (0..graph.n_edges())
            .map(|edge| {
                if external[edge].is_some() {
                    self.external_pull_scales[&self.anchors[&edge]]
                } else {
                    1.0
                }
            })
            .collect();
        // `graph.style` measures each label box before layout, in graph units.
        let labels = (0..graph.n_edges())
            .map(|edge| {
                let statements = &graph[EdgeIndex(edge)].statements;
                let size = |key: &str| {
                    statements
                        .get(key)
                        .and_then(|value| value.trim().parse::<f64>().ok())
                        .filter(|value| value.is_finite() && *value > 0.0)
                };
                match (size("label-width"), size("label-height")) {
                    (Some(width), Some(height)) if external[edge].is_none() => Some(EdgeLabel {
                        extents: [width / 2.0, height / 2.0],
                        center: None,
                        pinned: false,
                        pin: None,
                    }),
                    _ => None,
                }
            })
            .collect();
        Ok(ImpredLayout {
            positions: keys
                .iter()
                .map(|key| {
                    let p = self.initial[key];
                    [p.x, p.y]
                })
                .collect(),
            routes: routes
                .values()
                .map(|route| route.iter().map(|key| ids[key]).collect())
                .collect(),
            node_points: self.nodes.values().map(|key| ids[key]).collect(),
            anchor_points: self.anchors.values().map(|key| ids[key]).collect(),
            external,
            constraints,
            pull_scales,
            samples: keys.iter().map(|key| key.starts_with("b:")).collect(),
            nodes_fixed: graph.layout_config.layout_nodes.nodes_are_fixed(),
            labels,
        })
    }

    /// Import solved routes into the original graph without interpreting their
    /// samples as interaction vertices or rebuilding the physical topology.
    pub fn apply_impred_layout(
        &self,
        layout: &linnet::half_edge::layout::impred::ImpredLayout,
    ) -> Result<TypstGraph> {
        let mut graph: TypstGraph = serde_json::from_value(self.template.clone())?;
        if layout.node_points.len() != graph.n_nodes()
            || layout.anchor_points.len() != graph.n_edges()
            || layout.routes.len() != graph.n_edges()
        {
            return Err("ImPrEd result changed physical node or edge identity".into());
        }
        let point = |index: usize| -> Result<Point2<f64>> {
            let p = layout
                .positions
                .get(index)
                .ok_or("invalid solved point identity")?;
            if !p[0].is_finite() || !p[1].is_finite() {
                return Err("nonfinite solved ImPrEd geometry".into());
            }
            Ok(Point2::new(p[0], p[1]))
        };
        let validate_pin = |key: &str, position: Point2<f64>| -> Result<()> {
            let reference = self.initial[key];
            for axis in 0..2 {
                if self.fixed_axes[key][axis] && (position[axis] - reference[axis]).abs() > 1e-9 {
                    return Err("ImPrEd result moved an authored fixed coordinate".into());
                }
            }
            Ok(())
        };
        for (node, &index) in layout.node_points.iter().enumerate() {
            let position = point(index)?;
            validate_pin(&self.nodes[&node], position)?;
            graph.graph[crate::NodeIndex(node)].pos = position;
        }
        for (edge, &index) in layout.anchor_points.iter().enumerate() {
            let position = point(index)?;
            validate_pin(&self.anchors[&edge], position)?;
            graph.graph[EdgeIndex(edge)].pos = position;
        }
        for (edge, route) in layout.routes.iter().enumerate() {
            let slots: Vec<_> = route
                .iter()
                .enumerate()
                .filter(|(_, index)| **index == layout.anchor_points[edge])
                .map(|(slot, _)| slot)
                .collect();
            if slots.len() != 1 || route.len() < 2 {
                return Err("solved route lost its physical anchor".into());
            }
            let slot = slots[0];
            match graph[&EdgeIndex(edge)].1 {
                HedgePair::Paired { source, sink } => {
                    if slot == 0
                        || slot + 1 == route.len()
                        || route[0] != layout.node_points[graph.node_id(source).0]
                        || *route.last().unwrap() != layout.node_points[graph.node_id(sink).0]
                    {
                        return Err("solved paired route changed physical incidence".into());
                    }
                    graph.graph[source].route_points = route[1..slot]
                        .iter()
                        .map(|&index| point(index))
                        .collect::<Result<_>>()?;
                    graph.graph[sink].route_points = route[slot + 1..route.len() - 1]
                        .iter()
                        .rev()
                        .map(|&index| point(index))
                        .collect::<Result<_>>()?;
                }
                HedgePair::Unpaired { hedge, flow } => {
                    let source = flow == Flow::Source;
                    if slot != if source { route.len() - 1 } else { 0 }
                        || (if source {
                            route[0]
                        } else {
                            *route.last().unwrap()
                        }) != layout.node_points[graph.node_id(hedge).0]
                    {
                        return Err("solved dangling route changed physical incidence".into());
                    }
                    let mut route: Vec<_> = route[1..route.len() - 1]
                        .iter()
                        .map(|&index| point(index))
                        .collect::<Result<_>>()?;
                    if !source {
                        route.reverse();
                    }
                    graph.graph[hedge].route_points = route;
                }
                HedgePair::Split { .. } => {
                    return Err("split edges must be separated before ImPrEd layout".into());
                }
            }
        }
        // In every label layout, paired labels follow their carriers: drawing
        // widens this zero gap to clear the painted line, and the annotation
        // search places the labels. ImPrEd's label point, beside the carrier,
        // becomes the label position whose side the search prefers.
        for edge in 0..graph.n_edges() {
            let paired = matches!(graph[&EdgeIndex(edge)].1, HedgePair::Paired { .. });
            let center = layout
                .labels
                .get(edge)
                .copied()
                .flatten()
                .and_then(|label| label.center);
            let data = &mut graph.graph[EdgeIndex(edge)];
            if paired {
                data.statements
                    .insert("layout-label-gap".into(), "0".into());
            } else {
                data.statements.remove("layout-label-gap");
            }
            if let Some([x, y]) = center {
                data.label_pos = Some(Point2::new(x, y));
            }
        }
        graph.project_coordinates(CoordinateProjection::Proposal)?;
        Ok(graph)
    }
}
