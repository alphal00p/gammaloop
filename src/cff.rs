use std::{
    fmt::{Debug, Display},
    hash::{Hash, Hasher},
};

use crate::{
    graph::{EdgeType, Graph},
    utils::FloatLike,
};
use ahash::{HashMap, HashMapExt, HashSet};
use color_eyre::Report;
use eyre::{eyre, Result};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use symbolica::representations::Atom;

use log::debug;

const MAX_VERTEX_COUNT: usize = 32;

#[derive(Serialize, Deserialize, Debug, Clone, Eq)]
pub struct Esurface {
    energies: Vec<usize>,
    sub_orientation: Vec<bool>,
    shift: Vec<usize>,
    shift_signature: bool,
}

impl PartialEq for Esurface {
    fn eq(&self, other: &Self) -> bool {
        self.energies == other.energies && self.sub_orientation == other.sub_orientation
    }
}

#[allow(unused)]
impl Esurface {
    fn to_atom(&self) -> Atom {
        let symbolic_energies = self
            .energies
            .iter()
            .map(|i| Atom::parse(&format!("E{}", i)).unwrap())
            .collect_vec();

        let symbolic_shift = self
            .shift
            .iter()
            .map(|i| Atom::parse(&format!("p{}", i)).unwrap())
            .collect_vec();

        let builder_atom = Atom::new();
        let energy_sum = symbolic_energies
            .iter()
            .fold(builder_atom, |acc, energy| acc + energy);

        let esurf = symbolic_shift.iter().fold(energy_sum, |acc, shift| {
            if self.shift_signature {
                acc + shift
            } else {
                -(-acc + shift)
            }
        });

        esurf
    }

    // the energy cache contains the energies of external edges as well as the virtual,
    // use the location in the supergraph to determine the index
    fn compute_value<T: FloatLike>(&self, energy_cache: &[T]) -> T {
        let energy_sum = self
            .energies
            .iter()
            .map(|index| energy_cache[*index])
            .sum::<T>();

        let shift_sum = self
            .shift
            .iter()
            .map(|index| energy_cache[*index])
            .sum::<T>();

        let shift_sign = match self.shift_signature {
            true => Into::<T>::into(1.),
            false => Into::<T>::into(-1.),
        };

        energy_sum + shift_sign * shift_sum
    }
}

#[derive(Debug, Clone)]
pub struct CFFTree {
    nodes: Vec<CFFTreeNode>,
    orientation: Orientation,
    term_id: usize,
    num_data_nodes: usize,
}

impl CFFTree {
    fn from_root_graph(
        graph: CFFIntermediateGraph,
        orientation: Orientation,
        term_id: usize,
    ) -> Self {
        let node = CFFTreeNodeData {
            node_id: 0,
            graph,
            children: vec![],
            _parent: None,
            esurface_id: None,
        };

        let tree = vec![CFFTreeNode::Data(node)];
        Self {
            nodes: tree,
            orientation,
            term_id,
            num_data_nodes: 1,
        }
    }

    fn insert_graph(&mut self, parent: usize, graph: CFFIntermediateGraph) -> Result<(), Report> {
        let node_id = self.nodes.len();

        let node = CFFTreeNodeData {
            node_id,
            graph,
            children: vec![],
            _parent: Some(parent),
            esurface_id: None,
        };

        self.nodes.push(CFFTreeNode::Data(node));
        match self.nodes[parent] {
            CFFTreeNode::Data(ref mut tree_node_data) => {
                tree_node_data.children.push(node_id);
                self.num_data_nodes += 1;
                Ok(())
            }
            CFFTreeNode::Pointer(_) => Err(eyre!("Parent node is a pointer")),
        }
    }

    fn insert_pointer(&mut self, parent: usize, pointer: CFFTreeNodePointer) -> Result<(), Report> {
        let node_id = self.nodes.len();
        self.nodes.push(CFFTreeNode::Pointer(pointer));

        match self.nodes[parent] {
            CFFTreeNode::Data(ref mut tree_node_data) => {
                tree_node_data.children.push(node_id);
                Ok(())
            }
            CFFTreeNode::Pointer(_) => Err(eyre!("Parent node is a pointer")),
        }
    }

    fn insert_esurface(&mut self, node_id: usize, esurface_id: usize) -> Result<(), Report> {
        match self.nodes[node_id] {
            CFFTreeNode::Data(ref mut tree_node_data) => {
                tree_node_data.esurface_id = Some(esurface_id);
                Ok(())
            }
            CFFTreeNode::Pointer(_) => Err(eyre!("Node is a pointer")),
        }
    }

    fn get_bottom_layer(&self) -> Vec<usize> {
        self.nodes
            .iter()
            .filter(|node| match node {
                CFFTreeNode::Data(tree_node_data) => tree_node_data.children.is_empty(),
                CFFTreeNode::Pointer(_) => false,
            })
            .map(|node| match node {
                CFFTreeNode::Data(tree_node_data) => tree_node_data.node_id,
                CFFTreeNode::Pointer(_) => unreachable!(),
            })
            .collect_vec()
    }

    fn to_serializable(&self) -> SerializableCFFTree {
        let nodes = self.nodes.iter().map(|n| n.to_serializable()).collect_vec();
        let orientation = self.orientation.to_serializable();
        let term_id = self.term_id;
        let num_data_nodes = self.num_data_nodes;

        SerializableCFFTree {
            nodes,
            orientation,
            term_id,
            num_data_nodes,
        }
    }

    fn from_serializable(serializable: SerializableCFFTree) -> Self {
        let nodes = serializable
            .nodes
            .into_iter()
            .map(CFFTreeNode::from_serializable)
            .collect_vec();

        let term_id = serializable.term_id;
        let num_data_nodes = serializable.num_data_nodes;

        let orientation = Orientation::from_serializable(serializable.orientation);

        Self {
            nodes,
            orientation,
            term_id,
            num_data_nodes,
        }
    }

    fn recursive_eval_from_node<T: FloatLike>(
        &self,
        node_id: usize,
        esurface_cache: &[T],
        node_cache: &mut Vec<Vec<Option<T>>>,
    ) -> T {
        match &self.nodes[node_id] {
            CFFTreeNode::Data(tree_node_data) => {
                let option_esurface_at_node = tree_node_data.esurface_id;
                match option_esurface_at_node {
                    None => {
                        let res = Into::<T>::into(1.);
                        node_cache[self.term_id][node_id] = Some(res);
                        res
                    }
                    Some(esurface_id) => {
                        let res = if !tree_node_data.children.is_empty() {
                            esurface_cache[esurface_id].inv()
                                * (tree_node_data
                                    .children
                                    .iter()
                                    .map(|child_index| {
                                        self.recursive_eval_from_node(
                                            *child_index,
                                            esurface_cache,
                                            node_cache,
                                        )
                                    })
                                    .sum::<T>())
                        } else {
                            esurface_cache[esurface_id].inv()
                        };
                        node_cache[self.term_id][node_id] = Some(res);
                        res
                    }
                }
            }
            CFFTreeNode::Pointer(pointer) => node_cache[pointer.term_id][pointer.node_id].unwrap(),
        }
    }

    fn evaluate_tree<T: FloatLike>(
        &self,
        esurface_cache: &[T],
        node_cache: &mut Vec<Vec<Option<T>>>,
    ) -> T {
        self.recursive_eval_from_node(0, esurface_cache, node_cache)
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct SerializableCFFTree {
    nodes: Vec<SerializableCFFTreeNode>,
    orientation: SerializableOrientation,
    term_id: usize,
    num_data_nodes: usize,
}

#[derive(Debug, Clone)]
struct CFFTreeNodeData {
    node_id: usize,
    graph: CFFIntermediateGraph,
    children: Vec<usize>,
    _parent: Option<usize>,
    esurface_id: Option<usize>,
}

impl CFFTreeNodeData {
    fn to_serializable(&self) -> SerializableCFFTreeNodeData {
        let children = self.children.clone();
        let graph = self.graph.to_serializable();
        let node_id = self.node_id;
        let esurface_id = self.esurface_id;
        let parent = self._parent;

        SerializableCFFTreeNodeData {
            node_id,
            graph,
            children,
            esurface_id,
            parent,
        }
    }

    fn from_serializable(serializable: SerializableCFFTreeNodeData) -> Self {
        let SerializableCFFTreeNodeData {
            node_id,
            graph,
            children,
            esurface_id,
            parent,
        } = serializable;

        Self {
            node_id,
            graph: CFFIntermediateGraph::from_serializable(graph),
            children,
            _parent: parent,
            esurface_id,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
struct CFFTreeNodePointer {
    term_id: usize,
    node_id: usize,
}

#[derive(Debug, Clone)]
enum CFFTreeNode {
    Data(CFFTreeNodeData),
    Pointer(CFFTreeNodePointer),
}

impl CFFTreeNode {
    fn to_serializable(&self) -> SerializableCFFTreeNode {
        match self {
            CFFTreeNode::Data(data) => SerializableCFFTreeNode::Data(data.to_serializable()),
            CFFTreeNode::Pointer(pointer) => SerializableCFFTreeNode::Pointer(*pointer),
        }
    }

    fn from_serializable(serializable: SerializableCFFTreeNode) -> Self {
        match serializable {
            SerializableCFFTreeNode::Data(data) => {
                CFFTreeNode::Data(CFFTreeNodeData::from_serializable(data))
            }
            SerializableCFFTreeNode::Pointer(pointer) => CFFTreeNode::Pointer(pointer),
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct SerializableCFFTreeNodeData {
    node_id: usize,
    graph: SerializableCFFIntermediateGraph,
    children: Vec<usize>,
    parent: Option<usize>,
    esurface_id: Option<usize>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
enum SerializableCFFTreeNode {
    Data(SerializableCFFTreeNodeData),
    Pointer(CFFTreeNodePointer),
}

#[derive(Debug, Clone)]
pub struct CFFExpression {
    pub terms: Vec<CFFTree>,
    pub esurfaces: Vec<Esurface>,
    inequivalent_nodes: HashMap<HashableCFFIntermediateGraph, (usize, usize)>,
}

impl CFFExpression {
    pub fn to_serializable(&self) -> SerializableCFFExpression {
        let terms = self.terms.iter().map(|t| t.to_serializable()).collect_vec();
        let esurfaces = self.esurfaces.clone();
        let inequivalent_nodes = self
            .inequivalent_nodes
            .iter()
            .map(|(k, v)| (k.clone(), v.0, v.1))
            .collect_vec();

        SerializableCFFExpression {
            terms,
            esurfaces,
            inequivalent_nodes,
        }
    }

    pub fn from_serializable(serializable: SerializableCFFExpression) -> Self {
        let terms = serializable
            .terms
            .into_iter()
            .map(CFFTree::from_serializable)
            .collect_vec();

        let esurfaces = serializable.esurfaces;

        let inequivalent_nodes = serializable
            .inequivalent_nodes
            .into_iter()
            .map(|(k, v1, v2)| (k, (v1, v2)))
            .collect();

        Self {
            terms,
            esurfaces,
            inequivalent_nodes,
        }
    }

    #[inline]
    pub fn evaluate_orientations<T: FloatLike>(&self, energy_cache: &[T]) -> Vec<T> {
        let esurface_cache = self.compute_esurface_cache(energy_cache);

        let mut node_cache = self
            .terms
            .iter()
            .map(|t| vec![None; t.nodes.len()])
            .collect_vec();

        self.terms
            .iter()
            .map(|tree| tree.evaluate_tree(&esurface_cache, &mut node_cache))
            .collect()
    }

    #[inline]
    pub fn evaluate<T: FloatLike>(&self, energy_cache: &[T]) -> T {
        self.evaluate_orientations(energy_cache)
            .into_iter()
            .sum::<T>()
    }

    #[inline]
    pub fn compute_esurface_cache<T: FloatLike>(&self, energy_cache: &[T]) -> Vec<T> {
        self.esurfaces
            .iter()
            .map(|e| e.compute_value(energy_cache))
            .collect()
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SerializableCFFExpression {
    terms: Vec<SerializableCFFTree>,
    esurfaces: Vec<Esurface>,
    inequivalent_nodes: Vec<(HashableCFFIntermediateGraph, usize, usize)>,
}
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq, Eq, Hash, Copy)]
enum CFFVertexType {
    Source,
    Sink,
    Both,
}

#[derive(Serialize, Deserialize, Clone, PartialEq, Eq, Hash, Copy)]
struct CFFVertex {
    identifier: [u8; MAX_VERTEX_COUNT],
    len: usize,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct ReadableCFFVertex {
    identifier: Vec<u8>,
}

impl CFFVertex {
    fn from_vec(vec: Vec<u8>) -> Self {
        if vec.len() > MAX_VERTEX_COUNT {
            panic!(
                "Current maximum number of supported vertices is {}",
                MAX_VERTEX_COUNT
            )
        }

        let mut identifier = [0; MAX_VERTEX_COUNT];
        let len = vec.len();

        for (index, value) in vec.into_iter().enumerate() {
            identifier[index] = value;
        }

        Self { identifier, len }
    }

    fn join(&self, other: &CFFVertex) -> Self {
        let mut new_identifier = self.identifier;
        let new_len = self.len + other.len;

        for index in 0..other.len {
            new_identifier[self.len + index] = other.identifier[index];
        }

        Self {
            identifier: new_identifier,
            len: new_len,
        }
    }

    fn sorted_vertex(&self) -> Self {
        let mut new_identifier: Vec<u8> = Vec::from(&self.identifier[0..self.len]);
        new_identifier.sort();
        Self::from_vec(new_identifier)
    }

    fn iter(&self) -> impl Iterator<Item = u8> + '_ {
        self.identifier.iter().take(self.len).copied()
    }

    fn to_serializable(self) -> ReadableCFFVertex {
        let identifier = (0..self.len).map(|i| self.identifier[i]).collect_vec();
        ReadableCFFVertex { identifier }
    }

    fn from_serializable(serializable: ReadableCFFVertex) -> Self {
        Self::from_vec(serializable.identifier)
    }
}

impl Display for CFFVertex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut string = String::from("v(");
        for index in 0..self.len {
            string.push_str(&format!("{},", self.identifier[index]));
        }
        string.pop();
        string.push(')');
        write!(f, "{}", string)
    }
}

impl Debug for CFFVertex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self)
    }
}

// Orientation is a bitset that represents the orientation of the edges in a graph.
// 0 means +, 1 means -, in this representation the original graph is represented by the number 0
#[derive(Debug, Clone, Copy)]
struct Orientation {
    identifier: usize,
    num_edges: usize,
}

impl Orientation {
    #[allow(unused)]
    fn default(num_edges: usize) -> Self {
        Self {
            identifier: 0,
            num_edges,
        }
    }

    fn to_serializable(self) -> SerializableOrientation {
        let orientation = self.into_iter().collect_vec();
        SerializableOrientation { orientation }
    }

    fn from_serializable(serializable: SerializableOrientation) -> Self {
        let num_edges = serializable.orientation.len();

        let identifier = serializable
            .orientation
            .into_iter()
            .enumerate()
            .map(|(index, value)| if value { 0 } else { 1 << index })
            .sum();

        Self {
            identifier,
            num_edges,
        }
    }
}

impl IntoIterator for Orientation {
    type Item = bool;
    type IntoIter = OrientationIterator;

    fn into_iter(self) -> Self::IntoIter {
        OrientationIterator {
            identifier: self.identifier,
            current_location: 0,
            num_edges: self.num_edges,
        }
    }
}

// OrientationIterator allows us to iterate over the edges in a graph, and
// view their orientation as a boolean
struct OrientationIterator {
    identifier: usize,
    current_location: usize,
    num_edges: usize,
}

impl Iterator for OrientationIterator {
    type Item = bool;
    fn next(&mut self) -> Option<Self::Item> {
        if self.current_location < self.num_edges {
            let result = self.identifier & (1 << self.current_location) == 0;
            self.current_location += 1;
            Some(result)
        } else {
            None
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct SerializableOrientation {
    orientation: Vec<bool>,
}

// This function returns an iterator over all possible orientations of a graph
fn iterate_possible_orientations(num_edges: usize) -> impl Iterator<Item = Orientation> {
    if num_edges > 64 {
        panic!("Maximum number of edges supported is currently 64")
    }

    let max_size = 2_usize.pow(num_edges as u32);
    (0..max_size).map(move |x| Orientation {
        identifier: x,
        num_edges,
    })
}

#[derive(Debug, Clone)]
struct CFFIntermediateGraph {
    // outgoing edges, incoming edges
    vertices: HashMap<CFFVertex, (Vec<usize>, Vec<usize>)>,
    edges: HashMap<usize, (CFFVertex, CFFVertex)>,
}

#[allow(unused)]
impl CFFIntermediateGraph {
    fn from_serializable(serializable: SerializableCFFIntermediateGraph) -> Self {
        let mut edges = HashMap::default();
        let mut vertices = HashMap::default();

        for vertex in serializable.vertices.into_iter() {
            let (outgoing, incoming) = (vertex.1, vertex.2);
            vertices.insert(CFFVertex::from_serializable(vertex.0), (outgoing, incoming));
        }

        for edge in serializable.edges.into_iter() {
            edges.insert(
                edge.0,
                (
                    CFFVertex::from_serializable(edge.1),
                    CFFVertex::from_serializable(edge.2),
                ),
            );
        }

        CFFIntermediateGraph { vertices, edges }
    }

    fn to_serializable(&self) -> SerializableCFFIntermediateGraph {
        let vertices = self
            .vertices
            .iter()
            .map(|(v, (o, i))| (v.to_serializable(), o.clone(), i.clone()))
            .collect_vec();

        let edges = self
            .edges
            .iter()
            .map(|(i, (l, r))| (*i, l.to_serializable(), r.to_serializable()))
            .collect_vec();

        SerializableCFFIntermediateGraph { vertices, edges }
    }

    // helper function for testing
    fn from_vec(edges: Vec<(usize, usize)>) -> Self {
        let num_edges = edges.len();
        let max_vertex_num = num_edges + 1;

        let mut vertex_map = HashMap::with_capacity(max_vertex_num);
        let mut edge_map = HashMap::with_capacity(num_edges);

        for (index, edge) in edges.into_iter().enumerate() {
            let left_vertex = CFFVertex::from_vec(vec![edge.0 as u8]);
            let right_vertex = CFFVertex::from_vec(vec![edge.1 as u8]);

            edge_map.insert(index, (left_vertex, right_vertex));
        }

        for (index, edge) in edge_map.iter() {
            let (left_vertex, right_vertex) = edge;
            let (outgoing_edges, incoming_edges) =
                vertex_map.entry(*left_vertex).or_insert((vec![], vec![]));

            outgoing_edges.push(*index);

            let (outgoing_edges, incoming_edges) =
                vertex_map.entry(*right_vertex).or_insert((vec![], vec![]));

            incoming_edges.push(*index);
        }

        CFFIntermediateGraph {
            vertices: vertex_map,
            edges: edge_map,
        }
    }

    // bfs to determine if the graph is connected when vertex is removed
    fn has_connected_complement(&self, vertex: &CFFVertex) -> Result<bool, Report> {
        let vertex_that_is_not_v = self.get_vertex_that_is_not_v(vertex)?;

        let mut current_vertices = HashSet::default();
        current_vertices.insert(vertex_that_is_not_v);

        let mut visited_vertices = HashSet::default();
        visited_vertices.insert(vertex_that_is_not_v);

        let mut vertices_found_in_previous_iteration = HashSet::default();
        let mut delta = 1;

        // stop when no new vertices are found in the previous iteration
        while delta > 0 {
            delta = 0;

            for current_vertex in current_vertices.iter() {
                for all_vertex in self.vertices.keys().filter(|all_v| {
                    self.are_adjacent(current_vertex, all_v).unwrap() && **all_v != *vertex
                }) {
                    // disallow passing along the giving vertex, emulating the graph with vertex removed.
                    if self.are_adjacent(current_vertex, all_vertex)?
                        && *all_vertex != *vertex
                        && visited_vertices.insert(all_vertex)
                    {
                        delta += 1;
                        vertices_found_in_previous_iteration.insert(all_vertex);
                    }
                }
            }

            current_vertices = vertices_found_in_previous_iteration.clone();
            vertices_found_in_previous_iteration.clear();
        }

        // The graph is connected if all vertices have been found
        Ok(visited_vertices.len() == self.vertices.len() - 1)
    }

    // helper fuction for has_connected_complement
    fn get_vertex_that_is_not_v(&self, vertex: &CFFVertex) -> Result<&CFFVertex, Report> {
        for (key, _) in self.vertices.iter() {
            if key != vertex {
                return Ok(key);
            }
        }

        Err(eyre!("Could not find vertex that is not v"))
    }

    fn get_source_sink_candidate_list(&self) -> Result<Vec<(CFFVertex, CFFVertexType)>, Report> {
        let mut res = vec![];
        for (vertex, (outgoing_edges, incoming_edges)) in self
            .vertices
            .iter()
            .filter(|(v, _)| self.has_connected_complement(v).unwrap())
        {
            if outgoing_edges.is_empty() && !incoming_edges.is_empty() {
                res.push((*vertex, CFFVertexType::Sink));
            } else if incoming_edges.is_empty() && !outgoing_edges.is_empty() {
                res.push((*vertex, CFFVertexType::Source));
            }
        }

        if res.is_empty() {
            return Err(eyre!("no source or sink found"));
        }
        Ok(res)
    }

    fn get_source_or_sink(&self) -> Result<(CFFVertex, CFFVertexType), Report> {
        let mut sources_and_sinks = self.get_source_sink_candidate_list()?;

        // sort by the first vertex in the vertex set
        sources_and_sinks.sort_by(|(a, _), (b, _)| a.identifier[0].cmp(&b.identifier[0]));
        Ok(sources_and_sinks[0])
    }

    fn are_adjacent(&self, vertex1: &CFFVertex, vertex2: &CFFVertex) -> Result<bool, Report> {
        let (outgoing, incoming) = self
            .vertices
            .get(vertex1)
            .ok_or_else(|| eyre!("vertex1 not in graph"))?;

        for edge_index in incoming.iter() {
            let (left_vertex, _) = self.edges.get(edge_index).unwrap();

            if *left_vertex == *vertex2 {
                return Ok(true);
            }
        }

        for edge_index in outgoing.iter() {
            let (_, right_vertex) = self.edges.get(edge_index).unwrap();
            if *right_vertex == *vertex2 {
                return Ok(true);
            }
        }

        Ok(false)
    }

    fn contract_from_vertex(
        &self,
        vertex: &CFFVertex,
        vertex_type: CFFVertexType,
    ) -> Result<Vec<(Self, CFFVertex)>, Report> {
        match vertex_type {
            CFFVertexType::Sink => {
                let (_, edges_of_vertex) = self
                    .vertices
                    .get(vertex)
                    .ok_or_else(|| eyre!("vertex not in graph"))?;

                let mut adjacent_vertices = HashMap::default();

                for edge_index in edges_of_vertex.iter() {
                    let (left_vertex, _) = self
                        .edges
                        .get(edge_index)
                        .ok_or_else(|| eyre!("edge not in graph"))?;

                    if !adjacent_vertices.contains_key(left_vertex) {
                        adjacent_vertices.insert(left_vertex, vec![edge_index]);
                    } else {
                        adjacent_vertices
                            .get_mut(left_vertex)
                            .unwrap()
                            .push(edge_index);
                    }
                }

                let mut res = vec![];

                for (adjacent_vertex, edges_to_be_deleted) in adjacent_vertices.iter() {
                    let new_vertex = vertex.join(adjacent_vertex);

                    let (outgoing_edges_of_adjacent_vertex, incoming_edges_of_adjacent_vertex) =
                        self.vertices
                            .get(adjacent_vertex)
                            .ok_or_else(|| eyre!("vertex not in graph"))?;

                    let mut new_vertices = self.vertices.clone();
                    let mut new_edges = self.edges.clone();

                    let outgoing_edges_of_new_vertex = outgoing_edges_of_adjacent_vertex
                        .iter()
                        .filter(|e| !edges_to_be_deleted.contains(e))
                        .copied()
                        .collect_vec();

                    let incoming_edges_of_new_vertex = incoming_edges_of_adjacent_vertex
                        .iter()
                        .chain(edges_of_vertex.iter())
                        .filter(|e| !edges_to_be_deleted.contains(e))
                        .copied()
                        .collect_vec();

                    new_vertices.remove(vertex);
                    new_vertices.remove(adjacent_vertex);
                    new_vertices.insert(
                        new_vertex,
                        (outgoing_edges_of_new_vertex, incoming_edges_of_new_vertex),
                    );

                    for edge_index in edges_to_be_deleted.iter() {
                        new_edges.remove(edge_index);
                    }

                    for (_edge_index, (left_vertex, right_vertex)) in new_edges.iter_mut() {
                        if *left_vertex == **adjacent_vertex || left_vertex == vertex {
                            *left_vertex = new_vertex;
                        }

                        if *right_vertex == **adjacent_vertex || right_vertex == vertex {
                            *right_vertex = new_vertex;
                        }
                    }

                    new_vertices.shrink_to_fit();
                    new_edges.shrink_to_fit();

                    let new_graph = CFFIntermediateGraph {
                        vertices: new_vertices,
                        edges: new_edges,
                    };

                    res.push((new_graph, new_vertex));
                }

                Ok(res)
            }
            CFFVertexType::Source => {
                let (edges_of_vertex, _) = self
                    .vertices
                    .get(vertex)
                    .ok_or_else(|| eyre!("vertex not in graph"))?;

                let mut adjacent_vertices = HashMap::default();

                for edge_index in edges_of_vertex.iter() {
                    let (_, right_vertex) = self
                        .edges
                        .get(edge_index)
                        .ok_or_else(|| eyre!("edge not in graph"))?;

                    if !adjacent_vertices.contains_key(right_vertex) {
                        adjacent_vertices.insert(right_vertex, vec![edge_index]);
                    } else {
                        adjacent_vertices
                            .get_mut(right_vertex)
                            .unwrap()
                            .push(edge_index);
                    }
                }

                let mut res = vec![];

                for (adjacent_vertex, edges_to_be_deleted) in adjacent_vertices.iter() {
                    let new_vertex = vertex.join(adjacent_vertex);

                    let (outgoing_edges_of_adjacent_vertex, incoming_edges_of_adjacent_vertex) =
                        self.vertices
                            .get(adjacent_vertex)
                            .ok_or_else(|| eyre!("vertex not in graph"))?;

                    let mut new_vertices = self.vertices.clone();
                    let mut new_edges = self.edges.clone();

                    let outgoing_edges_of_new_vertex = outgoing_edges_of_adjacent_vertex
                        .iter()
                        .chain(edges_of_vertex.iter())
                        .filter(|e| !edges_to_be_deleted.contains(e))
                        .copied()
                        .collect_vec();

                    let incoming_edges_of_new_vertex = incoming_edges_of_adjacent_vertex
                        .iter()
                        .filter(|e| !edges_to_be_deleted.contains(e))
                        .copied()
                        .collect_vec();

                    new_vertices.remove(vertex);
                    new_vertices.remove(adjacent_vertex);
                    new_vertices.insert(
                        new_vertex,
                        (outgoing_edges_of_new_vertex, incoming_edges_of_new_vertex),
                    );

                    for edge_index in edges_to_be_deleted.iter() {
                        new_edges.remove(edge_index);
                    }

                    for (_edge_index, (left_vertex, right_vertex)) in new_edges.iter_mut() {
                        if *left_vertex == **adjacent_vertex || left_vertex == vertex {
                            *left_vertex = new_vertex;
                        }

                        if *right_vertex == **adjacent_vertex || right_vertex == vertex {
                            *right_vertex = new_vertex;
                        }
                    }

                    new_vertices.shrink_to_fit();
                    new_edges.shrink_to_fit();

                    let new_graph = CFFIntermediateGraph {
                        vertices: new_vertices,
                        edges: new_edges,
                    };

                    res.push((new_graph, new_vertex));
                }

                Ok(res)
            }
            CFFVertexType::Both => Err(eyre!("vertex is not a source or sink")),
        }
    }

    fn has_directed_cycle(&self, seed_vertex: &CFFVertex) -> Result<bool, Report> {
        if !self.vertices.contains_key(seed_vertex) {
            return Err(eyre!("seed vertex not in graph"));
        }

        let mut visited = HashSet::default();
        let mut stack = vec![];
        self.dfs(seed_vertex, &mut visited, &mut stack)
    }

    fn dfs(
        &self,
        vertex: &CFFVertex,
        visited: &mut HashSet<CFFVertex>,
        stack: &mut Vec<CFFVertex>,
    ) -> Result<bool, Report> {
        if visited.contains(vertex) {
            return Ok(stack.contains(vertex));
        }

        visited.insert(*vertex);
        stack.push(*vertex);

        // guard in has_directed_cycle should prevent unwrap() from failing
        let neighbours = self.get_directed_neighbours(vertex)?;
        for neighbour in neighbours.iter() {
            if self.dfs(neighbour, visited, stack)? {
                return Ok(true);
            }
        }

        stack.pop();
        Ok(false)
    }

    fn get_directed_neighbours(&self, vertex: &CFFVertex) -> Result<Vec<CFFVertex>, Report> {
        let (outgoing_edges, _) = self
            .vertices
            .get(vertex)
            .ok_or_else(|| eyre!("vertex not in graph, from directed_neighbours"))?;
        Ok(outgoing_edges
            .iter()
            .map(|e| self.edges.get(e).unwrap().1)
            .collect())
    }

    fn generate_children(
        &self,
        position_map: &HashMap<usize, usize>,
        external_data: &HashMap<usize, Vec<usize>>, // (external vertex, external edges)
        global_orientation: &Orientation,
    ) -> Result<(Option<Vec<Self>>, Esurface), Report> {
        let (vertex_to_contract_from, vertex_type) = self.get_source_or_sink()?;
        // no need to perform the contraction if we only have two vertices, we only need
        // to extract the final esurface.
        let children = if self.vertices.len() > 2 {
            Some(
                self.contract_from_vertex(&vertex_to_contract_from, vertex_type)?
                    .into_iter()
                    .filter(|(graph, new_vertex)| {
                        if graph.get_vertex_type(new_vertex).unwrap() != vertex_type {
                            return !graph.has_directed_cycle(new_vertex).unwrap();
                        }
                        true
                    })
                    .map(|(graph, _)| graph)
                    .collect_vec(),
            )
        } else {
            None
        };

        let mut energies = match vertex_type {
            CFFVertexType::Sink => self
                .vertices
                .get(&vertex_to_contract_from)
                .unwrap()
                .1
                .clone(),
            CFFVertexType::Source => self
                .vertices
                .get(&vertex_to_contract_from)
                .unwrap()
                .0
                .clone(),
            CFFVertexType::Both => unreachable!(),
        };

        energies.sort();
        let global_orientation_vec = global_orientation.into_iter().collect_vec();

        let mut sub_orientation = vec![];
        for energy in energies.iter() {
            sub_orientation.push(global_orientation_vec[*position_map.get(energy).unwrap()])
        }

        let shift = vertex_to_contract_from
            .iter()
            .filter(|v| external_data.contains_key(&(*v as usize)))
            .flat_map(|v| external_data.get(&(v as usize)).unwrap())
            .copied()
            .sorted()
            .collect_vec();

        let shift_signature = match vertex_type {
            CFFVertexType::Source => false,
            CFFVertexType::Sink => true,
            CFFVertexType::Both => unreachable!(),
        };

        let esurface = Esurface {
            energies,
            sub_orientation,
            shift,
            shift_signature,
        };

        Ok((children, esurface))
    }

    fn get_vertex_type(&self, vertex: &CFFVertex) -> Result<CFFVertexType, String> {
        let (outgoing, ingoing) = self
            .vertices
            .get(vertex)
            .ok_or_else(|| String::from("vertex not in graph"))?;
        if outgoing.is_empty() {
            Ok(CFFVertexType::Sink)
        } else if ingoing.is_empty() {
            Ok(CFFVertexType::Source)
        } else {
            Ok(CFFVertexType::Both)
        }
    }

    fn has_directed_cycle_initial(&self) -> Result<bool, Report> {
        // in the first iteration we can not just check one vertex.
        // this can probably be improved by skipping the visited vertices in the next iteration
        for vertex in self.vertices.keys() {
            if self.has_directed_cycle(vertex)? {
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn get_sorted_edge_list(&self) -> Vec<usize> {
        self.edges.keys().copied().sorted().collect_vec()
    }

    fn to_hashable(&self) -> HashableCFFIntermediateGraph {
        let sorted_keys = self.edges.keys().copied().sorted().collect_vec();
        let sorted_edges = sorted_keys
            .into_iter()
            .map(|k| {
                let (left_vertex, right_vertex) = self.edges.get(&k).unwrap();

                (k, left_vertex.sorted_vertex(), right_vertex.sorted_vertex())
            })
            .collect();
        HashableCFFIntermediateGraph {
            edges: sorted_edges,
        }
    }
}

impl Display for CFFIntermediateGraph {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut string = String::from("\nvertices: \n");
        for (vertex, (outgoing_edges, incoming_edges)) in self.vertices.iter() {
            string.push_str(&format!(
                "\t{}: outgoing: {:?}, incoming: {:?}\n",
                vertex, outgoing_edges, incoming_edges
            ));
        }

        string.push_str("edges: \n");
        for (edge, (left_vertex, right_vertex)) in self.edges.iter() {
            string.push_str(&format!(
                "\t{}: left: {}, right: {}\n",
                edge, left_vertex, right_vertex
            ));
        }

        write!(f, "{}", string)
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct SerializableCFFIntermediateGraph {
    edges: Vec<(usize, ReadableCFFVertex, ReadableCFFVertex)>,
    vertices: Vec<(ReadableCFFVertex, Vec<usize>, Vec<usize>)>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct HashableCFFIntermediateGraph {
    edges: Vec<(usize, CFFVertex, CFFVertex)>,
}

impl Hash for HashableCFFIntermediateGraph {
    fn hash<H: Hasher>(&self, state: &mut H) {
        for e in &self.edges {
            e.0.hash(state);
            e.1.sorted_vertex().hash(state);
            e.2.sorted_vertex().hash(state);
        }
    }
}

impl PartialEq for HashableCFFIntermediateGraph {
    fn eq(&self, other: &Self) -> bool {
        if self.edges.len() != other.edges.len() {
            return false;
        }

        for (self_edge, other_edge) in self.edges.iter().zip(other.edges.iter()) {
            if self_edge.0 != other_edge.0
                || self_edge.1.sorted_vertex() != other_edge.1.sorted_vertex()
                || self_edge.2.sorted_vertex() != other_edge.2.sorted_vertex()
            {
                return false;
            }
        }

        true
    }
}

impl Eq for HashableCFFIntermediateGraph {}

fn get_orientations(
    graph: &Graph,
) -> (
    Vec<(Orientation, CFFIntermediateGraph)>,
    HashMap<usize, usize>,
) {
    let virtual_edges = graph
        .get_virtual_edges_iterator()
        .map(|(index, _graph_edge)| index)
        .sorted()
        .collect_vec();

    let num_edges = virtual_edges.len();

    // maps the position of the edge in the complete Graph to the position in the CFFIntermediateGraph
    let mut position_map = HashMap::default();
    for (index, edge) in virtual_edges.iter().enumerate() {
        position_map.insert(*edge, index);
    }

    (
        iterate_possible_orientations(num_edges)
            .map(|orientation| {
                // map orientations onto graph
                let mut vertices: HashMap<CFFVertex, (Vec<usize>, Vec<usize>)> = HashMap::default();
                let mut edges = HashMap::default();

                for (edge_orientation, edge_id) in orientation.into_iter().zip(virtual_edges.iter())
                {
                    let edge_vertices = graph.edges[*edge_id].vertices;
                    let left_vertex;
                    let right_vertex;

                    if edge_orientation {
                        left_vertex = CFFVertex::from_vec(vec![edge_vertices[0] as u8]);
                        right_vertex = CFFVertex::from_vec(vec![edge_vertices[1] as u8]);
                    } else {
                        left_vertex = CFFVertex::from_vec(vec![edge_vertices[1] as u8]);
                        right_vertex = CFFVertex::from_vec(vec![edge_vertices[0] as u8]);
                    }

                    edges.insert(*edge_id, (left_vertex, right_vertex));

                    if let std::collections::hash_map::Entry::Vacant(e) =
                        vertices.entry(left_vertex)
                    {
                        e.insert((vec![*edge_id], vec![]));
                    } else {
                        vertices.get_mut(&left_vertex).unwrap().0.push(*edge_id);
                    }

                    if let std::collections::hash_map::Entry::Vacant(e) =
                        vertices.entry(right_vertex)
                    {
                        e.insert((vec![], vec![*edge_id]));
                    } else {
                        vertices.get_mut(&right_vertex).unwrap().1.push(*edge_id);
                    }
                }

                let graph = CFFIntermediateGraph { vertices, edges };
                (orientation, graph)
            })
            .collect_vec(),
        position_map,
    )
}

pub fn generate_cff_expression(graph: &Graph) -> Result<CFFExpression, Report> {
    // construct a hashmap that contains as keys all vertices that connect to external edges
    // and as values those external edges that it connects to
    let mut external_data: HashMap<usize, Vec<usize>> = HashMap::default();

    for external_edge in graph.edges.iter() {
        let edge_position = graph.get_edge_position(&external_edge.name).unwrap();
        match external_edge.edge_type {
            EdgeType::Incoming => {
                if let std::collections::hash_map::Entry::Vacant(e) =
                    external_data.entry(external_edge.vertices[1])
                {
                    e.insert(vec![edge_position]);
                } else {
                    external_data
                        .get_mut(&external_edge.vertices[1])
                        .unwrap_or_else(|| unreachable!())
                        .push(edge_position);
                }
            }
            EdgeType::Outgoing => {
                if let std::collections::hash_map::Entry::Vacant(e) =
                    external_data.entry(external_edge.vertices[0])
                {
                    e.insert(vec![edge_position]);
                } else {
                    external_data
                        .get_mut(&external_edge.vertices[0])
                        .unwrap_or_else(|| unreachable!())
                        .push(edge_position);
                }
            }

            EdgeType::Virtual => (),
        }
    }

    let (orientations, position_map) = get_orientations(graph);
    debug!("generating cff for graph: {}", graph.name);
    debug!("number of orientations: {}", orientations.len());

    generate_cff_from_orientations(orientations, &position_map, &external_data)
}

fn generate_cff_from_orientations(
    orientations_and_graphs: Vec<(Orientation, CFFIntermediateGraph)>,
    position_map: &HashMap<usize, usize>,
    external_data: &HashMap<usize, Vec<usize>>,
) -> Result<CFFExpression, Report> {
    let mut cff_expression = CFFExpression {
        terms: vec![],
        esurfaces: vec![],
        inequivalent_nodes: HashMap::default(),
    };

    // filter cyclic orientations beforehand
    let acyclic_orientations_and_graphs = orientations_and_graphs
        .into_iter()
        .filter(|(_, graph)| !graph.has_directed_cycle_initial().unwrap())
        .collect_vec();

    debug!(
        "number of acyclic orientations: {}",
        acyclic_orientations_and_graphs.len()
    );

    let mut cache_hits = 0;
    let mut non_cache_hits = 0;
    for (term_id, (orientation, graph)) in acyclic_orientations_and_graphs.into_iter().enumerate() {
        let mut tree = CFFTree::from_root_graph(graph, orientation, term_id);
        let mut tree_done = false;

        while !tree_done {
            let bottom_layer = tree.get_bottom_layer();
            if bottom_layer.is_empty() {
                break;
            }

            for node_id in bottom_layer.into_iter() {
                let node = match &tree.nodes[node_id] {
                    CFFTreeNode::Data(tree_node_data) => tree_node_data,
                    _ => unreachable!(), // this is impossible by definition of get_bottom_layer()
                };

                let (option_children, esurface) =
                    node.graph
                        .generate_children(position_map, external_data, &orientation)?;

                if let Some(esurface_id) =
                    cff_expression.esurfaces.iter().position(|e| e == &esurface)
                {
                    tree.insert_esurface(node_id, esurface_id)?;
                } else {
                    tree.insert_esurface(node_id, cff_expression.esurfaces.len())?;
                    cff_expression.esurfaces.push(esurface);
                }

                if let Some(children) = option_children {
                    for child in children.into_iter() {
                        let hashable_child = child.to_hashable();
                        if let Some((cff_expression_term_id, cff_expression_node_id)) =
                            cff_expression.inequivalent_nodes.get(&hashable_child)
                        {
                            let new_pointer = CFFTreeNodePointer {
                                term_id: *cff_expression_term_id,
                                node_id: *cff_expression_node_id,
                            };
                            tree.insert_pointer(node_id, new_pointer)?;
                            cache_hits += 1;
                        } else {
                            let child_node_id = tree.nodes.len();

                            cff_expression
                                .inequivalent_nodes
                                .insert(hashable_child, (term_id, child_node_id));
                            tree.insert_graph(node_id, child)?;
                            non_cache_hits += 1;
                        }
                    }
                } else {
                    tree_done = true;
                }
            }
        }

        cff_expression.terms.push(tree);
    }

    debug!("number of cache hits: {}", cache_hits);
    debug!(
        "percentage of cache hits: {:.1}%",
        cache_hits as f64 / (cache_hits + non_cache_hits) as f64 * 100.0
    );

    Ok(cff_expression)
}

#[cfg(test)]
mod tests_cff {
    use lorentz_vector::LorentzVector;
    use num::traits::Inv;

    use super::*;

    // helper function to do some quick tests
    #[allow(unused)]
    fn generate_orientations_for_testing(
        edges: Vec<(usize, usize)>,
    ) -> Vec<(Orientation, CFFIntermediateGraph)> {
        let num_edges = edges.len();

        iterate_possible_orientations(num_edges)
            .map(|or| {
                let mut new_edges = edges.clone();
                for (edge_id, edge_orientation) in or.into_iter().enumerate() {
                    if edge_orientation {
                        new_edges[edge_id] = edges[edge_id];
                    } else {
                        let rotated_edge = (edges[edge_id].1, edges[edge_id].0);
                        new_edges[edge_id] = rotated_edge;
                    }
                }

                let new_graph = CFFIntermediateGraph::from_vec(new_edges);
                (or, new_graph)
            })
            .filter(|(or, graph)| !graph.has_directed_cycle_initial().unwrap())
            .collect_vec()
    }

    #[allow(unused)]
    fn compute_one_loop_energy<T: FloatLike>(k: LorentzVector<T>, p: LorentzVector<T>, m: T) -> T {
        ((k + p).spatial_squared() + m * m).sqrt()
    }

    #[test]
    fn test_from_vec() {
        let triangle_edge_vec = vec![(0, 1), (1, 2), (2, 0)];
        let test_struct = CFFIntermediateGraph::from_vec(triangle_edge_vec);

        // test if test_struct contains edges, 0, 1, 2 and if they have the correct vertices
        for index in 0..3 {
            let (left_vertex, right_vertex) = test_struct.edges.get(&index).unwrap();
            assert_eq!(left_vertex, &CFFVertex::from_vec(vec![index as u8]));
            assert_eq!(
                right_vertex,
                &CFFVertex::from_vec(vec![((index + 1) % 3) as u8])
            );
        }

        // test if test_struct contains vertices 0, 1, 2 and if they have the correct edges
        for index in 0..3 {
            let (outgoing_edges, incoming_edges) = test_struct
                .vertices
                .get(&CFFVertex::from_vec(vec![index]))
                .unwrap();

            assert_eq!(outgoing_edges.len(), 1);
            assert_eq!(incoming_edges.len(), 1);

            assert_eq!(outgoing_edges[0] as u8, index);
            assert_eq!(incoming_edges[0] as u8, (index + 2) % 3);
        }
    }

    #[test]
    fn test_vertex_that_is_not_v() {
        let triangle_edge_vec = vec![(0, 1), (1, 2), (2, 0)];

        let cff_test_struct = CFFIntermediateGraph::from_vec(triangle_edge_vec);

        let vertex = CFFVertex::from_vec(vec![0]);

        let vertex_that_is_not_v = cff_test_struct.get_vertex_that_is_not_v(&vertex).unwrap();
        assert_ne!(*vertex_that_is_not_v, vertex);
    }

    #[test]
    fn test_has_connected_complement() {
        let triangle_edge_vec = vec![(0, 1), (1, 2), (2, 0)];
        let cff_test_struct1 = CFFIntermediateGraph::from_vec(triangle_edge_vec);
        assert!(cff_test_struct1
            .has_connected_complement(&CFFVertex::from_vec(vec![0]))
            .unwrap());

        let double_bubble_vec = vec![(0, 1), (0, 1), (1, 2), (1, 2)];
        let cff_test_struct2 = CFFIntermediateGraph::from_vec(double_bubble_vec);

        assert!(!cff_test_struct2
            .has_connected_complement(&CFFVertex::from_vec(vec![1]))
            .unwrap());
    }

    #[test]
    fn test_are_adjacent() {
        let triangle_edge_vec = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        let cff_test_struct1 = CFFIntermediateGraph::from_vec(triangle_edge_vec);

        assert!(cff_test_struct1
            .are_adjacent(&CFFVertex::from_vec(vec![0]), &CFFVertex::from_vec(vec![1]))
            .unwrap());

        assert!(!cff_test_struct1
            .are_adjacent(&CFFVertex::from_vec(vec![0]), &CFFVertex::from_vec(vec![2]))
            .unwrap());
    }

    #[test]
    fn test_source_sink_candidate_list() {
        let bubble_edge_vec = vec![(0, 1), (0, 1)];

        let cff_test_struct1 = CFFIntermediateGraph::from_vec(bubble_edge_vec);

        let source_sink_canditate_list = cff_test_struct1.get_source_sink_candidate_list().unwrap();

        assert_eq!(source_sink_canditate_list.len(), 2);
        assert!(source_sink_canditate_list
            .contains(&(CFFVertex::from_vec(vec![0]), CFFVertexType::Source)));
        assert!(source_sink_canditate_list
            .contains(&(CFFVertex::from_vec(vec![1]), CFFVertexType::Sink)));

        let triangle_edge_vec = vec![(0, 1), (1, 2), (2, 0)];

        let source_sink_canditate_list =
            CFFIntermediateGraph::from_vec(triangle_edge_vec).get_source_sink_candidate_list();

        assert!(
            source_sink_canditate_list.is_err(),
            "{:?}",
            source_sink_canditate_list
        );

        let double_bubble = vec![(0, 1), (0, 1), (2, 1), (2, 1)];
        let cff_test_struct3 = CFFIntermediateGraph::from_vec(double_bubble);
        let list = cff_test_struct3.get_source_sink_candidate_list().unwrap();
        assert!(list.contains(&(CFFVertex::from_vec(vec![0]), CFFVertexType::Source)));
        assert!(list.contains(&(CFFVertex::from_vec(vec![2]), CFFVertexType::Source)));
        assert!(!list.contains(&(CFFVertex::from_vec(vec![1]), CFFVertexType::Sink)));
    }

    #[test]
    fn test_has_directed_cycle() {
        let triangle_edge_vec = vec![(0, 1), (1, 2), (2, 0)];

        let cff_test_struct1 = CFFIntermediateGraph::from_vec(triangle_edge_vec);

        assert!(cff_test_struct1
            .has_directed_cycle(&CFFVertex::from_vec(vec![0]))
            .unwrap());

        let triangle_edge_vec = vec![(0, 1), (1, 2), (0, 2)];

        let cff_test_struct1 = CFFIntermediateGraph::from_vec(triangle_edge_vec);

        assert!(!cff_test_struct1
            .has_directed_cycle(&CFFVertex::from_vec(vec![0]))
            .unwrap());
    }

    #[test]
    fn test_contract_from_vertex() {
        let triangle_edge_vec = vec![(0, 1), (1, 2), (0, 2)];

        let cff_test_struct1 = CFFIntermediateGraph::from_vec(triangle_edge_vec);

        let res = cff_test_struct1
            .contract_from_vertex(&CFFVertex::from_vec(vec![0]), CFFVertexType::Source)
            .unwrap();

        assert_eq!(res.len(), 2);
        for (result, _new_vertex) in res.iter() {
            assert_eq!(result.vertices.len(), 2);
            assert_eq!(result.edges.len(), 2);
            assert!(
                (result
                    .vertices
                    .contains_key(&CFFVertex::from_vec(vec![0, 2]))
                    && result.vertices.contains_key(&CFFVertex::from_vec(vec![1])))
                    || (result
                        .vertices
                        .contains_key(&CFFVertex::from_vec(vec![0, 1]))
                        && result.vertices.contains_key(&CFFVertex::from_vec(vec![2])))
            );
        }

        let edge_vec = vec![(0, 1), (0, 1), (1, 2), (2, 3), (3, 1)];

        let cff_struct_2 = CFFIntermediateGraph::from_vec(edge_vec);

        let res = cff_struct_2
            .contract_from_vertex(&CFFVertex::from_vec(vec![0]), CFFVertexType::Source)
            .unwrap();

        assert_eq!(res.len(), 1);
        let (res_0, _new_vertex) = res[0].clone();
        assert_eq!(res_0.vertices.len(), 3);
        assert_eq!(res_0.edges.len(), 3);
        assert!(res_0
            .vertices
            .contains_key(&CFFVertex::from_vec(vec![0, 1])));

        assert!(res_0.vertices.contains_key(&CFFVertex::from_vec(vec![2])));
        assert!(res_0.vertices.contains_key(&CFFVertex::from_vec(vec![3])));

        let edge_vec = vec![(1, 0), (1, 0), (1, 2), (2, 3), (3, 1)];

        let cff_struct_2 = CFFIntermediateGraph::from_vec(edge_vec);

        let res = cff_struct_2
            .contract_from_vertex(&CFFVertex::from_vec(vec![0]), CFFVertexType::Sink)
            .unwrap();

        assert_eq!(res.len(), 1);
        let (res_0, _new_vertex) = res[0].clone();
        assert_eq!(res_0.vertices.len(), 3);
        assert_eq!(res_0.edges.len(), 3);
        assert!(res_0
            .vertices
            .contains_key(&CFFVertex::from_vec(vec![0, 1])));

        assert!(res_0.vertices.contains_key(&CFFVertex::from_vec(vec![2])));
        assert!(res_0.vertices.contains_key(&CFFVertex::from_vec(vec![3])));
    }

    #[test]
    fn test_orientation_struct() {
        let orientations = iterate_possible_orientations(3).collect_vec();
        assert_eq!(orientations.len(), 8);

        let orientation1 = orientations[0].into_iter().collect_vec();
        assert_eq!(orientation1, vec![true, true, true]);

        let orientation2 = orientations[1].into_iter().collect_vec();
        assert_eq!(orientation2, vec![false, true, true]);

        let orientation3 = orientations[2].into_iter().collect_vec();
        assert_eq!(orientation3, vec![true, false, true]);

        let orientation4 = orientations[3].into_iter().collect_vec();
        assert_eq!(orientation4, vec![false, false, true]);

        let orientation5 = orientations[4].into_iter().collect_vec();
        assert_eq!(orientation5, vec![true, true, false]);

        let orientation6 = orientations[5].into_iter().collect_vec();
        assert_eq!(orientation6, vec![false, true, false]);

        let orientation7 = orientations[6].into_iter().collect_vec();
        assert_eq!(orientation7, vec![true, false, false]);

        let orientation8 = orientations[7].into_iter().collect_vec();
        assert_eq!(orientation8, vec![false, false, false]);
    }

    #[test]
    fn test_tree_structure() {
        let triangle_edge_vec = vec![(0, 1), (1, 2), (0, 2)];

        let stupid_test_graph = CFFIntermediateGraph::from_vec(triangle_edge_vec);

        let test_tree = CFFTree::from_root_graph(stupid_test_graph, Orientation::default(3), 0);
        let bottom_layer = test_tree.get_bottom_layer();
        assert_eq!(bottom_layer.len(), 1);
        // needs more tests
    }

    #[test]
    fn test_get_directed_neighbours() {
        let triangle_edge_vec = vec![(0, 1), (1, 2), (0, 2)];
        let test_graph = CFFIntermediateGraph::from_vec(triangle_edge_vec);

        let neighbours = test_graph
            .get_directed_neighbours(&CFFVertex::from_vec(vec![0]))
            .unwrap();

        assert!(neighbours.len() == 2);
        assert!(neighbours.contains(&CFFVertex::from_vec(vec![1])));
        assert!(neighbours.contains(&CFFVertex::from_vec(vec![2])));

        let neighbours = test_graph
            .get_directed_neighbours(&CFFVertex::from_vec(vec![1]))
            .unwrap();

        assert!(neighbours.len() == 1);
        assert!(neighbours.contains(&CFFVertex::from_vec(vec![2])));

        let neighbours = test_graph
            .get_directed_neighbours(&CFFVertex::from_vec(vec![2]))
            .unwrap();

        assert!(neighbours.is_empty());
    }

    #[test]
    fn test_graph_equality() {
        let triangle_edge_vec = vec![(0, 1), (1, 2), (0, 2)];
        let triangle_edge_vec2 = vec![(0, 1), (1, 2), (0, 2)];

        let test_graph1 = CFFIntermediateGraph::from_vec(triangle_edge_vec);
        let test_graph2 = CFFIntermediateGraph::from_vec(triangle_edge_vec2);

        assert_eq!(test_graph1.to_hashable(), test_graph2.to_hashable());

        let triangle_edge_vec3 = vec![(0, 1), (1, 2), (2, 0)];
        let test_graph3 = CFFIntermediateGraph::from_vec(triangle_edge_vec3);
        assert_ne!(test_graph1.to_hashable(), test_graph3.to_hashable());

        let more_edge_vec = vec![(0, 1), (2, 1), (0, 2)];
        let test_graph1 = CFFIntermediateGraph::from_vec(more_edge_vec);
        let test_graph_1_contract_from_0 = test_graph1
            .contract_from_vertex(&CFFVertex::from_vec(vec![0]), CFFVertexType::Source)
            .unwrap()
            .into_iter()
            .map(|(g, _e)| g)
            .filter(|g| match g.has_directed_cycle_initial() {
                Ok(b) => !b,
                Err(err) => panic!("{:?}", err),
            })
            .collect_vec();

        let triangle_edge_vec4 = vec![(0, 1), (2, 1), (2, 0)];
        let test_graph4 = CFFIntermediateGraph::from_vec(triangle_edge_vec4);
        let test_graph4_contract_from_2 = test_graph4
            .contract_from_vertex(&CFFVertex::from_vec(vec![2]), CFFVertexType::Source)
            .unwrap()
            .into_iter()
            .map(|(g, _e)| g)
            .filter(|g| match g.has_directed_cycle_initial() {
                Ok(b) => !b,
                Err(err) => panic!("{:?}", err),
            })
            .collect_vec();

        assert_eq!(test_graph_1_contract_from_0.len(), 1);
        assert_eq!(test_graph4_contract_from_2.len(), 1);
        assert_eq!(
            test_graph_1_contract_from_0[0].to_hashable(),
            test_graph4_contract_from_2[0].to_hashable(),
        );
    }

    #[test]
    fn test_esurface() {
        let energies_cache = [1., 2., 3., 4., 5.];
        let shift = vec![3, 4];
        let energies = vec![0, 1, 2];
        let shift_signature = true;
        let sub_orientation = vec![true, false, true];

        let esurface = Esurface {
            sub_orientation: sub_orientation.clone(),
            energies,
            shift,
            shift_signature,
        };

        let _edge_types = [
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Incoming,
            EdgeType::Incoming,
        ];

        let res = esurface.compute_value(&energies_cache);
        assert_eq!(res, 15.);

        let shift_signature = false;
        let energies = vec![0, 2];
        let shift = vec![1];

        let esurface = Esurface {
            energies,
            sub_orientation: sub_orientation.clone(),
            shift,
            shift_signature,
        };

        let res = esurface.compute_value(&energies_cache);
        assert_eq!(res, 2.);
    }

    #[test]
    fn test_hashable_graph() {
        let triangle_edgfes = vec![(0, 1), (1, 2), (2, 0)];
        let cff_struct = CFFIntermediateGraph::from_vec(triangle_edgfes);

        println!("{:?}", cff_struct);
        println!("{:?}", cff_struct.to_hashable());

        let hashable = cff_struct.to_hashable();
        let vertex1 = CFFVertex::from_vec(vec![0]);
        let vertex2 = CFFVertex::from_vec(vec![1]);
        let vertex3 = CFFVertex::from_vec(vec![2]);

        assert_eq!(vertex1.sorted_vertex(), vertex1);
        assert_eq!(hashable.edges.len(), 3);
        assert_eq!(hashable.edges[0], (0, vertex1, vertex2));
        assert_eq!(hashable.edges[1], (1, vertex2, vertex3));
        assert_eq!(hashable.edges[2], (2, vertex3, vertex1));
    }
    #[test]
    fn test_cff_generation_triangle() {
        let triangle = vec![(2, 0), (0, 1), (1, 2)];
        let orientations = generate_orientations_for_testing(triangle);
        let mut external_data = HashMap::default();
        external_data.insert(0, vec![3]);
        external_data.insert(1, vec![4]);
        external_data.insert(2, vec![5]);
        assert_eq!(orientations.len(), 6);

        let mut position_map = HashMap::default();
        for i in 0..3 {
            position_map.insert(i, i);
        }

        let cff =
            generate_cff_from_orientations(orientations, &position_map, &external_data).unwrap();
        assert_eq!(cff.esurfaces.len(), 6);

        let p1 = LorentzVector::from_args(1., 3., 4., 5.);
        let p2 = LorentzVector::from_args(1., 6., 7., 8.);
        let p3 = -p1 - p2;
        let zero = LorentzVector::from_args(0., 0., 0., 0.);
        let m = 0.;

        let k = LorentzVector::from_args(0., 1., 2., 3.);

        let virtual_energy_cache = [
            compute_one_loop_energy(k, zero, m),
            compute_one_loop_energy(k, p1, m),
            compute_one_loop_energy(k, p1 + p2, m),
        ];

        let external_energy_cache = [p1.t, p2.t, p3.t];

        // combine the virtual and external energies
        let mut energy_cache = virtual_energy_cache.to_vec();
        energy_cache.extend(external_energy_cache);

        let _edge_types = [
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Incoming,
            EdgeType::Incoming,
            EdgeType::Incoming,
        ];

        let energy_prefactor = virtual_energy_cache
            .iter()
            .map(|e| (2. * e).inv())
            .product::<f64>();

        let cff_res: f64 =
            energy_prefactor * cff.evaluate(&energy_cache) * (2. * std::f64::consts::PI).powi(-3);

        let target_res = 6.333_549_225_536_17e-9_f64;
        let absolute_error: f64 = cff_res - target_res;
        let relative_error = absolute_error.abs() / cff_res.abs();

        assert!(
            relative_error.abs() < 1.0e-15,
            "relative error: {:+e} (ground truth: {:+e} vs reproduced: {:+e})",
            relative_error,
            target_res,
            cff_res
        );
    } //

    #[test]
    fn test_cff_test_double_triangle() {
        let double_triangle_edges = vec![(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)];

        let mut external_data = HashMap::default();
        external_data.insert(0, vec![5]);
        external_data.insert(3, vec![6]);

        let mut position_map = HashMap::default();
        for i in 0..5 {
            position_map.insert(i, i);
        }

        let orientations = generate_orientations_for_testing(double_triangle_edges);
        let cff =
            generate_cff_from_orientations(orientations, &position_map, &external_data).unwrap();

        for tree in cff.terms.iter() {
            for node in tree.nodes.iter() {
                println!("{:?}", node);
            }
        }

        let _edge_types = [
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Incoming,
            EdgeType::Incoming,
        ];

        let q = LorentzVector::from_args(1., 2., 3., 4.);
        let zero = LorentzVector::from_args(0., 0., 0., 0.);

        let k = LorentzVector::from_args(3., 6., 23., 9.);
        let l = LorentzVector::from_args(0., 3., 12., 34.);

        let virtual_energy_cache = [
            compute_one_loop_energy(k, zero, 0.),
            compute_one_loop_energy(q - k, zero, 0.),
            compute_one_loop_energy(k - l, zero, 0.),
            compute_one_loop_energy(l, zero, 0.),
            compute_one_loop_energy(q - l, zero, 0.),
        ];

        let external_energy_cache = [q.t, -q.t];

        let mut energy_cache = virtual_energy_cache.to_vec();
        energy_cache.extend(external_energy_cache);

        let energy_prefactor = virtual_energy_cache
            .iter()
            .map(|e| (2. * e).inv())
            .product::<f64>();
        let cff_res = energy_prefactor * cff.evaluate(&energy_cache);

        let target = 1.0794792137096797e-13;
        let absolute_error = cff_res - target;
        let relative_error = absolute_error / cff_res;

        assert!(
            relative_error.abs() < 1.0e-15,
            "relative error: {:+e}, target: {:+e}, result: {:+e}",
            relative_error,
            target,
            cff_res
        );
    }

    #[test]
    fn test_cff_tbt() {
        let tbt_edges = vec![
            (0, 1),
            (2, 0),
            (1, 2),
            (1, 3),
            (2, 4),
            (3, 4),
            (3, 5),
            (5, 4),
        ];

        let _edge_types = [
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Incoming,
            EdgeType::Incoming,
        ];

        let mut external_data = HashMap::default();
        external_data.insert(0, vec![8]);
        external_data.insert(5, vec![9]);

        let mut position_map = HashMap::default();
        for i in 0..tbt_edges.len() {
            position_map.insert(i, i);
        }

        let orientataions = generate_orientations_for_testing(tbt_edges);
        let cff =
            generate_cff_from_orientations(orientataions, &position_map, &external_data).unwrap();

        let q = LorentzVector::from_args(1.0, 2.0, 3.0, 4.0);
        let zero_vector = LorentzVector::from_args(0., 0., 0., 0.);

        let p0 = q;
        let p5 = -q;

        let k = LorentzVector::from_args(3., 6., 23., 9.);
        let l = LorentzVector::from_args(0., 3., 12., 34.);
        let m = LorentzVector::from_args(0., 7., 24., 1.);

        let mass = 0.;

        let energies_cache = [
            compute_one_loop_energy(k, zero_vector, mass),
            compute_one_loop_energy(k - q, zero_vector, mass),
            compute_one_loop_energy(k - l, zero_vector, mass),
            compute_one_loop_energy(l, zero_vector, mass),
            compute_one_loop_energy(q - l, zero_vector, mass),
            compute_one_loop_energy(l - m, zero_vector, mass),
            compute_one_loop_energy(m, zero_vector, mass),
            compute_one_loop_energy(m - q, zero_vector, mass),
            p0.t,
            p5.t,
        ];

        let virtual_energy_cache = energies_cache[0..8].to_vec();

        let energy_prefactor = virtual_energy_cache
            .iter()
            .map(|e| (2. * e).inv())
            .product::<f64>();
        let res = cff.evaluate(&energies_cache) * energy_prefactor;

        let absolute_error = res - 1.2625322619777278e-21;
        let relative_error = absolute_error / res;
        assert!(
            relative_error.abs() < 1.0e-15,
            "relative error: {:+e}",
            relative_error
        );
    }

    #[test]
    #[ignore]
    fn fishnet2b2() {
        let edges = vec![
            (0, 1),
            (1, 2),
            (3, 4),
            (4, 5),
            (6, 7),
            (7, 8),
            (0, 3),
            (1, 4),
            (2, 5),
            (3, 6),
            (4, 7),
            (5, 8),
        ];

        let _edge_types = [
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Virtual,
            EdgeType::Incoming,
            EdgeType::Incoming,
            EdgeType::Incoming,
            EdgeType::Incoming,
        ];

        let mut external_data = HashMap::default();
        external_data.insert(0, vec![12]);
        external_data.insert(2, vec![13]);
        external_data.insert(6, vec![14]);
        external_data.insert(8, vec![15]);

        let mut position_map = HashMap::default();
        for i in 0..edges.len() {
            position_map.insert(i, i);
        }

        let orientations = generate_orientations_for_testing(edges);

        // get time before cff generation
        let start = std::time::Instant::now();

        let cff =
            generate_cff_from_orientations(orientations, &position_map, &external_data).unwrap();
        let num_terms = cff.terms.len();
        let num_nodes = cff.terms.iter().map(|t| t.nodes.len()).sum::<usize>();

        let finish = std::time::Instant::now();
        println!("time to generate cff: {:?}", finish - start);
        println!("number of cff terms: {}", num_terms);
        println!("number of nodes: {}", num_nodes);
    }

    #[test]
    #[ignore]
    fn cube() {
        let edges = vec![
            (0, 1),
            (1, 3),
            (3, 2),
            (2, 0),
            (4, 5),
            (5, 7),
            (7, 6),
            (6, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];

        let mut external_data = HashMap::default();
        for v in 0..8 {
            external_data.insert(v, vec![12 + v]);
        }

        let mut position_map = HashMap::default();
        for i in 0..edges.len() {
            position_map.insert(i, i);
        }

        let orientations = generate_orientations_for_testing(edges);

        // get time before cff generation
        let _start = std::time::Instant::now();

        let cff =
            generate_cff_from_orientations(orientations, &position_map, &external_data).unwrap();
        let _num_terms = cff.terms.len();
        let _num_nodes = cff.terms.iter().map(|t| t.nodes.len()).sum::<usize>();

        let _finish = std::time::Instant::now();
    }

    #[test]
    #[ignore] // this is now in the python tests
    fn fishnet2b3() {
        let edges = vec![
            (0, 1),
            (1, 2),
            (3, 4),
            (4, 5),
            (6, 7),
            (7, 8),
            (9, 10),
            (10, 11),
            (0, 3),
            (3, 6),
            (6, 9),
            (1, 4),
            (4, 7),
            (7, 10),
            (2, 5),
            (5, 8),
            (8, 11),
        ];

        // vector with 17 virtuals and 4 incoming
        let mut edge_types = vec![];
        for _i in 0..17 {
            edge_types.push(EdgeType::Virtual);
        }

        for _i in 0..4 {
            edge_types.push(EdgeType::Incoming);
        }

        let external_data = HashMap::default();
        let mut position_map = HashMap::default();
        for i in 0..edges.len() {
            position_map.insert(i, i);
        }
        let _energy_cache = [3.0; 17];

        println!("generating orientations");
        let orientations = generate_orientations_for_testing(edges);
        println!("orientations generated");
        let cff =
            generate_cff_from_orientations(orientations, &position_map, &external_data).unwrap();
        println!("cff terms = {}", cff.terms.len());
    }
}
