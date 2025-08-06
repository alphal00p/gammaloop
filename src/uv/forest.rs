use crate::{
    cff::{
        expression::{GraphOrientation, OrientationData, OrientationID},
        generation::ShiftRewrite,
    },
    new_graph::{Edge, Vertex},
    uv::approx::CFFapprox,
};
use bitvec::vec::BitVec;
use symbolica::atom::Atom;

use linnet::half_edge::{involution::EdgeIndex, subgraph::SubGraph, HedgeGraph};

use typed_index_collections::TiVec;

use super::{
    approx::Approximation,
    poset::{DagNode, DAG},
    uv_graph::UVE,
    UltravioletGraph,
};

pub struct Forest {
    pub dag: DAG<Approximation, DagNode, ()>,
}

impl Forest {
    pub fn n_terms(&self) -> usize {
        self.dag.nodes.len()
    }

    pub fn compute<
        E: UVE,
        V,
        H,
        G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>,
        S: SubGraph,
        OID: OrientationID,
        O: GraphOrientation,
    >(
        &mut self,
        graph: &G,
        amplitude_subgraph: &S,
        orientations: &TiVec<OID, O>,
        canonize_esurface: &Option<ShiftRewrite>,
        cut_edges: &[EdgeIndex],
    ) {
        let order = self.dag.compute_topological_order();

        for n in order {
            match self.dag.nodes[n].parents.len() {
                0 => {
                    self.dag.nodes[n].data.root(
                        graph.as_ref(),
                        amplitude_subgraph,
                        canonize_esurface,
                        orientations,
                        cut_edges,
                    );
                }
                1 => {
                    let [ref mut current, ref mut parent] = &mut self
                        .dag
                        .nodes
                        .get_disjoint_mut([n, self.dag.nodes[n].parents[0]])
                        .unwrap();

                    assert!(matches!(parent.data.local_3d, CFFapprox::Dependent { .. }));
                    current
                        .data
                        .compute(graph, amplitude_subgraph, &parent.data);
                    current.data.compute_3d(
                        graph,
                        amplitude_subgraph,
                        canonize_esurface,
                        orientations,
                        cut_edges,
                        &parent.data,
                    );
                }
                _ => {
                    unimplemented!("Union not implemented");
                    // let (local_4d, simple, local_3d) = {
                    //     let current = &self.dag.nodes[n];
                    //     let mut dependents =
                    //         current.parents.iter().map(|p| &self.dag.nodes[*p].data);

                    //     let mut dep_vec: Vec<&Approximation> = Vec::new();
                    //     dep_vec.push(dependents.next().unwrap());

                    //     let Some((mut approx, _)) = dep_vec[0].local_3d.expr() else {
                    //         panic!("local 3d not computed for parents")
                    //     };

                    //     let mut iterative_approx = dep_vec[0].clone();

                    //     for p in dependents {
                    //         dep_vec.push(p);

                    //         approx = p.local_3d(&iterative_approx, graph, approx);

                    //         iterative_approx = p.clone();
                    //     }
                    //     (
                    //         ApproxOp::union(&dep_vec).unwrap(),
                    //         SimpleApprox::union(
                    //             self.dag.nodes[n].data.subgraph.clone(),
                    //             dep_vec.iter().map(|s| s.simple_approx.as_ref().unwrap()),
                    //         ),
                    //         iterative_approx.local_3d,
                    //     )
                    // };

                    // self.dag.nodes[n].data.local_4d = local_4d;
                    // self.dag.nodes[n].data.simple_approx = Some(simple);
                    // self.dag.nodes[n].data.local_3d = local_3d;
                }
            }
        }
    }

    pub fn local_expr<G, H, S: SubGraph<Base = BitVec>, OID: OrientationID, O: GraphOrientation>(
        &self,
        graph: &G,
        amplitude: &S,
        orientation: &OrientationData,
        canonize_esurface: &Option<ShiftRewrite>,
        orientations: &TiVec<OID, O>,
        cut_edges: &[EdgeIndex],
    ) -> Atom
    where
        G: UltravioletGraph + AsRef<HedgeGraph<Edge, Vertex, H>>,
    {
        let mut sum = Atom::Zero;

        for (_, n) in &self.dag.nodes {
            sum += n
                .data
                .final_cff(
                    graph,
                    amplitude,
                    orientation,
                    canonize_esurface,
                    orientations,
                    cut_edges,
                )
                .unwrap();
        }

        sum
    }

    // pub fn simple_expr(
    //     &self,
    //     graph: &UVGraph,
    //     amplitude: &InternalSubGraph,
    // ) -> Option<SerializableAtom> {
    //     let mut sum = Atom::num(0).into();
    //     for (_, n) in &self.dag.nodes {
    //         sum = sum + n.data.simple_expr(graph, amplitude)?;
    //     }

    //     Some(sum)
    // }

    pub fn graphs(&self) -> String {
        self.dag
            .to_dot_impl(&|n| format!("label=S_{}", n.data.subgraph.string_label()))
    }

    // pub fn show_structure(&self, graph: &UVGraph, amplitude: &InternalSubGraph) -> Option<String> {
    //     let mut out = String::new();

    //     out.push_str(
    //         &self
    //             .simple_expr(graph, amplitude)?
    //             .0
    //             .printer(PrintOptions {
    //                 terms_on_new_line: true,
    //                 ..Default::default()
    //             })
    //             .to_string(),
    //     );
    //     Some(out)
    // }
}
