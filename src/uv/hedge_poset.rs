use std::fmt::Display;

use ahash::AHashSet;
use linnet::half_edge::{
    HedgeGraph, NoData, NodeIndex,
    algorithms::trace_unfold::{HiddenData, Independence, TraceKey, TraceUnfold},
    involution::{EdgeIndex, Flow, HedgePair},
    nodestore::{NodeStorageOps, NodeStorageVec},
    subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps},
};

#[derive(Clone, Debug, Hash, Eq)]
pub struct Spinney {
    pub subgraph: SuBitGraph,
    // pub dod: i32,
    // pub lmb: LoopMomentumBasis,
}

impl PartialEq for Spinney {
    fn eq(&self, other: &Self) -> bool {
        self.subgraph == other.subgraph
    }
}

impl PartialOrd for Spinney {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.subgraph.includes(&other.subgraph) {
            Some(std::cmp::Ordering::Greater)
        } else if other.subgraph.includes(&self.subgraph) {
            Some(std::cmp::Ordering::Less)
        } else {
            None
        }
    }
}

pub struct SpinneyWood {
    pub graph: HedgeGraph<SuBitGraph, Spinney>,
    pub root: NodeIndex,
}

impl Independence<HiddenData<SuBitGraph, EdgeIndex>> for SpinneyWood {
    fn independent(
        &self,
        a: &HiddenData<SuBitGraph, EdgeIndex>,
        b: &HiddenData<SuBitGraph, EdgeIndex>,
    ) -> bool {
        !a.order.intersects(&b.order)
    }
}

impl TraceUnfold<SuBitGraph> for SpinneyWood {
    type EdgeData = SuBitGraph;
    type HedgeData = NoData;
    type NodeData = Spinney;
    type NodeStorage = NodeStorageVec<Spinney>;

    fn graph(
        &self,
    ) -> &HedgeGraph<Self::EdgeData, Self::NodeData, Self::HedgeData, Self::NodeStorage> {
        &self.graph
    }

    fn key(&self, e: EdgeIndex) -> SuBitGraph {
        self.graph[e].clone()
    }
}

impl SpinneyWood {
    pub(crate) fn from_spinneys<E, V, H, I: IntoIterator<Item = Spinney>>(
        s: I,
        graph: impl AsRef<HedgeGraph<E, V, H>>,
    ) -> Self {
        let mut spinneyset: AHashSet<_> = s.into_iter().collect();
        spinneyset.insert(Spinney {
            subgraph: graph.as_ref().empty_subgraph(),
        });

        let mut unions = AHashSet::new();
        let g: HedgeGraph<_, _> = HedgeGraph::poset(spinneyset);
        let mut poset = g.map(
            |_, _, v| v,
            |_, n, pair, _, e| {
                let HedgePair::Paired { source, sink } = pair else {
                    return e.map(|_| graph.as_ref().empty_subgraph());
                };
                let nsource = n.node_id_ref(source);
                let nsink = n.node_id_ref(sink);

                let source_subgraph = &n.get_node_data(nsource).subgraph;
                let sink_subgraph = &n.get_node_data(nsink).subgraph;
                let reduced_subgraph = sink_subgraph.subtract(source_subgraph);

                // println!("Reduced: {}", graph.as_ref().dot(&reduced_subgraph));
                // println!(
                //     "Bridges: {}",
                //     graph
                //         .as_ref()
                //         .dot(&graph.as_ref().bridges_of(&reduced_subgraph))
                // );
                //
                //
                let hairy_source = graph.as_ref().full_crown(source_subgraph);

                if graph.as_ref().bridges_of(&reduced_subgraph).is_empty()
                    && !hairy_source.intersects(&reduced_subgraph)
                {
                    // if the reduced graph is bridgless,and has no node overlap with the source, then the source subgraph is cycle independent of reduced subgraph

                    // if the source is not empty then this is a disjoint union
                    if !source_subgraph.is_empty() {
                        unions.insert(nsink);

                        // println!("//{nsink}:{}", reduced_subgraph.string_label());
                    }
                    return e.map(|_| reduced_subgraph);
                } else {
                    return e.map(|_| sink_subgraph.clone());
                }
            },
            |_, d| d,
        );

        let mut to_remove: SuBitGraph = poset.empty_subgraph();

        //Not quite transitive closure. We need to remove all edges that point to unions that are not single components
        for u in unions {
            // println!("//{u}:{}", poset[u].subgraph.string_label());
            let mut comps: AHashSet<_> = graph
                .as_ref()
                .connected_components(&poset[u].subgraph)
                .into_iter()
                .collect();
            for c in poset.iter_crown(u) {
                let Flow::Sink = poset.flow(c) else {
                    continue;
                };

                let Some(nid) = poset.involved_node_id(c) else {
                    continue;
                };

                // if poset[nid].subgraph.is_empty() {
                //     continue;
                // }
                if comps.contains(&poset[nid].subgraph) {
                    comps.remove(&poset[nid].subgraph);
                } else {
                    to_remove.add(c);
                    to_remove.add(poset.inv(c));
                }
            }
        }

        println!("Removing: {}", poset.dot(&to_remove));

        poset.delete_hedges(&to_remove);
        let root = poset
            .iter_nodes()
            .find(|(_, _, s)| s.subgraph.is_empty())
            .map(|(n, _, _)| n);

        SpinneyWood {
            graph: poset,
            root: root.expect("no empty spinney found"),
        }
    }

    pub fn unfold(&self) -> SpinneyForest {
        SpinneyForest {
            graph: self.trace_unfold::<NodeStorageVec<_>>(self.root),
            root: self.root.clone(),
        }
    }
}

impl Display for SpinneyWood {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.graph.dot_impl_fmt(
            f,
            &self.graph.full_filter(),
            "start=2;\n",
            &|_| None,
            &|a| Some(format!("label=\"{}\"", a.string_label())),
            &|v| Some(format!("label=\"{}\"", v.subgraph.string_label())),
        )
    }
}
// pub type SpinneyGraph = HedgeGraph<SuBitGraph, Spinney>;

// impl Key<&SuBitGraph> for SpinneyGraph {
//     fn key(&self, e: linnet::half_edge::involution::EdgeIndex) -> K {
//         &self[e]
//     }
// }

pub struct SpinneyForest {
    pub graph: HedgeGraph<NoData, TraceKey<SuBitGraph, EdgeIndex>>,
    pub root: NodeIndex,
}

impl SpinneyForest {
    pub fn walk(&self) {
        self.graph
            .topo_sort_kahn()
            .unwrap()
            .iter()
            .for_each(|nidx| {
                let trace_key = &self.graph[*nidx];
                println!("Node {}: {}", nidx, trace_key_label(trace_key));
            });
    }
}

pub fn trace_key_label(v: &TraceKey<SuBitGraph, EdgeIndex>) -> String {
    if v.levels.is_empty() {
        "∅".to_string()
    } else {
        let mut s = String::new();
        let mut acc: Option<SuBitGraph> = None;

        for (i, level) in v.levels.iter().enumerate() {
            if i > 0 {
                s.push_str(" · ");
            }
            s.push_str("{");
            for (j, op) in level.iter().enumerate() {
                if j > 0 {
                    s.push_str(",");
                }
                if let Some(a) = &mut acc {
                    a.union_with(&op.order);
                } else {
                    acc = Some(op.order.clone());
                }
                s.push_str(&op.order.string_label());
            }
            s.push_str("}");
        }
        if let Some(a) = acc {
            s.push_str(" => ");
            s.push_str(&a.string_label());
        }

        s
    }
}

impl Display for SpinneyForest {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.graph.dot_impl_fmt(
            f,
            &self.graph.full_filter(),
            "start=2;\n",
            &|_| None,
            &|_| None,
            &|v| Some(format!("label=\"{}\"", trace_key_label(v))),
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        dot,
        graph::{Graph, parse::IntoGraph},
        initialisation::test_initialise,
        uv::{UltravioletGraph, Wood},
    };

    use super::*;
    use color_eyre::Result;

    #[test]
    fn dumbells() -> Result<()> {
        test_initialise().unwrap();
        let dumbell: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v2 -> v2;
                v1 -> v1;
            },"scalars"
        )?;

        let spinneys = dumbell.spinneys(&dumbell.full_filter());
        let f = SpinneyWood::from_spinneys(
            dumbell
                .spinneys(&dumbell.full_filter())
                .into_iter()
                .map(|a| Spinney { subgraph: a.filter }),
            &dumbell,
        );

        println!("{}", f);

        insta::assert_snapshot!(
            f.graph.n_nodes(),
            @"4",
            // format!("Wood does not have correct number of spinneys: \n{}",f)
        );

        for (_, _, d) in f.graph.iter_nodes() {
            println!(
                "//Node {}: \n{}",
                d.subgraph.string_label(),
                dumbell.dot(&d.subgraph)
            );
        }
        let ff = Wood::from_spinneys(spinneys, &dumbell); //.unfold(&g, &g.loop_momentum_basis);

        println!("{}", ff.dot(&dumbell));

        let f = f.unfold();
        insta::assert_snapshot!(
            f.graph.n_nodes(),
            @"4");

        Ok(())
    }

    #[test]
    fn bugblatter() -> Result<()> {
        test_initialise().unwrap();

        match dot!(
            digraph G{
                A1 -> A2 [particle="t"];
                A2 -> A3 [particle="t"];
                A3 -> A1 [particle="t"];
                B1 -> B2 [particle="t"];
                B2 -> B3 [particle="t"];
                B3 -> B1 [particle="t"];
                A1 -> B1 [particle="a"];
                A2 -> B2 [particle="a"];
                A3 -> B3 [particle="a"];
            },"sm"
        ) {
            Ok(g) => {
                let g: Graph = g;
                let spinneys = g.spinneys(&g.full_filter());
                let f = SpinneyWood::from_spinneys(
                    g.spinneys(&g.full_filter())
                        .into_iter()
                        .map(|a| Spinney { subgraph: a.filter }),
                    &g,
                );

                println!("{}", f);

                assert_eq!(
                    20,
                    f.graph.n_nodes(),
                    "Wood does not have correct number of spinneys: \n{}",
                    f
                );

                for (_, _, d) in f.graph.iter_nodes() {
                    println!(
                        "//Node {}: \n{}",
                        d.subgraph.string_label(),
                        g.dot(&d.subgraph)
                    );
                }
                let ff = Wood::from_spinneys(spinneys, &g); //.unfold(&g, &g.loop_momentum_basis);

                println!("{}", ff.dot(&g));

                let f = f.unfold();
                f.walk();
                assert_eq!(
                    152,
                    f.graph.n_nodes(),
                    "Forest unfolds into the wrong number of terms :\n{}",
                    f
                );

                // println!("{}", f)
            }
            Err(e) => {
                eprintln!("{}", e);
            }
        }

        // let f = SpinneyWood::from_spinneys(
        //     g.spinneys(&g.full_filter()).into_iter().map(|a| a.filter),
        //     &g,
        // )
        // .unfold();

        // println!("{}", f.graph.base_dot());
        Ok(())
    }

    #[test]
    fn mercedes() -> Result<()> {
        test_initialise().unwrap();

        let mercedes: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v2 -> v3;
                v3 -> v1;
                v1 -> v4;
                v2 -> v4;
                v3 -> v4;
            },"scalars"
        )?;

        let f = SpinneyWood::from_spinneys(
            mercedes
                .spinneys(&mercedes.full_filter())
                .into_iter()
                .map(|a| Spinney { subgraph: a.filter }),
            &mercedes,
        );
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"2",
        );
        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"2",
         );

        Ok(())
    }

    #[test]
    fn spectacles() -> Result<()> {
        test_initialise().unwrap();

        let spectacles: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v1 -> v2;

                v3 -> v4;
                v3 -> v4;

                v2 -> v3;
                v1 -> v4;
            },"scalars"
        )?;

        // let spinneys = spectacles.spinneys(&spectacles.full_filter());
        let f = SpinneyWood::from_spinneys(
            spectacles
                .spinneys(&spectacles.full_filter())
                .into_iter()
                .map(|a| Spinney { subgraph: a.filter }),
            &spectacles,
        );
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"5",
        );
        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"8",
         );

        Ok(())
    }

    #[test]
    fn basketball() -> Result<()> {
        test_initialise().unwrap();

        let basketball: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v1 -> v2;
                v1 -> v2;
                v1 -> v2;
            },"scalars"
        )?;

        let f = SpinneyWood::from_spinneys(
            basketball
                .spinneys(&basketball.full_filter())
                .into_iter()
                .map(|a| Spinney { subgraph: a.filter }),
            &basketball,
        );
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"12",
              );
        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"46");
        Ok(())
    }

    #[test]
    fn fourloop_b() -> Result<()> {
        test_initialise().unwrap();

        let fourloop_b: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v1 -> v2;

                v3 -> v4;
                v3 -> v4;

                v2 -> v3;
                v1 -> v3;
                v1 -> v4;
            },"scalars"
        )?;

        // let spinneys = fourloop_b.spinneys(&fourloop_b.full_filter());
        let f = SpinneyWood::from_spinneys(
            fourloop_b
                .spinneys(&fourloop_b.full_filter())
                .into_iter()
                .map(|a| Spinney { subgraph: a.filter }),
            &fourloop_b,
        );
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"14",
        );
        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"80",
         );

        Ok(())
    }

    #[test]
    fn four_loop_a() -> Result<()> {
        test_initialise().unwrap();

        let four_loop_a: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v1 -> v2;
                v2 -> v3;
                v3 -> v1;
                v1 -> v4;
                v2 -> v4;
                v3 -> v4;
            },"scalars"
        )?;

        let f = SpinneyWood::from_spinneys(
            four_loop_a
                .spinneys(&four_loop_a.full_filter())
                .into_iter()
                .map(|a| Spinney { subgraph: a.filter }),
            &four_loop_a,
        );
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"12",
        );
        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"60",
         );

        Ok(())
    }
}
