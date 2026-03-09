use std::cmp::Ordering;
use std::collections::VecDeque;
use std::fmt::{self, Display};
use std::hash::{Hash, Hasher};

use indexmap::set::MutableValues;
use indexmap::IndexSet;

use crate::half_edge::builder::HedgeGraphBuilder;
use crate::half_edge::involution::{EdgeIndex, Flow};
use crate::half_edge::nodestore::NodeStorageOps;
use crate::half_edge::{HedgeGraph, NoData, NodeIndex};

/// Ops must be owned, hashable, and totally ordered for canonicalization.
pub trait Op: Clone + Eq + Ord {}
impl<T: Clone + Eq + Ord> Op for T {}

pub trait Independence<O> {
    fn independent(&self, a: &O, b: &O) -> bool;
}

pub trait Key<K> {
    fn key(&self, e: EdgeIndex) -> K;
}

pub trait TraceUnfold<Key>: Independence<HiddenData<Key, EdgeIndex>> + Sized
where
    Key: Eq + Hash + Clone + Ord,
{
    type EdgeData: Ord + Hash;
    type NodeData;
    type NodeStorage: NodeStorageOps<NodeData = Self::NodeData>;
    type HedgeData;

    fn graph(
        &self,
    ) -> &HedgeGraph<Self::EdgeData, Self::NodeData, Self::HedgeData, Self::NodeStorage>;

    fn key(&self, e: EdgeIndex) -> Key;

    fn trace_unfold<M: NodeStorageOps<NodeData = usize>>(
        &self,
        start: NodeIndex,
    ) -> HedgeGraph<NoData, TraceKey<Key, EdgeIndex>, NoData, M::OpStorage<TraceKey<Key, EdgeIndex>>>
    {
        let root = (start, TraceKey::empty());
        let mut q = VecDeque::new();
        let mut traces: IndexSet<(NodeIndex, TraceKey<Key, EdgeIndex>)> = IndexSet::new();
        let mut builder: HedgeGraphBuilder<NoData, usize, NoData> = HedgeGraphBuilder::new();
        let (ind, _) = traces.insert_full(root);
        let nid = builder.add_node(ind);
        q.push_back((nid, start, TraceKey::empty()));
        let g = self.graph();

        while let Some((bnid, nid, key)) = q.pop_front() {
            for hedge in g.iter_crown(nid) {
                if g.flow(hedge) == Flow::Source {
                    if let Some(to_node) = g.involved_node_id(hedge) {
                        let edge = self.key(g[&hedge]);

                        let new_key: TraceKey<Key, EdgeIndex> = key.push(
                            self,
                            HiddenData {
                                order: edge,
                                data: g[&hedge],
                            },
                        );

                        let (ind, is_new) = traces.insert_full((to_node, new_key.clone()));
                        if is_new {
                            let bbnid = builder.add_node(ind);
                            debug_assert_eq!(bbnid.0, ind);
                            q.push_back((bbnid, to_node, new_key.clone()))
                        }

                        builder.add_edge(bnid, NodeIndex(ind), NoData {}, true);
                    }
                }
            }
        }

        builder.build::<M>().map(
            |_, _, v| {
                let mut trace: TraceKey<Key, EdgeIndex> = TraceKey::empty();
                let v = &mut traces.get_index_mut2(v).unwrap().1;
                std::mem::swap(v, &mut trace);
                trace
            },
            |_, _, _, _, a| a,
            |_, h| h,
        )
    }
}

impl<E: Ord + Hash, H, V, N: NodeStorageOps<NodeData = V>, K: Eq + Hash + Clone + Ord>
    TraceUnfold<K> for HedgeGraph<E, V, H, N>
where
    Self: Independence<HiddenData<K, EdgeIndex>>,
    Self: Key<K>,
{
    type EdgeData = E;
    type NodeData = V;
    type NodeStorage = N;
    type HedgeData = H;

    fn graph(
        &self,
    ) -> &HedgeGraph<Self::EdgeData, Self::NodeData, Self::HedgeData, Self::NodeStorage> {
        self
    }

    fn key(&self, e: EdgeIndex) -> K {
        Key::key(self, e)
    }
}

// impl<E, H, V> Independence<HedgePair> for HedgeGraph<E, V, H>
// where
//     Self: Independence<E> + Independence<V> + Independence<H>,
// {
//     fn independent(&self, a: &HedgePair, b: &HedgePair) -> bool {
//         match (a, b) {
//             (
//                 HedgePair::Paired { source, sink },
//                 HedgePair::Paired {
//                     source: so,
//                     sink: sk,
//                 },
//             ) => {

//                 let a = self[source];
//                 let v = self.node_id(*source);
//                 self.independent(&self[source], &self[so])
//                     && self.independent(&self[sink], &self[sk]) && self.independent(&self[self[]], b)
//             }
//             _ => {}
//         }

//         false
//     }
// }

#[derive(Clone, Debug)]
pub struct HiddenData<O, D> {
    pub order: O, // the semantic label used for identity + ordering
    pub data: D,  // extra context used only for independence tests
}

// Equality/hash only by `order` so merging works.
impl<O: PartialEq, D> PartialEq for HiddenData<O, D> {
    fn eq(&self, other: &Self) -> bool {
        self.order == other.order
    }
}
impl<O: Eq, D> Eq for HiddenData<O, D> {}

impl<O: Hash, D> Hash for HiddenData<O, D> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.order.hash(state);
    }
}

// Sorting only by `order` so canonicalization works.
impl<O: Ord, D> Ord for HiddenData<O, D> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.order.cmp(&other.order)
    }
}
impl<O: PartialOrd, D> PartialOrd for HiddenData<O, D> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.order.partial_cmp(&other.order)
    }
}

#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct TraceKey<O, D> {
    pub levels: Vec<Vec<HiddenData<O, D>>>,
}

impl<O, D> TraceKey<O, D> {
    pub fn empty() -> Self {
        Self { levels: Vec::new() }
    }

    pub fn push<I>(&self, indep: &I, op: HiddenData<O, D>) -> Self
    where
        I: Independence<HiddenData<O, D>>,
        O: Clone + Eq + Ord,
        D: Clone,
    {
        let mut levels = self.levels.clone();

        for i in (0..levels.len()).rev() {
            let ok = levels[i].iter().all(|b| indep.independent(&op, b));
            if ok {
                levels[i].push(op);
                levels[i].sort(); // uses O only
                return Self { levels };
            }
        }

        levels.push(vec![op]);
        Self { levels }
    }

    pub fn write_foata_like<W: fmt::Write>(
        &self,
        f: &mut W,
        mut map: impl FnMut(&O) -> String,
    ) -> fmt::Result
    where
        O: Op,
    {
        if self.levels.is_empty() {
            return write!(f, "∅");
        }

        for (i, level) in self.levels.iter().enumerate() {
            if i > 0 {
                write!(f, " · ")?;
            }
            write!(f, "{{")?;
            for (j, op) in level.iter().enumerate() {
                if j > 0 {
                    write!(f, ",")?;
                }
                write!(f, "{}", map(&op.order))?;
            }
            write!(f, "}}")?;
        }
        Ok(())
    }
}

/// Display as Foata-like levels: {A,C} · {B} · {D,E}; empty = ∅
impl<O, D> Display for TraceKey<O, D>
where
    O: Op + Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.write_foata_like(f, |o| o.to_string())
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use crate::{
        dot,
        half_edge::{involution::HedgePair, nodestore::DefaultNodeStore},
        parser::{DotEdgeData, DotGraph, DotHedgeData, DotVertexData},
    };

    use super::*;
    use ahash::HashSet;
    use insta::assert_snapshot;

    #[derive(Clone, Eq, PartialEq, Hash, Ord, PartialOrd, Debug)]
    enum MyOp {
        A,
        B,
        C,
    }

    impl Display for MyOp {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            let s = match self {
                MyOp::A => "A",
                MyOp::B => "B",
                MyOp::C => "C",
            };
            write!(f, "{}", s)
        }
    }

    // Independence where A ⟂ C only; B overlaps with both.
    struct IndepACOnly;
    impl Independence<MyOp> for IndepACOnly {
        fn independent(&self, a: &MyOp, b: &MyOp) -> bool {
            matches!((a, b), (MyOp::A, MyOp::C) | (MyOp::C, MyOp::A))
        }
    }

    struct OverlapIndep;
    impl Independence<String> for OverlapIndep {
        fn independent(&self, a: &String, b: &String) -> bool {
            !has_overlap(a, b)
        }
    }

    fn has_overlap(a: &str, b: &str) -> bool {
        println!("Checking overlap between '{a}' and '{b}'");
        let set_a: BTreeSet<char> = a.chars().collect();
        b.chars().any(|c| set_a.contains(&c))
    }

    fn remove_chars(source: &str, remove: &str) -> String {
        let to_remove: HashSet<char> = remove.chars().collect();
        source.chars().filter(|c| !to_remove.contains(c)).collect()
    }

    impl Key<String> for HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData> {
        fn key(&self, e: EdgeIndex) -> String {
            let HedgePair::Paired { source, sink } = self[&e].1 else {
                return String::new();
            };

            let source_name = self[self.node_id(source)].name.clone().unwrap_or_default();
            let sink_name = self[self.node_id(sink)].name.clone().unwrap_or_default();

            // Remove "empty" from the name for key purposes
            remove_chars(&sink_name, &source_name)
        }
    }

    impl Independence<HiddenData<String, EdgeIndex>>
        for HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData>
    {
        fn independent(
            &self,
            a: &HiddenData<String, EdgeIndex>,
            b: &HiddenData<String, EdgeIndex>,
        ) -> bool {
            let HedgePair::Paired {
                source: sourcea, ..
            } = self[&a.data].1
            else {
                return true;
            };

            let HedgePair::Paired {
                source: sourceb, ..
            } = self[&b.data].1
            else {
                return true;
            };

            let sourcan = self[self.node_id(sourcea)].name.clone().unwrap_or_default();
            let sourcbn = self[self.node_id(sourceb)].name.clone().unwrap_or_default();

            sourcan.as_str() == "empty"
                || sourcbn.as_str() == "empty"
                || !has_overlap(&sourcan, &sourcbn)
        }
    }

    // #[test]
    // fn overlappping_trace() {
    //     let indep = OverlapIndep;
    //     let t_ac = TraceKey::empty()
    //         .push(&indep, "A".to_string())
    //         .push(&indep, "B".to_string());
    //     let t_ca = TraceKey::empty()
    //         .push(&indep, "B".to_string())
    //         .push(&indep, "A".to_string());

    //     // Same trace key => same display
    //     assert_snapshot!(t_ac.to_string(), @"{A,B}");
    //     assert_snapshot!(t_ca.to_string(), @"{A,B}");
    // }
    #[test]
    fn test_tbt() {
        let graph: DotGraph = dot!(digraph{
            empty [id=0 name=""]
            empty -> AB [id=0 ];
            empty -> CD [id=1];
            AB -> ABCD [id=2]
            AB -> ABCEF [id=3 ]
            CD -> BCDEF [id=4 ]
            CD -> ABCD [id=5]
            ABCD -> ABCDEF [id=6]
            BCDEF -> ABCDEF [id=7]
            ABCEF -> ABCDEF [id=8]
        })
        .unwrap();

        let g: HedgeGraph<NoData, TraceKey<String, EdgeIndex>> = graph
            .graph
            .transitive_closure()
            .unwrap()
            .trace_unfold::<DefaultNodeStore<usize>>(NodeIndex(0));
        let mut output = String::new();
        g.dot_impl_fmt(
            &mut output,
            &g.full_filter(),
            "start=2;\n",
            &|_| None,
            &|a| None,
            &|v| Some(format!("label=\"{v}\"")),
        )
        .unwrap();

        insta::assert_snapshot!(output, @r#"
        digraph {
          node	 [shape=circle,height=0.1,label=""];
          overlap = "scale";
          layout = "neato";
          start=2;
          0	 [label="∅"];
          1	 [label="{AB}"];
          2	 [label="{CD}"];
          3	 [label="{ABCEF}"];
          4	 [label="{ABCD}"];
          5	 [label="{BCDEF}"];
          6	 [label="{ABCDEF}"];
          7	 [label="{AB,CD}"];
          8	 [label="{AB,CEF}"];
          9	 [label="{AB,CDEF}"];
          10	 [label="{BEF,CD}"];
          11	 [label="{ABEF,CD}"];
          12	 [label="{ABCEF,D}"];
          13	 [label="{ABCD,EF}"];
          14	 [label="{A,BCDEF}"];
          15	 [label="{AB,CD} · {EF}"];
          16	 [label="{AB,CEF} · {D}"];
          17	 [label="{BEF,CD} · {A}"];
          0:0:s	-> 1:1:s	 [id=0  color="red:blue;0.5"];
          0:2:s	-> 2:3:s	 [id=1  color="red:blue;0.5"];
          0:4:s	-> 3:5:s	 [id=2  color="red:blue;0.5"];
          0:6:s	-> 4:7:s	 [id=3  color="red:blue;0.5"];
          0:8:s	-> 5:9:s	 [id=4  color="red:blue;0.5"];
          0:10:s	-> 6:11:s	 [id=5  color="red:blue;0.5"];
          1:12:s	-> 7:13:s	 [id=6  color="red:blue;0.5"];
          1:14:s	-> 8:15:s	 [id=7  color="red:blue;0.5"];
          1:16:s	-> 9:17:s	 [id=8  color="red:blue;0.5"];
          2:18:s	-> 7:19:s	 [id=9  color="red:blue;0.5"];
          2:20:s	-> 10:21:s	 [id=10  color="red:blue;0.5"];
          2:22:s	-> 11:23:s	 [id=11  color="red:blue;0.5"];
          3:24:s	-> 12:25:s	 [id=12  color="red:blue;0.5"];
          4:26:s	-> 13:27:s	 [id=13  color="red:blue;0.5"];
          5:28:s	-> 14:29:s	 [id=14  color="red:blue;0.5"];
          7:30:s	-> 15:31:s	 [id=15  color="red:blue;0.5"];
          8:32:s	-> 16:33:s	 [id=16  color="red:blue;0.5"];
          10:34:s	-> 17:35:s	 [id=17  color="red:blue;0.5"];
        }
        "#);
    }

    // #[test]
    // fn snapshot_empty() {
    //     let t = TraceKey::<MyOp>::empty();
    //     assert_snapshot!(t.to_string(), @"∅");
    // }

    // #[test]
    // fn snapshot_commutation_merges_ac_and_ca() {
    //     let indep = IndepACOnly;
    //     let t_ac = TraceKey::empty()
    //         .push(&indep, MyOp::A)
    //         .push(&indep, MyOp::C);
    //     let t_ca = TraceKey::empty()
    //         .push(&indep, MyOp::C)
    //         .push(&indep, MyOp::A);

    //     // Same trace key => same display
    //     assert_snapshot!(t_ac.to_string(), @"{A,C}");
    //     assert_snapshot!(t_ca.to_string(), @"{A,C}");
    // }

    // #[test]
    // fn snapshot_dependent_orders_differ() {
    //     let indep = IndepACOnly;

    //     let t_ab = TraceKey::empty()
    //         .push(&indep, MyOp::A)
    //         .push(&indep, MyOp::B);
    //     let t_ba = TraceKey::empty()
    //         .push(&indep, MyOp::B)
    //         .push(&indep, MyOp::A);

    //     assert_snapshot!(t_ab.to_string(), @"{A} · {B}");
    //     assert_snapshot!(t_ba.to_string(), @"{B} · {A}");
    // }

    // #[test]
    // fn snapshot_ac_then_b_creates_new_level() {
    //     let indep = IndepACOnly;

    //     let t = TraceKey::empty()
    //         .push(&indep, MyOp::A)
    //         .push(&indep, MyOp::C)
    //         .push(&indep, MyOp::B);

    //     assert_snapshot!(t.to_string(), @"{A,C} · {B}");
    // }

    // #[test]
    // fn snapshot_place_as_late_as_possible() {
    //     let indep = IndepACOnly;

    //     // Start with {B} · {C}; then push A.
    //     // A can't join {B} but can join {C} => {B} · {A,C}
    //     let t = TraceKey::empty()
    //         .push(&indep, MyOp::B)
    //         .push(&indep, MyOp::C)
    //         .push(&indep, MyOp::A);

    //     assert_snapshot!(t.to_string(), @"{B} · {A,C}");
    // }
}
