use crate::half_edge::subgraph::{ModifySubSet, SuBitGraph, SubSetLike};

use super::{Flow, Involution};

#[test]
fn invextract() {
    let mut a = Involution::new();

    let _s1 = a.add_identity(0, false, Flow::Sink);
    let _s2 = a.add_identity(1, false, Flow::Sink);
    let _s3 = a.add_identity(2, false, Flow::Sink);
    let h1 = a.add_identity(3, false, Flow::Sink);
    let s4 = a.add_identity(4, false, Flow::Sink);
    let h3 = a.add_identity(5, false, Flow::Sink);
    let h2 = a.add_identity(6, false, Flow::Sink);
    let h4 = a.add_identity(7, false, Flow::Sink);

    a.connect_identities(h1, h2, |f, d, _, _| (f, d)).unwrap();

    a.connect_identities(h3, h4, |f, d, _, _| (f, d)).unwrap();

    let mut subgraph = SuBitGraph::empty(8);
    subgraph.add(s4);
    subgraph.add(h2);
    subgraph.add(h3);
    subgraph.add(h4);

    println!("{a}");
    let extracted = a.extract(&subgraph, |a| a.map(Clone::clone), |a| a);

    println!("{a}");
    println!("{extracted}");
}
