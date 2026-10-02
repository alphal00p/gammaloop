use crate::half_edge::subgraph::{ModifySubSet, SuBitGraph, SubSetLike};
use crate::half_edge::swap::Swap;

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

#[test]
fn involution_iterators_keep_exact_remaining_cardinality() {
    fn check<I: ExactSizeIterator<Item = (super::Hedge, T)> + std::iter::FusedIterator, T>(
        mut iter: I,
        count: usize,
    ) {
        for position in 0..count {
            assert_eq!(iter.len(), count - position);
            assert_eq!(iter.size_hint(), (count - position, Some(count - position)));
            assert_eq!(iter.next().unwrap().0, super::Hedge(position));
        }
        for _ in 0..3 {
            assert_eq!(iter.len(), 0);
            assert_eq!(iter.size_hint(), (0, Some(0)));
            assert!(iter.next().is_none());
        }
    }
    for count in 0..32 {
        let mut source = Involution::new();
        for id in 0..count {
            if id % 3 == 0 {
                source.add_pair(id, id % 2 == 0);
            } else {
                source.add_identity(id, id % 2 == 0, Flow::Source);
            }
        }
        let hedges: super::Hedge = source.len();
        check((&source).into_iter(), hedges.0);
        check((&mut source).into_iter(), hedges.0);
        check(source.into_iter(), hedges.0);
    }
}

#[test]
fn draining_involution_keeps_unconsumed_nonclone_drop_order() {
    use std::{cell::RefCell, rc::Rc};
    struct Payload(usize, Rc<RefCell<Vec<usize>>>);
    impl Drop for Payload {
        fn drop(&mut self) {
            self.1.borrow_mut().push(self.0);
        }
    }
    for consumed in 0..=8 {
        let drops = Rc::new(RefCell::new(Vec::new()));
        let mut source = Involution::new();
        for id in 0..8 {
            source.add_identity(Payload(id, drops.clone()), true, Flow::Source);
        }
        let mut iter = source.into_iter();
        assert!(drops.borrow().is_empty());
        for id in 0..consumed {
            drop(iter.next().unwrap());
            assert_eq!(&*drops.borrow(), &(0..=id).collect::<Vec<_>>());
            assert_eq!(iter.size_hint(), (7 - id, Some(7 - id)));
        }
        drop(iter);
        assert_eq!(&*drops.borrow(), &(0..8).collect::<Vec<_>>());
    }
}
