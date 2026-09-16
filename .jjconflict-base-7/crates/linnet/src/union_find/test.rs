use super::*;
// use crate::union_find::bitvec_find::*;

// Helper functions for testing
fn sum_merge(a: i32, b: i32) -> i32 {
    a + b
}
// fn concat_merge(a: String, b: String) -> String {
//     format!("{}{}", a, b)
// }

#[test]
fn test_basic_union_find() {
    let data = vec![10, 20, 30, 40, 50];
    let mut uf = UnionFind::new(data);

    // Test initial state
    assert_eq!(*uf.find_data(Hedge(0)), 10);
    assert_eq!(*uf.find_data(Hedge(1)), 20);

    // Test union
    let _ = uf.union(Hedge(0), Hedge(1), sum_merge);
    assert_eq!(*uf.find_data(Hedge(0)), 30); // 10 + 20
    assert_eq!(uf.find(Hedge(0)), uf.find(Hedge(1)));
}

#[test]
fn test_path_compression() {
    let data = vec![10, 20, 30, 40];
    let mut uf = UnionFind::new(data);

    // Create a chain: 3->2->1->0
    uf.union(Hedge(0), Hedge(1), sum_merge);
    uf.union(Hedge(1), Hedge(2), sum_merge);
    uf.union(Hedge(2), Hedge(3), sum_merge);

    // Find should compress the path
    let root = uf.find(Hedge(3));
    assert_eq!(uf.find(Hedge(2)), root);
    assert_eq!(uf.find(Hedge(1)), root);
}
