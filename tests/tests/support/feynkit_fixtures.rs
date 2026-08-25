#![allow(dead_code)]

use gammalooprs::{feyngen::diagram_generator::evaluate_overall_factor, graph::Graph};
use symbolica::atom::AtomCore;

pub const SCALAR_TREE: &[u64] = &[0x277aa3ff567f454a, 0xe59200ed98b77b6f, 0xe98eb2023fbd4eb4];
pub const STANDARD_MODEL_TREE: &[u64] = &[0x945be7540b7f5dfb];
pub const STANDARD_MODEL_CROSS_SECTION: &[u64] = &[0x6a4905355b86b7c4];
pub const STANDARD_MODEL_VACUUM: &[u64] = &[
    0x6d2e66ecdb7415e4,
    0x72836e116d8ae276,
    0x97a8e5d13fbc4cc4,
    0xa777b2d714039a1a,
    0xd56bf1a35fa06154,
];
pub const SCALAR_SELECTED_PREFACTOR: &[u64] = &[0x16d6cf7f341dad62];
pub const SCALAR_ONE_LOOP: &[u64] = &[
    0x11cde9ac401439af,
    0x127c7739e22ce313,
    0x19ba1292eee536a7,
    0x20ba89c68f1e6086,
    0x29143ccc33a122ff,
    0x3ca1d4455534cce1,
    0x4dc4953dae7f5af7,
    0x4e799f3073589510,
    0x57dce85e9a918b6b,
    0x5dd08ddfc88d3eb5,
    0x6ed1bd64b3732c09,
    0x74b7bdc95a3a01ec,
    0xac556f197fd77958,
    0xb2bb68505ab31df6,
    0xb599d2bae760bc52,
    0xbb67fbe6c0a1a725,
    0xbb79618d09ad1ff0,
    0xc148a82e837393bb,
    0xd5a7ac1ed4612d0d,
    0xe31a1ce90066ec31,
    0xe3f1cf81781f4249,
    0xf207fb2207ec648d,
    0xf96acd583253824d,
    0xff6e6f1cb740596d,
];
pub const SCALAR_ONE_LOOP_FILTERED: &[u64] = &[
    0x127c7739e22ce313,
    0x19ba1292eee536a7,
    0x29143ccc33a122ff,
    0x3ca1d4455534cce1,
    0x4dc4953dae7f5af7,
    0x5dd08ddfc88d3eb5,
    0xb2bb68505ab31df6,
    0xb599d2bae760bc52,
    0xbb79618d09ad1ff0,
    0xc148a82e837393bb,
    0xe31a1ce90066ec31,
    0xe3f1cf81781f4249,
];
pub const FERMION_AMPLITUDES: &[&[u64]] = &[
    &[0x945be7540b7f5dfb],
    &[0x69f7874c4eaabbb1],
    &[0x6259010cfb39ea5d],
    &[0xe05d113e4e87b55f],
];
pub const FERMION_CROSS_SECTION: &[u64] = &[0x67347256bc48119a];
pub const EMPTY: &[u64] = &[];

pub fn graph_fingerprints(graphs: &[Graph]) -> Vec<u64> {
    let mut fingerprints = graphs
        .iter()
        .map(|graph| {
            let mut normalized = graph.clone();
            normalized.name.clear();
            normalized.overall_factor =
                evaluate_overall_factor(graph.overall_factor.as_view()).expand();
            normalized
                .debug_dot()
                .bytes()
                .fold(0xcbf29ce484222325, |hash, byte| {
                    (hash ^ u64::from(byte)).wrapping_mul(0x100000001b3)
                })
        })
        .collect::<Vec<_>>();
    fingerprints.sort_unstable();
    fingerprints
}

pub fn assert_graph_fixture(graphs: &[Graph], expected: &[u64]) {
    assert_eq!(graph_fingerprints(graphs), expected);
}
