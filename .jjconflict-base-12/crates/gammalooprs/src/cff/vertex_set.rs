use serde::{Deserialize, Serialize};

const MAX_VERTEX_COUNT: usize = 64;

#[derive(
    Clone,
    Copy,
    Debug,
    PartialEq,
    Eq,
    Hash,
    Serialize,
    Deserialize,
    bincode::Encode,
    bincode::Decode,
)]
pub struct VertexSet {
    vertex_set: u64,
}

impl VertexSet {
    pub(crate) fn from_usize(id: usize) -> Self {
        assert!(id < MAX_VERTEX_COUNT, "Vertex ID out of bounds");

        VertexSet {
            vertex_set: 1 << id,
        }
    }

    pub(crate) fn join(&self, other: &VertexSet) -> VertexSet {
        VertexSet {
            vertex_set: self.vertex_set | other.vertex_set,
        }
    }

    pub(crate) fn dummy() -> Self {
        VertexSet { vertex_set: 0 }
    }
}
