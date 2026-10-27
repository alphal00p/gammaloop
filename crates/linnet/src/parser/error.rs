use std::fmt::Display;

use thiserror::Error;

use super::DotSyntaxError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExplicitIdKind {
    Node,
    Edge,
    Hedge,
}

impl Display for ExplicitIdKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Node => write!(f, "node"),
            Self::Edge => write!(f, "edge"),
            Self::Hedge => write!(f, "hedge"),
        }
    }
}

#[derive(Debug, Error)]
pub enum DotParseError {
    #[error("invalid edge direction `{value}`; expected `forward`, `back`, or `none`")]
    InvalidEdgeDirection { value: String },

    #[error("external source endpoint cannot carry half-edge data")]
    ExternalSourceEndpointData,

    #[error("external sink endpoint cannot carry half-edge data")]
    ExternalSinkEndpointData,

    #[error("an edge cannot connect two external nodes")]
    EdgeBetweenExternalNodes,

    #[error("invalid explicit {kind} ID `{value}`: {source}")]
    InvalidExplicitId {
        kind: ExplicitIdKind,
        value: String,
        source: std::num::ParseIntError,
    },

    #[error("duplicate explicit {kind} ID {id}")]
    DuplicateExplicitId { kind: ExplicitIdKind, id: usize },

    #[error("explicit {kind} ID {id} is out of bounds for {len} {kind}s")]
    ExplicitIdOutOfBounds {
        kind: ExplicitIdKind,
        id: usize,
        len: usize,
    },

    #[error("failed to construct the explicit {kind} ID permutation")]
    InvalidExplicitIdPermutation { kind: ExplicitIdKind },
}

#[derive(Debug, Error)]
pub enum HedgeParseError<E, V, H, G> {
    DotVertexDataError(V),

    /// The input is not valid DOT.
    Syntax(DotSyntaxError),

    /// The DOT file could not be read.
    Io(std::io::Error),

    DotEdgeDataError(E),

    HedgeDataFromStringError(H),

    GraphDataFromStringError(G),

    DotParseError(#[source] DotParseError),
}

impl Display for HedgeParseError<(), (), (), ()> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HedgeParseError::Syntax(a) => write!(f, "{a}")?,
            HedgeParseError::Io(a) => write!(f, "{a}")?,
            HedgeParseError::DotParseError(a) => write!(f, "{a}")?,
            _ => unreachable!(),
        }
        Ok(())
    }
}

impl<E, V, H, G> From<DotParseError> for HedgeParseError<E, V, H, G> {
    fn from(error: DotParseError) -> Self {
        Self::DotParseError(error)
    }
}

impl<E, V, H, G> From<DotSyntaxError> for HedgeParseError<E, V, H, G> {
    fn from(error: DotSyntaxError) -> Self {
        Self::Syntax(error)
    }
}

impl<E, V, H, G> From<std::io::Error> for HedgeParseError<E, V, H, G> {
    fn from(error: std::io::Error) -> Self {
        Self::Io(error)
    }
}

pub trait MapGlobal<E, V, H, D> {
    fn map_global<G>(self) -> Result<D, HedgeParseError<E, V, H, G>>;
}

impl<E, V, H, D> MapGlobal<E, V, H, D> for Result<D, HedgeParseError<E, V, H, ()>> {
    fn map_global<G>(self) -> Result<D, HedgeParseError<E, V, H, G>> {
        self.map_err(|e| e.map_global())
    }
}

impl<E, V, H> HedgeParseError<E, V, H, ()> {
    pub fn map_global<G>(self) -> HedgeParseError<E, V, H, G> {
        match self {
            HedgeParseError::GraphDataFromStringError(_) => unreachable!(),
            HedgeParseError::DotEdgeDataError(a) => HedgeParseError::DotEdgeDataError(a),
            HedgeParseError::DotVertexDataError(a) => HedgeParseError::DotVertexDataError(a),
            HedgeParseError::HedgeDataFromStringError(a) => {
                HedgeParseError::HedgeDataFromStringError(a)
            }
            HedgeParseError::Syntax(a) => HedgeParseError::Syntax(a),
            HedgeParseError::Io(a) => HedgeParseError::Io(a),
            HedgeParseError::DotParseError(a) => HedgeParseError::DotParseError(a),
        }
    }
}

pub trait HedgeParseExt<D> {
    type Error;

    fn global<E, V, H>(self) -> Result<D, HedgeParseError<E, V, H, Self::Error>>;
    fn hedge<V, E, G>(self) -> Result<D, HedgeParseError<V, E, Self::Error, G>>;
    fn edge<V, H, G>(self) -> Result<D, HedgeParseError<Self::Error, V, H, G>>;
    fn vertex<E, H, G>(self) -> Result<D, HedgeParseError<E, Self::Error, H, G>>;
}

impl<Err, D> HedgeParseExt<D> for Result<D, Err> {
    type Error = Err;

    fn global<E, V, H>(self) -> Result<D, HedgeParseError<E, V, H, Self::Error>> {
        self.map_err(|e| HedgeParseError::GraphDataFromStringError(e))
    }

    fn hedge<V, E, G>(self) -> Result<D, HedgeParseError<V, E, Self::Error, G>> {
        self.map_err(|e| HedgeParseError::HedgeDataFromStringError(e))
    }

    fn edge<V, H, G>(self) -> Result<D, HedgeParseError<Self::Error, V, H, G>> {
        self.map_err(|e| HedgeParseError::DotEdgeDataError(e))
    }
    fn vertex<E, H, G>(self) -> Result<D, HedgeParseError<E, Self::Error, H, G>> {
        self.map_err(|e| HedgeParseError::DotVertexDataError(e))
    }
}
