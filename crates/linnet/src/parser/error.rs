use std::fmt::Display;

use dot_parser::ast::{GraphFromFileError, PestError};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum HedgeParseError<'a, E, V, H, G> {
    DotVertexDataError(V),

    GraphFromFile(GraphFromFileError<'a>),

    ParseError(PestError),

    DotEdgeDataError(E),

    HedgeDataFromStringError(H),

    GraphDataFromStringError(G),
}

impl<'a> Display for HedgeParseError<'a, (), (), (), ()> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HedgeParseError::ParseError(a) => write!(f, "{a}")?,
            HedgeParseError::GraphFromFile(a) => write!(f, "{a}")?,
            _ => unreachable!(),
        }
        Ok(())
    }
}

impl<'a, E, V, H, G> From<GraphFromFileError<'a>> for HedgeParseError<'a, E, V, H, G> {
    fn from(e: GraphFromFileError<'a>) -> Self {
        HedgeParseError::GraphFromFile(e)
    }
}

impl<E, V, H, G> From<PestError> for HedgeParseError<'_, E, V, H, G> {
    fn from(e: PestError) -> Self {
        HedgeParseError::ParseError(e)
    }
}

pub trait MapGlobal<'a, E, V, H, D> {
    #[allow(clippy::result_large_err)]
    fn map_global<G>(self) -> Result<D, HedgeParseError<'a, E, V, H, G>>;
}

impl<'a, E, V, H, D> MapGlobal<'a, E, V, H, D> for Result<D, HedgeParseError<'a, E, V, H, ()>> {
    fn map_global<G>(self) -> Result<D, HedgeParseError<'a, E, V, H, G>> {
        self.map_err(|e| e.map_global())
    }
}

impl<'a, E, V, H> HedgeParseError<'a, E, V, H, ()> {
    pub fn map_global<G>(self) -> HedgeParseError<'a, E, V, H, G> {
        match self {
            HedgeParseError::GraphDataFromStringError(_) => unreachable!(),
            HedgeParseError::DotEdgeDataError(a) => HedgeParseError::DotEdgeDataError(a),
            HedgeParseError::DotVertexDataError(a) => HedgeParseError::DotVertexDataError(a),
            HedgeParseError::HedgeDataFromStringError(a) => {
                HedgeParseError::HedgeDataFromStringError(a)
            }
            HedgeParseError::GraphFromFile(a) => HedgeParseError::GraphFromFile(a),
            HedgeParseError::ParseError(a) => HedgeParseError::ParseError(a),
        }
    }
}

#[allow(clippy::result_large_err)]
pub trait HedgeParseExt<'a, D> {
    type Error;

    fn global<E, V, H>(self) -> Result<D, HedgeParseError<'a, E, V, H, Self::Error>>;
    fn hedge<V, E, G>(self) -> Result<D, HedgeParseError<'a, V, E, Self::Error, G>>;
    fn edge<V, H, G>(self) -> Result<D, HedgeParseError<'a, Self::Error, V, H, G>>;
    fn vertex<E, H, G>(self) -> Result<D, HedgeParseError<'a, E, Self::Error, H, G>>;
}

impl<'a, Err, D> HedgeParseExt<'a, D> for Result<D, Err> {
    type Error = Err;

    fn global<E, V, H>(self) -> Result<D, HedgeParseError<'a, E, V, H, Self::Error>> {
        self.map_err(|e| HedgeParseError::GraphDataFromStringError(e))
    }

    fn hedge<V, E, G>(self) -> Result<D, HedgeParseError<'a, V, E, Self::Error, G>> {
        self.map_err(|e| HedgeParseError::HedgeDataFromStringError(e))
    }

    fn edge<V, H, G>(self) -> Result<D, HedgeParseError<'a, Self::Error, V, H, G>> {
        self.map_err(|e| HedgeParseError::DotEdgeDataError(e))
    }
    fn vertex<E, H, G>(self) -> Result<D, HedgeParseError<'a, E, Self::Error, H, G>> {
        self.map_err(|e| HedgeParseError::DotVertexDataError(e))
    }
}
