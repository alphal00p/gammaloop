use serde::{Deserialize, Serialize};

use crate::{ExpressionTree, GraphOrientation, OrientationId, SurfaceId};

/// One CFF denominator tree and its acyclic graph orientation.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct OrientationExpression {
    pub id: OrientationId,
    pub orientation: GraphOrientation,
    pub expression: ExpressionTree<SurfaceId>,
}

impl OrientationExpression {
    /// Unfold the tree into denominator products, one vector per additive term.
    pub fn denominator_products(&self) -> Vec<Vec<SurfaceId>> {
        self.expression
            .root_to_leaf_paths()
            .into_iter()
            .map(|path| path.into_iter().copied().collect())
            .collect()
    }
}

/// Sum of the expressions associated with all accepted acyclic orientations.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct CffExpression {
    orientations: Vec<OrientationExpression>,
}

impl CffExpression {
    pub(crate) fn new(orientations: Vec<OrientationExpression>) -> Self {
        Self { orientations }
    }

    pub fn orientations(&self) -> &[OrientationExpression] {
        &self.orientations
    }

    pub fn orientation(&self, id: OrientationId) -> Option<&OrientationExpression> {
        self.orientations.get(id.index())
    }

    pub fn unfolded_term_count(&self) -> usize {
        self.orientations
            .iter()
            .map(|orientation| orientation.expression.leaves().count())
            .sum()
    }

    pub fn is_empty(&self) -> bool {
        self.orientations.is_empty()
    }
}
