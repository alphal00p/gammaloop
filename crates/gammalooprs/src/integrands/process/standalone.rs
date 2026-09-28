use bincode_trait_derive::{Decode, Encode};
use eyre::{Result, eyre};
use serde::{Deserialize, Serialize};
use symbolica::prelude::Rational;

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneParametricResidueRows<A> {
    pub(crate) parameters: Vec<A>,
    pub(crate) rows: Vec<Vec<Rational>>,
}

impl<A> StandaloneParametricResidueRows<A> {
    pub(crate) fn validate(&self) -> Result<()> {
        if self.rows.is_empty() {
            return Err(eyre!(
                "Standalone parametric residue catalog must contain at least one row"
            ));
        }
        if self
            .rows
            .iter()
            .any(|row| row.len() != self.parameters.len())
        {
            return Err(eyre!(
                "Standalone parametric residue row width differs from its parameter count"
            ));
        }
        Ok(())
    }
}
