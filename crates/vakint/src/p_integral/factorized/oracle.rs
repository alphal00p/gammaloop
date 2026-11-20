//! A naive contraction of a [`FactorizedNumerator`] with idenso, the R*
//! production chain documented by the beta thread
//! (`RStarIntegrand::simplify_view_certified`): cook (reversible encoding with
//! representation payloads) → `infer_indexed` → `simplify_algebra` (gamma
//! `canonical()`, epsilon on, colour off unless asked) → `to_dots` → uncook.
//!
//! It multiplies all factors first and hands the product to idenso: numerator
//! algebra is done by idenso / spenso only.

use idenso::{
    CookMode, CookSettings,
    color::ColorSimplifySettings,
    dirac::GammaSimplifySettings,
    tensor::{AlgebraSettings, SymbolicTensor},
};
use symbolica::atom::{Atom, AtomCore};

use super::{FactorizedError, FactorizedNumerator, Result};

/// Contract a whole expression in source notation to a scalar in idenso's dot
/// notation (`to_dots`), with the production settings. `colour` additionally
/// reduces colour (SU(3) numbers, as the beta's root callback).
pub fn contract_with_idenso(expression: &Atom, colour: bool) -> Result<Atom> {
    let _step = crate::timing::contract();
    Ok(contract_with_idenso_unexpanded(expression, colour)?.expand())
}

/// [`contract_with_idenso`] without its final `expand`: the uncooked
/// closed scalar exactly as idenso's `to_dots` returns it (the R* producer's
/// `simplify_view_certified` stops there too).
pub(crate) fn contract_with_idenso_unexpanded(expression: &Atom, colour: bool) -> Result<Atom> {
    let _step = crate::timing::contract();
    let err = |detail: String| FactorizedError::Oracle { detail };
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let cooked = cooking
        .try_cook(expression.as_view())
        .map_err(|e| err(format!("cooking: {e}")))?;
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::canonical()),
        epsilon: true,
        color: colour.then(|| ColorSimplifySettings::default().with_cof_dimension_invariants()),
        ..Default::default()
    };
    let reduced = SymbolicTensor::infer_indexed(cooked)
        .and_then(|t| t.simplify_algebra(&settings))
        .and_then(|t| t.to_dots())
        .map_err(|e| err(format!("idenso: {e}")))?;
    Ok(cooking.uncook(reduced.expression().as_view()))
}

impl FactorizedNumerator {
    /// `projector.tensor · Π factors`, contracted by [`contract_with_idenso`]
    /// (colour off: the factors are colour-free). Does NOT include the
    /// prefactor nor `(p²)^p2_power`. Slow: a reference oracle only.
    pub fn contract_naively(&self) -> Result<Atom> {
        let mut product = self.factor_product();
        if let Some(projector) = &self.projector {
            product *= &projector.tensor;
        }
        contract_with_idenso(&product, false)
    }
}
