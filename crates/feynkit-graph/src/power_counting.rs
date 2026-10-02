use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    function, symbol,
};

use crate::DiagramError;

/// Momentum power counting shared by FeynKit diagrams and GammaLoop UV analysis.
///
/// The momentum head is explicit so consumers share the algorithm without
/// depending on each other's symbol registries.
pub trait DOD: AtomCore {
    /// Rescales momentum of edge `eid`, and computes the leading scaling.
    fn edge_dod(&self, momentum: Symbol, eid: usize) -> Result<i32, DiagramError> {
        let arguments = symbol!("feynkit_graph::uv_args___");
        let scale = symbol!("feynkit_graph::uv_rescale");
        self.replace(function!(momentum, eid, arguments))
            .with(function!(momentum, eid, arguments) / scale)
            .trailing_exponent(scale)
    }

    /// Rescales all momenta, and computes the leading scaling.
    fn all_dod(&self, momentum: Symbol) -> Result<i32, DiagramError> {
        let arguments = symbol!("feynkit_graph::uv_args___");
        let scale = symbol!("feynkit_graph::uv_rescale");
        self.replace(function!(momentum, arguments))
            .with(function!(momentum, arguments) / scale)
            .replace(function!(symbol!("UFO::P"), arguments))
            .with(function!(symbol!("UFO::P"), arguments) / scale)
            .replace(function!(symbol!("UFO::PSlash"), arguments))
            .with(function!(symbol!("UFO::PSlash"), arguments) / scale)
            .trailing_exponent(scale)
    }

    /// Return the negative leading exponent as an inverse UV scale tends to zero.
    fn trailing_exponent(&self, scale: Symbol) -> Result<i32, DiagramError> {
        let series = self
            .series(scale, Atom::Zero, 1)
            .map_err(|error| DiagramError::UvPowerCounting(error.to_string()))?;
        let degree = series.get_trailing_exponent();
        if degree.is_integer()
            && let Some(value) = degree.numerator().to_i64().and_then(i64::checked_neg)
            && let Ok(value) = i32::try_from(value)
        {
            return Ok(value);
        }
        Err(DiagramError::UvPowerCounting(format!(
            "non-integral or out-of-range scaling exponent {degree}"
        )))
    }
}

impl<T: AtomCore> DOD for T {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dod() {
        let momentum = crate::momentum_symbol();
        let (e1, e2) = (1, 2);
        let q1 = function!(momentum, e1, Atom::Zero);
        let q2 = function!(momentum, e2, Atom::Zero);
        let atom = (&q1 * &q2 + &q2) / (&q1 * &q1);
        let atom2 = Atom::num(1) / (&q1 * &q1 + symbol!("feynkit_graph::m"));
        let atom3 = &q1 * &q2 + &q2;

        assert_eq!(-1, atom.edge_dod(momentum, e1).unwrap());
        assert_eq!(1, atom.edge_dod(momentum, e2).unwrap());
        assert_eq!(-2, atom2.edge_dod(momentum, e1).unwrap());
        assert_eq!(1, atom3.edge_dod(momentum, e1).unwrap());
        assert_eq!(1, atom3.edge_dod(momentum, e2).unwrap());
        assert_eq!(2, atom3.all_dod(momentum).unwrap());
        assert_eq!(0, (&q1 / &q1).all_dod(momentum).unwrap());
        let ufo = function!(symbol!("UFO::P"), 1, 1) * function!(symbol!("UFO::PSlash"), 2, 1, 1);
        assert_eq!(2, ufo.all_dod(momentum).unwrap());
    }
}
