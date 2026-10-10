//! Parametric scalar coefficients retained between residue Taylor operations.
//!
//! Native series owns Laurent valuations and precision. Physical functions remain
//! complete actual arguments; templates contain only arithmetic in scalar formals.

use std::{collections::BTreeMap, sync::Arc};

#[cfg(test)]
use color_eyre::eyre::{Result, ensure};
use spenso::network::parsing::{AtomStructureExt, StrictTensorFilter};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    domains::rational::Rational,
    evaluate::InliningPolicy,
    symbol,
};

use crate::{
    cff::expression::OrientationID, integrands::process::param_builder::FnMapEntry, utils::GS,
    uv::Integrands,
};

pub(super) struct ScalarCoefficients {
    variable: Symbol,
    scope: Atom,
    calls: BTreeMap<Atom, Atom>,
    definitions: Vec<Arc<FnMapEntry>>,
}

impl ScalarCoefficients {
    pub(super) fn new(variable: Symbol, scope: Atom) -> Self {
        Self {
            variable,
            scope,
            calls: BTreeMap::new(),
            definitions: Vec::new(),
        }
    }

    fn head() -> Symbol {
        Integrands::scalar_symbol()
    }

    fn scalar(&self, atom: AtomView<'_>) -> bool {
        !atom.contains_symbol(Self::head())
            && !atom.contains_symbol(self.variable)
            && !atom.contains_symbol(symbol!("gammalooprs::uv::numerator_family"))
            && !atom.is_tensorial(StrictTensorFilter::ContainsReps)
    }

    fn call(&mut self, atom: Atom) -> Atom {
        if let Some(call) = self.calls.get(&atom) {
            return call.clone();
        }
        // Small coefficients cost less than a signature and a separate body.
        if atom.as_view().get_byte_size() < 256
            || !(atom.contains_symbol(GS.on_shell_energy) || atom.contains_symbol(GS.esurface))
        {
            return atom;
        }
        let mut arguments = BTreeMap::new();
        let mut values = Vec::new();
        let mut parameters = Vec::new();
        let body = atom.replace_map(|view, _, out| {
            if matches!(view, AtomView::Var(_) | AtomView::Fun(_)) {
                let parameter = arguments.entry(view.to_owned()).or_insert_with(|| {
                    let parameter = symbol!(&format!(
                        "gammalooprs::uv::scalar_coefficient_p{}",
                        values.len()
                    ));
                    values.push(view.to_owned());
                    parameters.push(parameter.into());
                    Atom::var(parameter)
                });
                **out = parameter.clone();
            }
        });
        // Energies remain complete owner/invariant calls in the arguments.
        // Other physical functions and variables are parameters too; the body
        // contains only arithmetic in formals, never hidden momentum or mass.
        let tags = vec![self.scope.clone(), Atom::num(self.definitions.len())];
        let call = Self::head().call_args(tags.iter().chain(&values));
        if call.as_view().get_byte_size() >= atom.as_view().get_byte_size() {
            return atom;
        }
        let lhs = Self::head().call_args(
            tags.iter()
                .cloned()
                .chain(parameters.iter().cloned().map(Atom::from)),
        );
        self.definitions.push(Arc::new(FnMapEntry {
            lhs,
            rhs: body,
            args: parameters,
            tags,
            inlining: InliningPolicy::Always,
            is_alias: false,
        }));
        self.calls.insert(atom, call.clone());
        call
    }

    pub(super) fn retain(&mut self, atom: &Atom) -> Atom {
        // Keep lazy evaluation boundaries intact. Existing calls stay opaque,
        // while new scalar factors elsewhere can enter the same flat store.
        if [
            OrientationID::symbol(),
            GS.theta,
            GS.orientation_delta,
            Symbol::IF,
        ]
        .into_iter()
        .any(|head| atom.contains_symbol(head))
        {
            return atom.clone();
        }
        atom.replace_map(|view, _, out| {
            if self.scalar(view) {
                **out = self.call(view.to_owned());
            } else if matches!(view, AtomView::Fun(_))
                || matches!(view, AtomView::Pow(power) if !Rational::try_from(power.get_exp())
                    .is_ok_and(|power| power.is_integer() && !power.is_negative()))
            {
                // Calls must be polynomial factors in the root. In particular,
                // never hide leading cancellations under an enclosing inverse.
                **out = view.to_owned();
            } else if let AtomView::Mul(product) = view {
                let (scalar, other): (Vec<_>, Vec<_>) =
                    product.iter().partition(|factor| self.scalar(*factor));
                let scalar = Atom::mul_many(scalar);
                let retained = self.call(scalar.clone());
                if retained != scalar {
                    **out = retained * Atom::mul_many(other);
                }
            }
        })
    }

    #[cfg(test)]
    pub(super) fn restore(&self, atom: &Atom) -> Result<Atom> {
        let replacements = self
            .definitions
            .iter()
            .map(|entry| entry.replacement())
            .collect::<Vec<_>>();
        let restored = atom.replace_multiple(&replacements);
        ensure!(
            !restored.contains_symbol(Self::head()),
            "unresolved scalar Taylor coefficient"
        );
        Ok(restored)
    }

    pub(super) fn definitions(&self) -> &[Arc<FnMapEntry>] {
        &self.definitions
    }

    pub(super) fn len(&self) -> usize {
        self.definitions.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{initialisation::test_initialise, numerator::symbolica_ext::NumeratorAtomExt};
    use symbolica::{function, parse_lit};

    fn coefficient() -> Atom {
        let q = parse_lit!(momentum(q));
        let mass = parse_lit!(mass);
        let energy = function!(GS.on_shell_energy, 7, q.pow(2) + mass.pow(2));
        let surface = GS.wrap_esurface(&(&energy + parse_lit!(shift)));
        Atom::add_many((1..=12).map(|power| energy.pow(-power) * surface.pow(-power - 1)))
    }

    #[test]
    fn shared_coefficients_are_parametric_and_keep_energy_owners_visible() -> Result<()> {
        test_initialise()?;
        let t = symbol!("coefficient_t");
        let source = coefficient();
        let mut coefficients = ScalarCoefficients::new(t, Atom::num(0));
        let retained = coefficients.retain(&source);
        assert!(retained.as_view().get_byte_size() < source.as_view().get_byte_size());
        assert_eq!(coefficients.len(), 1);
        assert_eq!(coefficients.retain(&source), retained);
        assert_eq!(coefficients.retain(&retained), retained);
        assert_eq!(coefficients.len(), 1);
        assert_eq!(coefficients.restore(&retained)?, source);
        assert!(retained.contains_symbol(GS.on_shell_energy));
        for (from, to) in [
            (
                parse_lit!(momentum(q)),
                parse_lit!(coefficient_t * momentum(q) + momentum(p)),
            ),
            (parse_lit!(mass), parse_lit!(vacuum_mass)),
            (parse_lit!(shift), parse_lit!(shift + external_energy)),
        ] {
            assert_eq!(
                coefficients.restore(&retained.replace(from.to_pattern()).with(to.to_pattern()))?,
                source.replace(from.to_pattern()).with(to.to_pattern())
            );
        }
        Ok(())
    }

    #[test]
    fn shared_coefficients_preserve_regular_composition_and_enclosing_taylor() -> Result<()> {
        test_initialise()?;
        let t = symbol!("coefficient_t");
        let u = symbol!("coefficient_u");
        let source = coefficient();
        let family = function!(symbol!("gammalooprs::uv::numerator_family"), 42);
        let root = &source * &family / Atom::var(t).pow(2);
        let mut coefficients = ScalarCoefficients::new(t, Atom::num(1));
        let retained = coefficients.retain(&root);
        let jet =
            parse_lit!(n0 + coefficient_t * n1 + coefficient_t ^ 2 * n2 + coefficient_t ^ 3 * n3);
        let finish = |root: &Atom| {
            root.replace(family.to_pattern())
                .with(jet.to_pattern())
                .series_preserving_factors(t, Atom::Zero.as_view(), 0, &[])
        };
        let result = coefficients.restore(&finish(&retained)?)?;
        let reference = finish(&root)?;
        // This scalar oracle contains no graph numerator to distribute.
        assert!((&result - &reference).expand().is_zero());
        // A later Taylor operation first receives resolved coefficients, so
        // derivatives act on the actual shifted momenta and vacuum masses.
        let deform = |atom: Atom| {
            atom.replace(parse_lit!(momentum(q)))
                .with(parse_lit!(momentum(q) + coefficient_u * momentum(p)))
        };
        let actual = deform(result).series(u, Atom::Zero, 1)?.to_atom();
        let expected = deform(reference).series(u, Atom::Zero, 1)?.to_atom();
        assert!((actual - expected).expand().is_zero());
        Ok(())
    }

    #[test]
    fn shared_coefficients_do_not_hide_taylor_dependence_or_lazy_guards() -> Result<()> {
        test_initialise()?;
        let t = symbol!("coefficient_t");
        let source = coefficient();
        let mut coefficients = ScalarCoefficients::new(t, Atom::num(2));
        let dependent = source
            .replace(parse_lit!(momentum(q)))
            .with(parse_lit!(coefficient_t * momentum(q)));
        assert_eq!(coefficients.retain(&dependent), dependent);
        let guarded = function!(Symbol::IF, symbol!("condition"), &source, 0);
        assert_eq!(coefficients.retain(&guarded), guarded);
        assert_eq!(coefficients.len(), 0);
        Ok(())
    }
}
