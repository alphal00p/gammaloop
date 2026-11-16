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

use three_dimensional_reps::ThermalDistributionCall;

use crate::{
    cff::expression::OrientationID, integrands::process::param_builder::FnMapEntry, utils::GS,
    uv::Integrands,
};

pub(super) struct ScalarCoefficients {
    variable: Symbol,
    scope: Atom,
    calls: BTreeMap<Atom, Atom>,
    definitions: Vec<Arc<FnMapEntry>>,
    /// Bodies of extracted homogeneous factors, keyed by body, with their tags.
    factors: BTreeMap<Atom, Atom>,
    factor_definitions: Vec<Arc<FnMapEntry>>,
}

/// Head of the tag that names one extracted homogeneous factor.
pub(super) fn factor_symbol() -> Symbol {
    Integrands::factor_symbol()
}

/// Whether `view` is a closed homogeneous-factor family call.
pub(super) fn is_factor_call(view: AtomView<'_>) -> bool {
    Integrands::factor_tag(view).is_some()
}

impl ScalarCoefficients {
    pub(super) fn new(variable: Symbol, scope: Atom) -> Self {
        Self {
            variable,
            scope,
            calls: BTreeMap::new(),
            definitions: Vec::new(),
            factors: BTreeMap::new(),
            factor_definitions: Vec::new(),
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

    /// Degree of `view` as an exactly homogeneous function of the Taylor
    /// variable, `G(t) = t^d G(1)` for positive `t`, or `None` when it is not
    /// certified. Lazy guards, markers, retained calls, tensors, and occupations
    /// that are not soft-compensated thermal distributions are never certified.
    fn degree(&self, view: AtomView<'_>) -> Option<Rational> {
        match view {
            AtomView::Num(_) => Some(Rational::from(0)),
            AtomView::Var(var) => {
                Some(Rational::from(i64::from(var.get_symbol() == self.variable)))
            }
            AtomView::Add(sum) => {
                let mut terms = sum.iter();
                let degree = self.degree(terms.next()?)?;
                terms
                    .all(|term| self.degree(term).as_ref() == Some(&degree))
                    .then_some(degree)
            }
            AtomView::Mul(product) => product
                .iter()
                .try_fold(Rational::from(0), |degree, factor| {
                    Some(degree + self.degree(factor)?)
                }),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let exponent = Rational::try_from(exponent).ok()?;
                Some(self.degree(base)? * exponent)
            }
            AtomView::Fun(call) => {
                let head = call.get_symbol();
                if [
                    OrientationID::symbol(),
                    GS.theta,
                    GS.orientation_delta,
                    Symbol::IF,
                    GS.ct_marker,
                    Self::head(),
                    symbol!("gammalooprs::uv::numerator_family"),
                    factor_symbol(),
                ]
                .contains(&head)
                    || view.is_tensorial(StrictTensorFilter::ContainsReps)
                {
                    return None;
                }
                if head == GS.thermal_distribution
                    && (!view.contains_symbol(self.variable)
                        || !ThermalDistributionCall::try_from(view)
                            .is_ok_and(|call| call.temperature_flag.is_one()))
                {
                    return None;
                }
                if head == GS.on_shell_energy && call.get_nargs() == 2 {
                    // A square root of its invariant; the owner is a tag.
                    return (self.degree(call.get(0))?.is_zero())
                        .then(|| self.degree(call.get(1)))
                        .flatten()
                        .map(|degree| degree / Rational::from(2));
                }
                if head == GS.esurface {
                    // An n-ary sum of its arguments.
                    let mut terms = call.iter();
                    let degree = self.degree(terms.next()?)?;
                    return terms
                        .all(|term| self.degree(term).as_ref() == Some(&degree))
                        .then_some(degree);
                }
                if head == GS.thermal_weight_wrapper && call.get_nargs() == 1 {
                    // Numerically the identity.
                    return self.degree(call.get(0));
                }
                call.iter()
                    .all(|argument| self.degree(argument).is_some_and(|d| d.is_zero()))
                    .then(|| Rational::from(0))
            }
        }
    }

    /// Replace each maximal exactly homogeneous, Taylor-dependent subtree in a
    /// polynomial position by one opaque t-free factor times its power of the
    /// Taylor variable. The rewrite is an identity for positive `t`, so native
    /// series of the result equal those of `atom`, while each factor body is
    /// expanded never and stored once. Factors never enter inverses or
    /// function arguments, so they remain polynomial in every coefficient.
    pub(super) fn retain_homogeneous(&mut self, atom: &Atom) -> Atom {
        match self.homogeneous(atom.as_view()) {
            (Some(degree), _) => self.factor(atom, &degree),
            (None, rewritten) => rewritten.unwrap_or_else(|| atom.clone()),
        }
    }

    /// The certified degree of `view`, or its rewrite when it is not
    /// certified and one of its polynomial parts changed.
    fn homogeneous(&mut self, view: AtomView<'_>) -> (Option<Rational>, Option<Atom>) {
        let parts = match view {
            AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
            AtomView::Mul(product) => product.iter().collect(),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let Ok(exponent) = Rational::try_from(exponent) else {
                    return (None, None);
                };
                if !exponent.is_integer() || exponent.is_negative() {
                    return (self.degree(view), None);
                }
                return match self.homogeneous(base) {
                    (Some(degree), _) => (Some(degree * exponent), None),
                    (None, rewritten) => (
                        None,
                        rewritten.map(|base| base.pow(Atom::num(exponent.clone()))),
                    ),
                };
            }
            _ => return (self.degree(view), None),
        };
        let walked = parts
            .iter()
            .map(|part| self.homogeneous(*part))
            .collect::<Vec<_>>();
        let is_sum = matches!(view, AtomView::Add(_));
        if let Some(degrees) = walked
            .iter()
            .map(|(degree, _)| degree.clone())
            .collect::<Option<Vec<_>>>()
        {
            if !is_sum {
                return (
                    Some(
                        degrees
                            .into_iter()
                            .fold(Rational::from(0), |total, degree| total + degree),
                    ),
                    None,
                );
            }
            if degrees.iter().all(|degree| *degree == degrees[0]) {
                return (degrees.into_iter().next(), None);
            }
        }
        // Mixed: certified Taylor-dependent parts become factors, grouped by
        // degree in a sum and multiplied together in a product. Taylor-free
        // certified parts stay visible to enclosing operations.
        let mut kept = Vec::new();
        let mut groups = BTreeMap::<Rational, Vec<Atom>>::new();
        let mut product_degree = Rational::from(0);
        let mut changed = false;
        for (part, (degree, rewritten)) in parts.into_iter().zip(walked) {
            match degree {
                Some(degree) if part.contains_symbol(self.variable) => {
                    let key = if is_sum {
                        degree
                    } else {
                        product_degree += degree;
                        Rational::from(0)
                    };
                    groups.entry(key).or_default().push(part.to_owned());
                }
                _ => {
                    changed |= rewritten.is_some();
                    kept.push(rewritten.unwrap_or_else(|| part.to_owned()));
                }
            }
        }
        for (degree, group) in groups {
            let (group, degree) = if is_sum {
                (Atom::add_many(group), degree)
            } else {
                (Atom::mul_many(group), product_degree.clone())
            };
            let factor = self.factor(&group, &degree);
            changed |= factor != group;
            kept.push(factor);
        }
        if !changed {
            return (None, None);
        }
        (
            None,
            Some(if is_sum {
                Atom::add_many(kept)
            } else {
                Atom::mul_many(kept)
            }),
        )
    }

    /// `G = t^degree G(1)`, with a large Taylor-dependent `G(1)` shared as one
    /// bare factor tag.
    fn factor(&mut self, atom: &Atom, degree: &Rational) -> Atom {
        if !atom.contains_symbol(self.variable) {
            return atom.clone();
        }
        let body = atom.replace(self.variable).with(Atom::num(1));
        if body.is_zero() {
            return Atom::Zero;
        }
        // Small factors cost less inline than a separate shared body.
        if body.as_view().get_byte_size() < 256 {
            return atom.clone();
        }
        let tag = if let Some(tag) = self.factors.get(&body) {
            tag.clone()
        } else {
            let tag = factor_symbol()
                .call_args([self.scope.clone(), Atom::num(self.factor_definitions.len())]);
            self.factor_definitions.push(Arc::new(FnMapEntry {
                lhs: symbol!("gammalooprs::uv::numerator_family").call_args([tag.clone()]),
                rhs: body.clone(),
                args: vec![],
                tags: vec![tag.clone()],
                inlining: InliningPolicy::Always,
                is_alias: false,
            }));
            self.factors.insert(body, tag.clone());
            tag
        };
        tag * Atom::var(self.variable).pow(Atom::num(degree.clone()))
    }

    /// Hand bare factor tags out as closed parameter-free families.
    pub(super) fn wrap_factors(&self, atom: &Atom) -> Atom {
        if self.factor_definitions.is_empty() {
            return atom.clone();
        }
        atom.replace_map(|view, _, out| {
            if let AtomView::Fun(call) = view
                && call.get_symbol() == factor_symbol()
            {
                **out = symbol!("gammalooprs::uv::numerator_family").call_args([view]);
            }
        })
    }

    pub(super) fn factor_definitions(&self) -> &[Arc<FnMapEntry>] {
        &self.factor_definitions
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
    use symbolica::{function, parse_lit, poly::series::SeriesDepth};

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

    /// A loop-momentum component after an enclosing rescaling by `t`.
    fn component(edge: i64, index: usize, t: &Atom) -> Atom {
        GS.emr_vec.call_args([Atom::num(edge), GS.cind(index)]) / t
    }

    /// A frozen soft energy after an enclosing rescaling: its invariant is
    /// homogeneous of degree zero, as a sum over explicit components.
    fn frozen_energy(t: &Atom, loop_edge: i64) -> Atom {
        let invariant = Atom::add_many((1..=3).map(|index| component(loop_edge, index, t).pow(2)))
            + parse_lit!(frozen_mass).pow(2) / t.pow(2);
        function!(
            GS.on_shell_energy,
            symbol!("S_frozen"),
            t.pow(2) * invariant
        )
    }

    /// Expand energy invariants, so equal energies compare equal.
    fn canonical_energies(atom: &Atom) -> Atom {
        atom.replace_map(|view, _, out| {
            if let AtomView::Fun(call) = view
                && call.get_symbol() == GS.on_shell_energy
                && call.get_nargs() == 2
            {
                **out = function!(GS.on_shell_energy, call.get(0), call.get(1).expand());
            }
        })
    }

    /// The MUV rearrangement of an ordinary energy, which is not homogeneous.
    fn muv_energy(edge: i64, t: &Atom) -> Atom {
        let invariant = -parse_lit!(mUV).pow(2)
            + (parse_lit!(mUV).pow(2) + GS.emr_vec.call_args([Atom::num(edge), GS.cind(1)]).pow(2))
                / t.pow(2);
        function!(GS.on_shell_energy, edge, t.pow(2) * invariant)
    }

    #[test]
    fn homogeneous_degrees_are_certified_only_for_exact_scalings() -> Result<()> {
        test_initialise()?;
        let t_symbol = symbol!("coefficient_t");
        let t = Atom::var(t_symbol);
        let coefficients = ScalarCoefficients::new(t_symbol, Atom::num(3));
        let degree = |atom: &Atom| coefficients.degree(atom.as_view());
        let frozen = frozen_energy(&t, 0);
        assert_eq!(degree(&frozen), Some(Rational::from(0)));
        assert_eq!(degree(&(&frozen / &t)), Some(Rational::from(-1)));
        let scaled = function!(GS.on_shell_energy, 4, t.pow(2) * parse_lit!(x));
        assert_eq!(degree(&scaled), Some(Rational::from(1)));
        // A soft-compensated occupation: vacuum + t^-d N(E/t).
        for order in 0..=2 {
            let energy = function!(GS.on_shell_energy, 4, t.pow(2) * parse_lit!(x)) / &t;
            let occupation = ThermalDistributionCall::atom_at_energy(4, order, 1, -1, 1, &energy);
            assert_eq!(degree(&occupation), Some(Rational::from(0)));
            assert_eq!(
                degree(&(t.pow(Atom::num(-order)) * &occupation)),
                Some(Rational::from(-order))
            );
        }
        // Zero-temperature steps and Taylor-free occupations are barriers.
        let energy = function!(GS.on_shell_energy, 4, t.pow(2) * parse_lit!(x)) / &t;
        assert_eq!(
            degree(&ThermalDistributionCall::atom_at_energy(
                4, 1, 0, -1, 1, &energy
            )),
            None
        );
        assert_eq!(
            degree(&ThermalDistributionCall::atom_at_energy(
                4,
                0,
                1,
                -1,
                1,
                parse_lit!(y)
            )),
            None
        );
        let surface = GS.wrap_esurface(&(&scaled + &t * parse_lit!(z)));
        assert_eq!(degree(&surface), Some(Rational::from(1)));
        assert_eq!(
            degree(&function!(GS.thermal_weight_wrapper, &scaled)),
            Some(Rational::from(1))
        );
        // The MUV rearrangement of an ordinary energy is not homogeneous.
        assert_eq!(degree(&muv_energy(0, &t)), None);
        assert_eq!(degree(&(&t + 1)), None);
        assert_eq!(degree(&function!(symbol!("f"), &t)), None);
        assert_eq!(degree(&t.pow(parse_lit!(y))), None);
        for barrier in [
            function!(Symbol::IF, symbol!("condition"), &t, 0),
            function!(GS.theta, &t),
            function!(symbol!("gammalooprs::uv::numerator_family"), 1),
            function!(factor_symbol(), 1, 2),
            function!(Integrands::scalar_symbol(), 1, 2),
        ] {
            assert_eq!(degree(&barrier), None, "{barrier}");
        }
        Ok(())
    }

    #[test]
    fn homogeneous_factors_are_exact_polynomial_and_series_transparent() -> Result<()> {
        test_initialise()?;
        let t_symbol = symbol!("coefficient_t");
        let t = Atom::var(t_symbol);
        let frozen = frozen_energy(&t, 0);
        // A frozen-like product of degree -1: the soft jet of a hard line.
        let dot = Atom::add_many((1..=3).map(|i| component(0, i, &t) * component(2, i, &t)));
        let square = Atom::add_many((1..=3).map(|i| component(2, i, &t).pow(2)));
        let frozen_product = t.pow(3) * dot.pow(2) / frozen.pow(3) + &t * &square / &frozen;
        // The cograph energy carries the MUV rearrangement; the eikonal
        // denominator mixes it with a frozen, homogeneous shift.
        let cograph = muv_energy(5, &t);
        let shift = &dot / &frozen;
        let eikonal = (&cograph / &t + &shift / &t).pow(-1);
        let small = &t * GS.emr_vec.call_args([Atom::num(5), GS.cind(1)]);
        let source = &eikonal * &frozen_product + &cograph * &small;

        let mut coefficients = ScalarCoefficients::new(t_symbol, Atom::num(4));
        let skeleton = coefficients.retain_homogeneous(&source);
        assert_eq!(coefficients.factor_definitions().len(), 1);
        assert!(skeleton.as_view().get_byte_size() < source.as_view().get_byte_size());
        // Factors only in polynomial positions.
        let mut hidden = false;
        skeleton.visitor(&mut |view| {
            match view {
                AtomView::Fun(call) if call.get_symbol() != factor_symbol() => {
                    hidden |= view.contains_symbol(factor_symbol());
                    return false;
                }
                AtomView::Pow(power) => {
                    let negative = Rational::try_from(power.get_base_exp().1)
                        .is_ok_and(|exponent| exponent.is_negative());
                    hidden |= negative && view.contains_symbol(factor_symbol());
                }
                _ => {}
            }
            true
        });
        assert!(!hidden, "{skeleton}");
        // The small Taylor-dependent cograph term stays as it was.
        assert!(skeleton.contains(&(&cograph * &small)));
        let wrapped = coefficients.wrap_factors(&skeleton);
        let inline = coefficients
            .factor_definitions()
            .iter()
            .map(|entry| entry.replacement())
            .collect::<Vec<_>>();
        let restored = wrapped.replace_multiple(&inline);
        assert!(!restored.contains_symbol(factor_symbol()));
        let difference = (canonical_energies(&restored) - canonical_energies(&source))
            .expand()
            .together()
            .cancel();
        assert!(difference.is_zero(), "{difference}");
        // Native series of the skeleton, with the factor bodies restored
        // afterwards, equal those of the source.
        for depth in [-1, 0, 2] {
            let expected = source
                .series(
                    t_symbol,
                    Atom::Zero,
                    SeriesDepth::absolute(Rational::from(depth)),
                )?
                .to_atom();
            let actual = wrapped
                .series(
                    t_symbol,
                    Atom::Zero,
                    SeriesDepth::absolute(Rational::from(depth)),
                )?
                .to_atom()
                .replace_multiple(&inline);
            let difference = (canonical_energies(&actual) - canonical_energies(&expected))
                .expand()
                .together()
                .cancel();
            assert!(difference.is_zero(), "depth {depth}: {difference}");
        }
        // Identical bodies share one factor; a second pass finds no new one.
        let again = coefficients.retain_homogeneous(&source);
        assert_eq!(again, skeleton);
        assert_eq!(coefficients.factor_definitions().len(), 1);
        assert_eq!(coefficients.retain_homogeneous(&skeleton), skeleton);
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
