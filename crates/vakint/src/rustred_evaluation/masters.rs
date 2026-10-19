//! Runtime master coefficients, with availability checked after exact cancellation.
//!
//! Input masters are unit-mass integrals in the package's FMFT/MATAD convention.
//! The caller restores target mass powers and measure normalization exactly once.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::sync::Arc;

use serde::Deserialize;
use sha2::{Digest, Sha256};
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::coefficient::CoefficientView;
use symbolica::domains::float::{Complex, Float};
use symbolica::domains::rational::Rational;
use symbolica::id::Replacement;
use symbolica::poly::series::SeriesDepth;

use super::data::{LoopData, checked_bytes, failure};
use super::{RustRedEvaluationError, RustRedMasterMode};
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{VakintError, VakintSettings, get_integer_from_atom};

#[derive(Deserialize)]
struct FileReference {
    path: String,
    sha256: String,
}

#[derive(Deserialize)]
struct CoefficientMetadata {
    epsilon_power: i64,
    symbolic: bool,
    numerical: bool,
    significant_digits: Option<u32>,
}

#[derive(Deserialize)]
struct MasterMetadata {
    master_id: String,
    loop_count: usize,
    family_sha256: String,
    powers: Vec<i64>,
    parent_graph: String,
    #[serde(default)]
    exact_basis_expression: Option<String>,
    coefficients: Vec<CoefficientMetadata>,
    files: BTreeMap<String, FileReference>,
}

#[derive(Debug)]
struct Coefficient {
    exact: Option<String>,
    numerical: Option<String>,
    significant_digits: Option<u32>,
}

#[derive(Debug)]
struct Master {
    identity: String,
    loop_count: usize,
    coefficients: BTreeMap<i64, Coefficient>,
    graph: String,
}

#[derive(Debug)]
pub(super) struct MasterCatalog {
    version: String,
    coefficient_tag: Atom,
    masters: BTreeMap<Atom, Master>,
    family_graphs: BTreeMap<Atom, String>,
    basis_aliases: BTreeMap<Atom, Atom>,
}

pub(crate) fn combine_exact(expression: Atom) -> Atom {
    let expression = expression.expand();
    let mut has_float = false;
    expression.visitor(&mut |part| {
        if let AtomView::Num(number) = part
            && matches!(number.get_coeff_view(), CoefficientView::Float(_, _))
        {
            has_float = true;
        }
        !has_float
    });
    // Symbolica's rational canonicalization requires exact coefficients.
    // In particular, `together` currently corrupts powers in sums with
    // floating coefficients, so numerical expressions stay expanded.
    if has_float {
        expression
    } else {
        expression.together()
    }
}

impl MasterCatalog {
    pub(super) fn load(
        loops: &BTreeMap<usize, Arc<LoopData>>,
        version: &str,
    ) -> Result<Self, RustRedEvaluationError> {
        let mut catalog = Self {
            version: version.into(),
            coefficient_tag: Atom::Zero,
            masters: BTreeMap::new(),
            family_graphs: BTreeMap::new(),
            basis_aliases: BTreeMap::new(),
        };
        let mut basis_expressions = Vec::new();
        let mut identity = Sha256::new();
        identity.update(version.as_bytes());
        for (&count, data) in loops {
            let mut fingerprints = BTreeMap::new();
            for family in &data.manifest.families {
                let fingerprint = family["fingerprint"]
                    .as_str()
                    .ok_or_else(|| failure("missing family fingerprint"))?;
                let digest = format!("{:x}", Sha256::digest(fingerprint.as_bytes()));
                if family["sha256"].as_str() != Some(digest.as_str()) {
                    return Err(failure(
                        "master family digest disagrees with its full fingerprint",
                    ));
                }
                fingerprints.insert(digest, fingerprint);
                let graph = family["graph"]
                    .as_str()
                    .ok_or_else(|| failure("missing family graph"))?;
                catalog.family_graphs.insert(
                    Self::family_symbol(fingerprint),
                    Self::graph_location(data, graph)?,
                );
            }
            for value in &data.manifest.masters {
                // A mutable local checkout can keep its version label while
                // changing precision or exactness metadata. Bind retained
                // coefficients to those capabilities as well as file hashes.
                identity.update(serde_json::to_vec(value).map_err(failure)?);
                let metadata: MasterMetadata =
                    serde_json::from_value(value.clone()).map_err(failure)?;
                if metadata.loop_count != count {
                    return Err(failure("master loop count differs from its package"));
                }
                let fingerprint = fingerprints
                    .get(&metadata.family_sha256)
                    .ok_or_else(|| failure("master references an undeclared family"))?;
                let symbol = catalog.symbol(fingerprint, &metadata.powers);
                if let Some(expression) = &metadata.exact_basis_expression {
                    basis_expressions.push((
                        count,
                        symbol.clone(),
                        vk_parse!(expression).map_err(failure)?,
                    ));
                }
                let symbolic = Self::read_coefficients(data, metadata.files.get("symbolic.json"))?;
                let numerical =
                    Self::read_coefficients(data, metadata.files.get("numerical.json"))?;
                let mut coefficients = BTreeMap::new();
                for entry in metadata.coefficients {
                    let exact =
                        if entry.symbolic {
                            Some(symbolic.get(&entry.epsilon_power).cloned().ok_or_else(|| {
                                failure("exact coefficient has no symbolic payload")
                            })?)
                        } else {
                            None
                        };
                    let number = if entry.numerical {
                        if entry.significant_digits.is_none_or(|digits| digits == 0) {
                            return Err(failure(
                                "numerical coefficient lacks validated significant digits",
                            ));
                        }
                        Some(
                            numerical
                                .get(&entry.epsilon_power)
                                .cloned()
                                .ok_or_else(|| failure("numerical coefficient has no payload"))?,
                        )
                    } else {
                        None
                    };
                    if exact.is_none() && number.is_none() {
                        return Err(failure("available coefficient has no representation"));
                    }
                    if coefficients
                        .insert(
                            entry.epsilon_power,
                            Coefficient {
                                exact,
                                numerical: number,
                                significant_digits: entry.significant_digits,
                            },
                        )
                        .is_some()
                    {
                        return Err(failure("duplicate master coefficient power"));
                    }
                }
                let master = Master {
                    identity: metadata.master_id,
                    loop_count: count,
                    coefficients,
                    graph: Self::graph_location(data, &metadata.parent_graph)?,
                };
                if catalog.masters.insert(symbol, master).is_some() {
                    return Err(failure("duplicate family-qualified master identity"));
                }
            }
        }
        let tag = format!("rustred_data_{:x}", identity.finalize());
        catalog.coefficient_tag = Atom::var(vk_symbol!(&tag));
        catalog.load_basis_aliases(basis_expressions);
        Ok(catalog)
    }

    fn load_basis_aliases(&mut self, expressions: Vec<(usize, Atom, Atom)>) {
        let epsilon = Atom::var(vk_symbol!("ep"));
        let mut identical = BTreeMap::new();
        let mut primitives = BTreeMap::new();
        for (loops, master, expression) in &expressions {
            if !self.masters[master].coefficients.is_empty()
                && Self::finite_precision(expression).is_none()
            {
                // Equality of complete exact basis expressions also certifies
                // aliases containing Gamma functions or products of periods.
                // No function evaluation or numerical comparison is needed.
                let representative = identical
                    .entry((*loops, expression.clone()))
                    .or_insert_with(|| master.clone());
                if representative != master {
                    self.basis_aliases
                        .insert(master.clone(), representative.clone());
                }
            }
            if matches!(expression.as_view(), AtomView::Var(_))
                && expression != &epsilon
                && !self.masters[master].coefficients.is_empty()
            {
                primitives
                    .entry((*loops, expression.clone()))
                    .or_insert_with(|| master.clone());
            }
        }
        for (loops, master, expression) in expressions {
            let mut exact_rational_relation = true;
            expression.visitor(&mut |part| {
                match part {
                    AtomView::Var(_) => {
                        exact_rational_relation &= part == epsilon.as_view()
                            || primitives.contains_key(&(loops, part.to_owned()));
                    }
                    AtomView::Fun(_) => exact_rational_relation = false,
                    AtomView::Num(number) => {
                        exact_rational_relation &=
                            !matches!(number.get_coeff_view(), CoefficientView::Float(_, _));
                    }
                    _ => {}
                }
                exact_rational_relation
            });
            if !exact_rational_relation {
                continue;
            }
            let replacements: Vec<_> = primitives
                .iter()
                .filter_map(|((count, primitive), target)| {
                    if *count == loops {
                        Some(Replacement::new(primitive.clone(), target.clone()))
                    } else {
                        None
                    }
                })
                .collect();
            let normalized = expression.replace_multiple(&replacements);
            if normalized != master {
                self.basis_aliases.insert(master, normalized);
            }
        }
    }

    fn read_coefficients(
        data: &LoopData,
        file: Option<&FileReference>,
    ) -> Result<BTreeMap<i64, String>, RustRedEvaluationError> {
        let Some(file) = file else {
            return Ok(BTreeMap::new());
        };
        let bytes = data.read(&file.path)?;
        checked_bytes(&bytes, &file.sha256)?;
        let rows: Vec<(i64, String)> = serde_json::from_slice(&bytes).map_err(failure)?;
        let count = rows.len();
        let result: BTreeMap<_, _> = rows.into_iter().collect();
        if count != result.len() {
            return Err(failure("duplicate coefficient in master payload"));
        }
        Ok(result)
    }

    fn graph_location(data: &LoopData, path: &str) -> Result<String, RustRedEvaluationError> {
        data.read(path)?;
        Ok(data
            .path(path)?
            .map(|path| path.display().to_string())
            .unwrap_or_else(|| {
                format!(
                    "https://github.com/alphal00p/vakint-external-data/blob/{}/{}L/{path}",
                    data.manifest.version, data.manifest.loop_count
                )
            }))
    }

    fn family_symbol(fingerprint: &str) -> Atom {
        Atom::var(vk_symbol!(format!(
            "rustred_family_{:x}",
            Sha256::digest(fingerprint.as_bytes())
        )))
    }

    pub(super) fn symbol(&self, family_fingerprint: &str, powers: &[i64]) -> Atom {
        let mut arguments = vec![Self::family_symbol(family_fingerprint)];
        arguments.extend(powers.iter().copied().map(Atom::num));
        vk_symbol!("RustRedMaster").call(&arguments)
    }

    fn coefficient_symbol(&self, master: &Atom, power: i64) -> Atom {
        vk_symbol!("RustRedCoefficient").call((&self.coefficient_tag, master, power))
    }

    fn referenced_masters(expression: &Atom) -> BTreeSet<Atom> {
        let mut result = BTreeSet::new();
        expression.visitor(&mut |part| {
            if let AtomView::Fun(function) = part {
                if function.get_symbol() == vk_symbol!("RustRedCoefficient") {
                    return false;
                }
                if function.get_symbol() == vk_symbol!("RustRedMaster") {
                    result.insert(part.to_owned());
                    return false;
                }
            }
            true
        });
        result
    }

    fn unresolved(&self, expression: &Atom) -> Vec<String> {
        Self::referenced_masters(expression)
            .iter()
            .filter_map(|symbol| match self.masters.get(symbol) {
                Some(master) if !master.coefficients.is_empty() => None,
                Some(master) => Some(format!(
                    "{}: declared terminal has no computed coefficients; parent graph {}",
                    master.identity, master.graph
                )),
                None => {
                    let graph = match symbol.as_view() {
                        AtomView::Fun(function) => function
                            .iter()
                            .next()
                            .and_then(|family| self.family_graphs.get(&family.to_owned())),
                        _ => None,
                    };
                    Some(format!(
                        "{}: undeclared residual master; parent graph {}",
                        symbol.to_canonical_string(),
                        graph
                            .map(String::as_str)
                            .unwrap_or("not present in this data snapshot")
                    ))
                }
            })
            .collect()
    }

    pub(super) fn finalize(
        &self,
        raw: Atom,
        settings: &VakintSettings,
        loop_count: usize,
        mode: RustRedMasterMode,
    ) -> Result<Atom, VakintError> {
        let direct = self.finalize_independent(raw.clone(), settings, loop_count, mode);
        if matches!(&direct, Err(error) if !matches!(error,
            VakintError::RustRedEvaluation(RustRedEvaluationError::MasterCapability { .. })))
            || mode == RustRedMasterMode::Retain
            || self.basis_aliases.is_empty()
        {
            return direct;
        }
        // Prefer authenticated exact relations even when independently rounded
        // coefficient rows are available: a large prefactor must not amplify
        // their rounding differences in a mathematically zero combination.
        // Keep direct data as a fallback when an alias has additional orders.
        // Corruption, parsing and algebra errors are never hidden by this path.
        let replacements: Vec<_> = self
            .basis_aliases
            .iter()
            .map(|(master, basis)| {
                let basis = basis
                    .replace(vk_parse!("ep").unwrap().to_pattern())
                    .with(Atom::var(vk_symbol!(&settings.epsilon_symbol)));
                Replacement::new(master.clone(), basis)
            })
            .collect();
        let normalized = combine_exact(raw.replace_multiple(&replacements));
        if normalized == raw {
            return direct;
        }
        match self.finalize_independent(normalized, settings, loop_count, mode) {
            Err(VakintError::RustRedEvaluation(RustRedEvaluationError::MasterCapability {
                ..
            })) => direct,
            result => result,
        }
    }

    fn finalize_independent(
        &self,
        raw: Atom,
        settings: &VakintSettings,
        loop_count: usize,
        mode: RustRedMasterMode,
    ) -> Result<Atom, VakintError> {
        let raw = combine_exact(raw);
        let unresolved = self.unresolved(&raw);
        if !unresolved.is_empty() {
            let detail = format!(
                "Unresolved RustRed masters (data {}): {}",
                self.version,
                unresolved.join("; ")
            );
            if mode == RustRedMasterMode::Numerical {
                return Err(failure(detail).into());
            }
            log::warn!("{detail}");
        }
        if mode == RustRedMasterMode::Retain {
            return Ok(raw);
        }
        let epsilon = vk_symbol!(&settings.epsilon_symbol);
        let target = settings
            .number_of_terms_in_epsilon_expansion
            .checked_sub(i64::try_from(loop_count).map_err(failure)?)
            .and_then(|value| value.checked_sub(1))
            .ok_or_else(|| failure("epsilon expansion depth overflow"))?;
        let symbols: Vec<_> = Self::referenced_masters(&raw).into_iter().collect();
        let terms = raw.coefficient_list::<i32>(&symbols);
        let mut required = BTreeMap::<Atom, i64>::new();
        let mut known = Atom::Zero;
        let mut retained = Atom::Zero;
        for (monomial, coefficient) in terms {
            let coefficient = combine_exact(coefficient);
            if coefficient.is_zero() {
                continue;
            }
            let mut multiplicities = BTreeMap::new();
            Self::master_multiplicities(monomial.as_view(), &mut multiplicities)?;
            if multiplicities.keys().any(|symbol| {
                self.masters
                    .get(symbol)
                    .is_none_or(|master| master.coefficients.is_empty())
            }) {
                retained += monomial * coefficient;
                continue;
            }
            let series = coefficient
                .series(epsilon, Atom::Zero.as_view(), SeriesDepth::relative(1))
                .map_err(failure)?;
            let exponent = series.get_trailing_exponent();
            if !exponent.is_integer() {
                return Err(failure("fractional epsilon power in reduction coefficient").into());
            }
            let valuation = exponent
                .numerator()
                .to_i64()
                .ok_or_else(|| failure("epsilon valuation exceeds integer range"))?;
            let minimum = multiplicities
                .iter()
                .try_fold(0_i64, |sum, (symbol, power)| {
                    let loops = i64::try_from(self.masters[symbol].loop_count).map_err(failure)?;
                    sum.checked_sub(
                        loops
                            .checked_mul(*power)
                            .ok_or_else(|| failure("master product order overflow"))?,
                    )
                    .ok_or_else(|| failure("master product order overflow"))
                })?;
            for symbol in multiplicities.keys() {
                let own_minimum =
                    -i64::try_from(self.masters[symbol].loop_count).map_err(failure)?;
                let last = target
                    .checked_sub(valuation)
                    .and_then(|value| value.checked_sub(minimum))
                    .and_then(|value| value.checked_add(own_minimum))
                    .ok_or_else(|| failure("required master epsilon order overflow"))?;
                required
                    .entry(symbol.clone())
                    .and_modify(|current| *current = (*current).max(last))
                    .or_insert(last);
            }
            known += monomial * coefficient;
        }
        let mut exact_values = BTreeMap::new();
        for (symbol, last) in required {
            let master = &self.masters[&symbol];
            let first = -(master.loop_count as i64);
            if last.saturating_sub(first) > 4096 {
                return Err(failure(format!("master {} requires epsilon powers {first}..={last}; expansion work limit is 4096 coefficients (data {})", master.identity, self.version)).into());
            }
            let mut expansion = Atom::Zero;
            for power in first..=last {
                let placeholder = self.coefficient_symbol(&symbol, power);
                if let Some(coefficient) = master.coefficients.get(&power)
                    && let Some(formula) = &coefficient.exact
                {
                    let exact = vk_parse!(formula).map_err(failure)?;
                    // Keep finite embedded constants behind a named coefficient
                    // until their precision has been checked in numerical mode.
                    if Self::finite_precision(&exact).is_none() {
                        exact_values.insert(
                            placeholder.clone(),
                            crate::hpl::normalize_constants(exact.as_view())?,
                        );
                    }
                }
                expansion += placeholder * Atom::var(epsilon).pow(Atom::num(power));
            }
            known = known.replace(symbol.to_pattern()).with(expansion);
        }
        // Exact zero coefficients and exact relations must cancel before an
        // absent/high-order or low-precision coefficient can cause an error.
        for (source, value) in exact_values {
            known = known.replace(source.to_pattern()).with(value);
        }
        let mut result = combine_exact(
            known
                .series(epsilon, Atom::Zero.as_view(), Rational::from(target))
                .map_err(failure)?
                .to_atom(),
        );
        self.check_coefficients(&result, settings, mode == RustRedMasterMode::Numerical)?;
        if mode == RustRedMasterMode::Numerical {
            result = self.resolve_numerical(result, settings)?;
        }
        Ok(result + retained)
    }

    fn master_multiplicities(
        value: AtomView,
        output: &mut BTreeMap<Atom, i64>,
    ) -> Result<(), RustRedEvaluationError> {
        match value {
            AtomView::Fun(function) if function.get_symbol() == vk_symbol!("RustRedMaster") => {
                *output.entry(value.to_owned()).or_default() += 1;
            }
            AtomView::Mul(product) => {
                for factor in product.iter() {
                    Self::master_multiplicities(factor, output)?;
                }
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let exponent = get_integer_from_atom(exponent)
                    .filter(|value| *value > 0)
                    .ok_or_else(|| {
                        failure("master expressions must be polynomials in terminal integrals")
                    })?;
                if !matches!(base, AtomView::Fun(function) if function.get_symbol() == vk_symbol!("RustRedMaster"))
                {
                    return Err(failure("unexpected factor in master monomial"));
                }
                *output.entry(base.to_owned()).or_default() += exponent;
            }
            AtomView::Num(_) if value == Atom::num(1).as_view() => {}
            _ => {
                return Err(failure(
                    "master expressions must be polynomials in terminal integrals",
                ));
            }
        }
        Ok(())
    }

    fn coefficient_references(expression: &Atom) -> BTreeSet<Atom> {
        let mut references = BTreeSet::new();
        expression.visitor(&mut |part| {
            if matches!(part, AtomView::Fun(function) if function.get_symbol() == vk_symbol!("RustRedCoefficient")) {
                references.insert(part.to_owned());
                false
            } else { true }
        });
        references
    }

    fn coefficient(
        &self,
        symbol: &Atom,
    ) -> Result<(&Master, &Coefficient, i64), RustRedEvaluationError> {
        let AtomView::Fun(function) = symbol.as_view() else {
            return Err(failure("malformed master coefficient"));
        };
        let arguments: Vec<_> = function.iter().collect();
        if arguments.len() != 3 || arguments[0] != self.coefficient_tag.as_view() {
            return Err(failure(format!(
                "coefficient {} belongs to another data snapshot; active version is {}",
                symbol.to_canonical_string(),
                self.version
            )));
        }
        let power = get_integer_from_atom(arguments[2])
            .ok_or_else(|| failure("malformed coefficient epsilon power"))?;
        let master = self
            .masters
            .get(&arguments[1].to_owned())
            .ok_or_else(|| failure(format!("unknown master {} in coefficient", arguments[1])))?;
        let coefficient = master.coefficients.get(&power).ok_or_else(|| RustRedEvaluationError::MasterCapability { detail: format!("master {} requires unavailable epsilon coefficient {power}; available powers {:?}; data {}", master.identity, master.coefficients.keys().collect::<Vec<_>>(), self.version) })?;
        Ok((master, coefficient, power))
    }

    fn finite_precision(expression: &Atom) -> Option<u32> {
        let mut bits = None;
        expression.visitor(&mut |part| {
            if let AtomView::Num(number) = part
                && let CoefficientView::Float(real, imaginary) = number.get_coeff_view()
            {
                for component in [real, imaginary] {
                    let value = component.to_float();
                    if !value.as_raw().is_zero() {
                        bits = Some(bits.map_or(value.prec(), |old: u32| old.min(value.prec())));
                    }
                }
            }
            true
        });
        bits.map(|value| (f64::from(value) * std::f64::consts::LOG10_2).floor() as u32)
    }

    fn check_coefficients(
        &self,
        expression: &Atom,
        settings: &VakintSettings,
        numerical: bool,
    ) -> Result<(), VakintError> {
        let mut issues = Vec::new();
        for symbol in Self::coefficient_references(expression) {
            match self.coefficient(&symbol) {
                Err(RustRedEvaluationError::MasterCapability { detail }) => issues.push(detail),
                Err(error) => return Err(error.into()),
                Ok((master, coefficient, power)) if numerical => {
                    let available = if let Some(exact) = &coefficient.exact {
                        if let Some(embedded) =
                            Self::finite_precision(&vk_parse!(exact).map_err(failure)?)
                        {
                            let validated = coefficient.significant_digits.filter(|digits| *digits > 0).ok_or_else(|| failure(format!("master {} coefficient epsilon^{power} embeds a finite constant without validated significant digits; data {}", master.identity, self.version)))?;
                            Some(embedded.min(validated))
                        } else {
                            None
                        }
                    } else {
                        coefficient.significant_digits
                    };
                    if available.is_some_and(|digits| settings.run_time_decimal_precision > digits)
                    {
                        issues.push(format!("master {} coefficient epsilon^{power}: requested {} significant decimal digits, available {}; data {}", master.identity, settings.run_time_decimal_precision, available.unwrap(), self.version));
                    }
                }
                _ => {}
            }
        }
        if issues.is_empty() {
            Ok(())
        } else {
            Err(RustRedEvaluationError::MasterCapability {
                detail: issues.join("; "),
            }
            .into())
        }
    }

    pub(super) fn resolve_numerical(
        &self,
        expression: Atom,
        settings: &VakintSettings,
    ) -> Result<Atom, VakintError> {
        let mut result = combine_exact(expression);
        let unresolved = self.unresolved(&result);
        if !unresolved.is_empty() {
            return Err(failure(format!(
                "Cannot numerically evaluate unresolved masters (data {}): {}",
                self.version,
                unresolved.join("; ")
            ))
            .into());
        }
        // Retained declared masters require the caller's loop-count-specific
        // Laurent expansion, not a hidden numerical evaluation convention.
        let retained = Self::referenced_masters(&result);
        if !retained.is_empty() {
            return Err(failure(format!(
                "Retained masters require symbolic or numerical master expansion first: {}",
                retained
                    .iter()
                    .map(Atom::to_canonical_string)
                    .collect::<Vec<_>>()
                    .join(", ")
            ))
            .into());
        }
        self.check_coefficients(&result, settings, true)?;
        // Period formulas can contain large cancellations. Keep guard digits
        // through materialization and subsequent symbolic coefficient sums;
        // capability checks above deliberately use the requested precision.
        let working_settings = VakintSettings {
            run_time_decimal_precision: settings
                .run_time_decimal_precision
                .checked_add(32)
                .ok_or_else(|| failure("numerical working precision overflow"))?,
            ..settings.clone()
        };
        let working_precision = working_settings.get_binary_precision();
        for symbol in Self::coefficient_references(&result) {
            let (master, coefficient, power) = self.coefficient(&symbol)?;
            let value = if let Some(exact) = &coefficient.exact {
                crate::hpl::normalize_constants(vk_parse!(exact).map_err(failure)?.as_view())?
            } else {
                let number = Float::parse(
                        coefficient.numerical.as_ref().ok_or_else(|| {
                            failure("coefficient has no numerical representation")
                        })?,
                        Some(working_precision),
                    )
                    .map_err(|error| failure(format!("master {} coefficient epsilon^{power} has an invalid numerical value: {error}; data {}", master.identity, self.version)))?;
                if !number.as_raw().is_finite() {
                    return Err(failure(format!("master {} coefficient epsilon^{power} has a non-finite numerical value; data {}", master.identity, self.version)).into());
                }
                Atom::num(number)
            };
            result = result.replace(symbol.to_pattern()).with(value);
        }
        // Tagged evaluatable constants deliberately contain metadata that the
        // generic `is_constant` predicate does not classify as numerical.
        // The shared bridge evaluates the whole HPL and never its word tags.
        result = crate::hpl::substitute_numerically(result.as_view(), &working_settings)?;
        let map = HashMap::<Atom, Complex<Float>>::new();
        let mut issue = None;
        result = result.replace_map(|part, _context, output| {
            // Keep numeric exponents and tensor indices exact. Evaluate whole
            // constant subexpressions before entering their children.
            if !matches!(part, AtomView::Num(_)) && part.is_constant() {
                match part.evaluate_with_prec(&map, working_precision) {
                    Ok(number) => **output = Atom::num(number),
                    Err(error) => issue = Some(error.to_string()),
                }
            }
        });
        if let Some(issue) = issue {
            return Err(failure(format!("master constant evaluation failed: {issue}")).into());
        }
        Ok(result)
    }
}

#[cfg(test)]
#[path = "masters_acceptance_tests.rs"]
mod acceptance_tests;

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::domains::float::Real;

    fn catalog(exact: Option<&str>, digits: u32) -> (MasterCatalog, Atom) {
        let mut catalog = MasterCatalog {
            version: "test-v1".into(),
            coefficient_tag: Atom::var(vk_symbol!("test_master_data")),
            masters: BTreeMap::new(),
            family_graphs: BTreeMap::new(),
            basis_aliases: BTreeMap::new(),
        };
        let symbol = catalog.symbol("test-family", &[1]);
        catalog.masters.insert(
            symbol.clone(),
            Master {
                identity: "test-master".into(),
                loop_count: 1,
                graph: "test-parent.dot".into(),
                coefficients: BTreeMap::from([
                    (
                        -1,
                        Coefficient {
                            exact: Some("1".into()),
                            numerical: Some("1".into()),
                            significant_digits: Some(digits),
                        },
                    ),
                    (
                        0,
                        Coefficient {
                            exact: exact.map(str::to_owned),
                            numerical: Some("2.5".into()),
                            significant_digits: Some(digits),
                        },
                    ),
                ]),
            },
        );
        (catalog, symbol)
    }

    fn settings(digits: u32, terms: i64) -> VakintSettings {
        VakintSettings {
            epsilon_symbol: "ep".into(),
            run_time_decimal_precision: digits,
            number_of_terms_in_epsilon_expansion: terms,
            ..Default::default()
        }
    }

    #[test]
    fn finite_coefficients_are_checked_only_when_needed() {
        let (catalog, master) = catalog(None, 30);
        let pole = catalog
            .finalize(
                master.clone(),
                &settings(100, 1),
                1,
                RustRedMasterMode::Numerical,
            )
            .unwrap();
        assert_eq!(pole, vk_parse!("ep^-1").unwrap());
        let error = catalog
            .finalize(master, &settings(31, 2), 1, RustRedMasterMode::Numerical)
            .unwrap_err()
            .to_string();
        assert!(error.contains("requested 31") && error.contains("available 30"));
    }

    #[test]
    fn symbolic_numerical_coefficients_resolve_with_active_data() {
        let (catalog, master) = catalog(None, 30);
        let symbolic = catalog
            .finalize(master, &settings(200, 2), 1, RustRedMasterMode::Symbolic)
            .unwrap();
        assert!(symbolic.contains_symbol(vk_symbol!("RustRedCoefficient")));
        assert!(
            catalog
                .resolve_numerical(symbolic.clone(), &settings(31, 2))
                .is_err()
        );
        let numerical = catalog
            .resolve_numerical(symbolic, &settings(30, 2))
            .unwrap();
        assert!(!numerical.contains_symbol(vk_symbol!("RustRedCoefficient")));
    }

    #[test]
    fn spurious_poles_require_deeper_master_coefficients() {
        let (catalog, master) = catalog(Some("5/2"), 30);
        let raw = master * vk_parse!("ep^-1").unwrap();
        let error = catalog
            .finalize(raw, &settings(200, 2), 1, RustRedMasterMode::Symbolic)
            .unwrap_err()
            .to_string();
        assert!(error.contains("unavailable epsilon coefficient 1"));
    }

    #[test]
    fn exact_cancellation_removes_unavailable_coefficients_and_unknowns() {
        let (catalog, master) = catalog(Some("5/2"), 30);
        let unknown = catalog.symbol("other-family", &[2]);
        let raw =
            (master.clone() - master + unknown.clone() - unknown) * vk_parse!("ep^-9").unwrap();
        assert!(
            catalog
                .finalize(raw, &settings(200, 2), 1, RustRedMasterMode::Numerical)
                .unwrap()
                .is_zero()
        );
    }

    #[test]
    fn missing_coefficients_are_not_assumed_zero() {
        let (catalog, master) = catalog(Some("5/2"), 30);
        assert!(
            catalog
                .finalize(
                    master.clone(),
                    &settings(200, 2),
                    1,
                    RustRedMasterMode::Symbolic
                )
                .is_ok()
        );
        assert!(
            catalog
                .finalize(master, &settings(200, 3), 1, RustRedMasterMode::Symbolic)
                .is_err()
        );
    }

    #[test]
    fn exact_leading_zero_avoids_unneeded_orders_in_master_products() {
        let (mut catalog, master) = catalog(Some("1"), 30);
        catalog
            .masters
            .get_mut(&master)
            .unwrap()
            .coefficients
            .get_mut(&-1)
            .unwrap()
            .exact = Some("0".into());
        let result = catalog
            .finalize(
                master.pow(Atom::num(2)),
                &settings(200, 3),
                2,
                RustRedMasterMode::Symbolic,
            )
            .unwrap();
        assert_eq!(result, Atom::num(1));
    }

    #[test]
    fn numerical_only_symbols_cannot_cross_data_snapshots() {
        let (mut catalog, master) = catalog(None, 30);
        let symbolic = catalog
            .finalize(master, &settings(200, 2), 1, RustRedMasterMode::Symbolic)
            .unwrap();
        catalog.coefficient_tag = Atom::var(vk_symbol!("replacement_data"));
        let error = catalog
            .resolve_numerical(symbolic, &settings(20, 2))
            .unwrap_err()
            .to_string();
        assert!(error.contains("another data snapshot"));
    }

    #[test]
    fn numerical_mode_evaluates_tagged_hpl_constants() {
        let (catalog, master) = catalog(Some("HI(1,0,om1)"), 30);
        let numerical = catalog
            .finalize(master, &settings(80, 2), 1, RustRedMasterMode::Numerical)
            .unwrap();
        assert!(!numerical.contains_symbol(vk_symbol!("symbolica::hpl_at_sixth_root")));
        let value = numerical
            .coefficient_list::<i32>(&[vk_parse!("ep").unwrap()])
            .into_iter()
            .find(|(power, _)| *power == Atom::num(1))
            .unwrap()
            .1;
        assert!(matches!(value.as_view(), AtomView::Num(_)));
        let number = value
            .evaluate_with_prec(&HashMap::<Atom, Complex<Float>>::new(), 256)
            .unwrap();
        assert!(number.re < Float::with_val(256, 0));
    }

    #[test]
    fn malformed_or_nonfinite_numerical_coefficients_report_the_master() {
        for literal in ["NaN", "inf", "not-a-number"] {
            let (mut catalog, master) = catalog(None, 30);
            catalog
                .masters
                .get_mut(&master)
                .unwrap()
                .coefficients
                .get_mut(&0)
                .unwrap()
                .numerical = Some(literal.into());
            let error = catalog
                .finalize(master, &settings(20, 2), 1, RustRedMasterMode::Numerical)
                .unwrap_err()
                .to_string();
            assert!(
                error.contains("test-master") && error.contains("epsilon^0"),
                "{error}"
            );
        }
    }

    #[test]
    fn embedded_float_annotations_do_not_override_validated_precision() {
        let (mut catalog, master) = catalog(Some("2.5`100"), 5);
        let error = catalog
            .finalize(
                master.clone(),
                &settings(20, 2),
                1,
                RustRedMasterMode::Numerical,
            )
            .unwrap_err();
        assert!(error.to_string().contains("available 5"), "{error}");
        catalog
            .masters
            .get_mut(&master)
            .unwrap()
            .coefficients
            .get_mut(&0)
            .unwrap()
            .significant_digits = None;
        let error = catalog
            .finalize(master, &settings(20, 2), 1, RustRedMasterMode::Numerical)
            .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("without validated significant digits"),
            "{error}"
        );
    }

    #[test]
    fn exact_rational_and_zeta_coefficients_are_not_capped_by_numeric_metadata() {
        for formula in ["5/2", "symbolica::zeta(3)"] {
            let (catalog, master) = catalog(Some(formula), 5);
            let settings = settings(80, 2);
            let numerical = catalog
                .finalize(master, &settings, 1, RustRedMasterMode::Numerical)
                .unwrap();
            let coefficient = numerical
                .coefficient_list::<i32>(&[vk_parse!("ep").unwrap()])
                .into_iter()
                .find(|(power, _)| *power == Atom::num(1))
                .unwrap()
                .1;
            let parameters = HashMap::<Atom, Complex<Float>>::new();
            let actual = coefficient
                .evaluate_with_prec(&parameters, settings.get_binary_precision())
                .unwrap();
            let expected = vk_parse!(formula)
                .unwrap()
                .evaluate_with_prec(&parameters, settings.get_binary_precision())
                .unwrap();
            let tolerance = Float::parse("1e-70", Some(settings.get_binary_precision())).unwrap();
            assert!((actual.re - expected.re).norm() < tolerance);
        }
    }

    #[test]
    fn finite_input_coefficients_preserve_constant_powers_and_laurent_orders() {
        let (catalog, master) = catalog(Some("2"), 80);
        let settings = settings(80, 2);
        let parameters = HashMap::<Atom, Complex<Float>>::new();
        let input = vk_parse!("1/(4608*𝜋^2)+0.0009101213175016047579806603153986/𝜋^4").unwrap();
        let expected = input.evaluate_with_prec(&parameters, 256).unwrap();
        let actual = catalog
            .resolve_numerical(input, &settings)
            .unwrap()
            .evaluate_with_prec(&parameters, 256)
            .unwrap();
        let tolerance = Float::parse("1e-60", Some(256)).unwrap();
        assert!((actual.re - expected.re).norm() < tolerance);

        let raw = vk_parse!("0.125/(1-ep)").unwrap() * &master;
        let result = catalog
            .finalize(raw, &settings, 1, RustRedMasterMode::Numerical)
            .unwrap();
        let coefficients = result.coefficient_list::<i32>(&[vk_parse!("ep").unwrap()]);
        assert_eq!(coefficients.len(), 2);
        for (power, coefficient) in coefficients {
            let expected = if power == vk_parse!("ep^-1").unwrap() {
                Float::with_val(256, 0.125)
            } else {
                assert_eq!(power, Atom::num(1));
                Float::with_val(256, 0.375)
            };
            let actual = coefficient.evaluate_with_prec(&parameters, 256).unwrap();
            assert!((actual.re - expected).norm() < tolerance);
        }
        let equivalent = catalog
            .finalize(
                vk_parse!("0.125*(1/(1-ep)-1-ep)/ep^2").unwrap() * &master,
                &settings,
                1,
                RustRedMasterMode::Numerical,
            )
            .unwrap();
        assert_eq!(equivalent, result);
        let zero = catalog
            .finalize(
                vk_parse!("0.125*(1/(ep*(ep-1))-1/(ep-1)+1/ep)").unwrap() * master,
                &settings,
                1,
                RustRedMasterMode::Numerical,
            )
            .unwrap();
        assert!(zero.is_zero());
    }

    #[test]
    fn declared_basis_relations_cancel_missing_orders_without_discarding_richer_alias_data() {
        let (mut catalog, master) = catalog(Some("2"), 80);
        let alias = catalog.symbol("test-family", &[2]);
        catalog.masters.insert(
            alias.clone(),
            Master {
                identity: "test-alias".into(),
                loop_count: 1,
                graph: "test-parent.dot".into(),
                coefficients: BTreeMap::from([
                    (
                        -1,
                        Coefficient {
                            exact: Some("3".into()),
                            numerical: Some("3".into()),
                            significant_digits: Some(80),
                        },
                    ),
                    (
                        0,
                        Coefficient {
                            exact: Some("2".into()),
                            numerical: Some("2".into()),
                            significant_digits: Some(80),
                        },
                    ),
                ]),
            },
        );
        catalog.load_basis_aliases(vec![
            (1, master.clone(), vk_parse!("test_basis").unwrap()),
            (1, alias.clone(), vk_parse!("(3-4*ep)*test_basis").unwrap()),
        ]);
        let raw =
            (alias.clone() - vk_parse!("3-4*ep").unwrap() * master) / vk_parse!("ep").unwrap();
        let retained = catalog
            .finalize(raw.clone(), &settings(80, 2), 1, RustRedMasterMode::Retain)
            .unwrap();
        assert!(MasterCatalog::referenced_masters(&retained).contains(&alias));
        assert!(
            catalog
                .finalize(raw, &settings(80, 2), 1, RustRedMasterMode::Numerical)
                .unwrap()
                .is_zero()
        );

        catalog
            .masters
            .get_mut(&alias)
            .unwrap()
            .coefficients
            .insert(
                1,
                Coefficient {
                    exact: Some("9".into()),
                    numerical: Some("9".into()),
                    significant_digits: Some(80),
                },
            );
        let result = catalog
            .finalize(
                alias.clone(),
                &settings(80, 3),
                1,
                RustRedMasterMode::Symbolic,
            )
            .unwrap();
        assert!(
            (result - vk_parse!("3/ep+2+9*ep").unwrap())
                .expand()
                .is_zero()
        );
        let corrupted = catalog
            .masters
            .get_mut(&alias)
            .unwrap()
            .coefficients
            .get_mut(&0)
            .unwrap();
        corrupted.exact = None;
        corrupted.numerical = Some("NaN".into());
        let error = catalog
            .finalize(alias, &settings(80, 2), 1, RustRedMasterMode::Numerical)
            .unwrap_err();
        assert!(matches!(
            error,
            VakintError::RustRedEvaluation(RustRedEvaluationError::ExternalData { .. })
        ));
        assert!(error.to_string().contains("test-alias"));
    }

    #[test]
    fn identical_function_basis_expressions_cancel_missing_alias_orders() {
        let (mut catalog, master) = catalog(Some("2"), 80);
        let (duplicate_catalog, _) = self::catalog(Some("2"), 80);
        let duplicate = catalog.symbol("test-family", &[3]);
        catalog.masters.insert(
            duplicate.clone(),
            Master {
                identity: "duplicate-function-master".into(),
                ..duplicate_catalog.masters.into_values().next().unwrap()
            },
        );
        let basis = vk_parse!("Gam(1,1)/(ep*(1-ep))").unwrap();
        catalog.load_basis_aliases(vec![
            (1, master.clone(), basis.clone()),
            (1, duplicate.clone(), basis),
        ]);
        let result = catalog
            .finalize(
                (master - duplicate) / vk_parse!("ep").unwrap(),
                &settings(80, 2),
                1,
                RustRedMasterMode::Numerical,
            )
            .unwrap();
        assert!(result.is_zero());
    }

    #[test]
    fn vendored_fallback_graphs_keep_their_package_version() {
        let vendored = super::super::data::RustRedData::vendored().unwrap();
        let package = Arc::clone(&vendored.loops[&3]);
        let catalog = MasterCatalog::load(
            &BTreeMap::from([(3, Arc::clone(&package))]),
            "local-mixed-release",
        )
        .unwrap();
        assert_ne!(catalog.version, package.manifest.version);
        let expected = format!(
            "https://github.com/alphal00p/vakint-external-data/blob/{}/3L/",
            package.manifest.version
        );
        assert!(!catalog.family_graphs.is_empty());
        assert!(!catalog.masters.is_empty());
        assert!(
            catalog
                .family_graphs
                .values()
                .all(|graph| graph.starts_with(&expected))
        );
        assert!(
            catalog
                .masters
                .values()
                .all(|master| master.graph.starts_with(&expected))
        );
    }

    #[test]
    fn exact_alias_cancellation_precedes_amplified_numerical_rounding() {
        let (mut catalog, master) = catalog(None, 50);
        catalog
            .masters
            .get_mut(&master)
            .unwrap()
            .coefficients
            .get_mut(&0)
            .unwrap()
            .numerical = Some("0.33333333333333333333333333333333333333333333333333".into());
        let (duplicate_catalog, _) = self::catalog(None, 50);
        let alias = catalog.symbol("test-family", &[4]);
        let mut alias_data = duplicate_catalog.masters.into_values().next().unwrap();
        alias_data.identity = "five-times-master".into();
        alias_data.coefficients.get_mut(&-1).unwrap().exact = Some("5".into());
        alias_data.coefficients.get_mut(&0).unwrap().numerical =
            Some("1.6666666666666666666666666666666666666666666666667".into());
        catalog.masters.insert(alias.clone(), alias_data);
        catalog.load_basis_aliases(vec![
            (1, master.clone(), vk_parse!("test_basis").unwrap()),
            (1, alias.clone(), vk_parse!("5*test_basis").unwrap()),
        ]);
        let raw = vk_parse!("10^60").unwrap() * (alias - Atom::num(5) * master);
        assert!(
            !catalog
                .finalize_independent(
                    raw.clone(),
                    &settings(32, 2),
                    1,
                    RustRedMasterMode::Numerical
                )
                .unwrap()
                .is_zero()
        );
        assert!(
            catalog
                .finalize(raw, &settings(32, 2), 1, RustRedMasterMode::Numerical)
                .unwrap()
                .is_zero()
        );
    }
}
