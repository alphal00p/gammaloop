//! Runtime master coefficients, with availability checked after exact cancellation.
//!
//! Input masters are unit-mass integrals in the package's FMFT/MATAD convention.
//! The caller restores target mass powers and measure normalization exactly once.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use serde::Deserialize;
use sha2::{Digest, Sha256};
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::coefficient::CoefficientView;
use symbolica::domains::rational::Rational;
use symbolica::id::Replacement;
use symbolica::poly::series::SeriesDepth;

use super::data::{FileReference, LoopData, checked_bytes, failure};
use super::{
    CandidateRoutedApplicationObserver, MonitorContext, RustRedEvaluationError, RustRedMasterMode,
};
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{VakintError, VakintSettings, get_integer_from_atom};

#[derive(Deserialize)]
struct CoefficientMetadata {
    epsilon_power: i64,
    symbolic: bool,
    numerical: bool,
    significant_digits: Option<u32>,
    #[serde(default)]
    absolute_error: Option<String>,
    #[serde(default)]
    error_quality: Option<String>,
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
pub(crate) struct Coefficient {
    pub(crate) exact: Option<String>,
    pub(crate) numerical: Option<String>,
    pub(crate) significant_digits: Option<u32>,
    pub(crate) absolute_error: Option<String>,
    pub(crate) error_is_bounded: bool,
}

#[derive(Debug)]
pub(crate) struct Master {
    pub(crate) identity: String,
    pub(crate) loop_count: usize,
    pub(crate) coefficients: BTreeMap<i64, Coefficient>,
    graph: String,
    pub(crate) metadata: serde_json::Value,
    pub(crate) data: Option<Arc<LoopData>>,
}

#[derive(Debug)]
pub(crate) struct MasterCatalog {
    pub(crate) version: String,
    coefficient_tag: Atom,
    pub(crate) masters: BTreeMap<Atom, Master>,
    pub(crate) constants: BTreeMap<Atom, crate::inspection::FiniteSource>,
    family_graphs: BTreeMap<Atom, String>,
    basis_aliases: BTreeMap<Atom, Atom>,
}

/// A declaration read from a native application's authenticated terminal inventory.
/// It supplies declarations only, never numerical coefficients or invented graphs.
#[derive(Clone, Debug)]
pub(crate) struct BareMasterDeclaration {
    pub(crate) application_path: std::path::PathBuf,
}

pub(crate) type BareMasterDeclarations = BTreeMap<Atom, BareMasterDeclaration>;

pub(crate) fn combine_exact(expression: Atom) -> Atom {
    let expression = expression.expand();
    // Symbolica's rational canonicalization requires exact coefficients.
    // In particular, `together` currently corrupts powers in sums with
    // floating coefficients, so numerical expressions stay expanded.
    if MasterCatalog::has_floating_coefficients(&expression) {
        expression
    } else {
        expression.together()
    }
}

impl MasterCatalog {
    fn has_floating_coefficients(expression: &Atom) -> bool {
        let mut has_float = false;
        expression.visitor(&mut |part| {
            if let AtomView::Num(number) = part
                && matches!(number.get_coeff_view(), CoefficientView::Float(_, _))
            {
                has_float = true;
            }
            !has_float
        });
        has_float
    }

    #[cfg(test)]
    pub(crate) fn test_catalog(exact: Option<&str>, digits: u32) -> (Self, Atom) {
        tests::catalog(exact, digits)
    }

    pub(super) fn load(
        loops: &BTreeMap<usize, Arc<LoopData>>,
        version: &str,
    ) -> Result<Self, RustRedEvaluationError> {
        crate::Vakint::initialize_vakint_symbols();
        let mut catalog = Self {
            version: version.into(),
            coefficient_tag: Atom::Zero,
            masters: BTreeMap::new(),
            constants: BTreeMap::new(),
            family_graphs: BTreeMap::new(),
            basis_aliases: BTreeMap::new(),
        };
        let mut basis_expressions = Vec::new();
        let mut identity = Sha256::new();
        identity.update(version.as_bytes());
        for (&count, data) in loops {
            catalog.load_constants(data, &mut identity)?;
            let mut fingerprints = BTreeMap::new();
            let mut parent_graphs = BTreeMap::<Atom, BTreeSet<String>>::new();
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
                    let error_is_bounded = match entry.error_quality.as_deref() {
                        None | Some("estimated") => false,
                        Some("bounded") if entry.absolute_error.is_some() => true,
                        Some("bounded") => {
                            return Err(failure("bounded coefficient lacks an absolute error"));
                        }
                        Some(quality) => {
                            return Err(failure(format!(
                                "unknown coefficient error quality {quality}"
                            )));
                        }
                    };
                    if let Some(bound) = &entry.absolute_error {
                        let bound = symbolica::domains::float::Float::parse(bound, Some(64))
                            .map_err(failure)?;
                        if !bound.as_raw().is_finite()
                            || (bound.as_raw().is_sign_negative() && !bound.as_raw().is_zero())
                        {
                            return Err(failure(
                                "coefficient absolute error must be finite and nonnegative",
                            ));
                        }
                    }
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
                                absolute_error: entry.absolute_error,
                                error_is_bounded,
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
                    metadata: value.clone(),
                    data: Some(Arc::clone(data)),
                };
                parent_graphs
                    .entry(Self::family_symbol(fingerprint))
                    .or_default()
                    .insert(master.graph.clone());
                if catalog.masters.insert(symbol, master).is_some() {
                    return Err(failure("duplicate family-qualified master identity"));
                }
            }
            // A common denominator family can cover several maximal parents.
            // Prefer those authenticated references to its representative graph;
            // the family graph remains the fallback when no parents are listed.
            for (family, graphs) in parent_graphs {
                catalog
                    .family_graphs
                    .insert(family, graphs.into_iter().collect::<Vec<_>>().join("; "));
            }
        }
        if basis_expressions
            .iter()
            .any(|(_, _, basis)| catalog.constants.contains_key(basis))
        {
            return Err(failure(
                "a finite constant shadows a declared master basis identity",
            ));
        }
        let tag = format!("rustred_data_{:x}", identity.finalize());
        catalog.coefficient_tag = Atom::var(vk_symbol!(&tag));
        catalog.load_basis_aliases(basis_expressions);
        Ok(catalog)
    }

    fn load_constants(
        &mut self,
        data: &LoopData,
        identity: &mut Sha256,
    ) -> Result<(), RustRedEvaluationError> {
        #[derive(Deserialize)]
        struct Constant {
            symbol: String,
            value: String,
            significant_digits: u32,
            absolute_error: String,
            error_quality: String,
            definition: String,
            provenance: serde_json::Value,
        }
        let Some(file) = &data.manifest.constant_bindings else {
            return Ok(());
        };
        let bytes = data.read(&file.path)?;
        checked_bytes(&bytes, &file.sha256)?;
        identity.update(&bytes);
        let entries: Vec<serde_json::Value> = serde_json::from_slice(&bytes).map_err(failure)?;
        for entry in entries {
            let constant: Constant = serde_json::from_value(entry.clone()).map_err(failure)?;
            let symbol = vk_parse!(&constant.symbol).map_err(failure)?;
            let reserved_head = symbol.get_symbol().is_some_and(|head| {
                head.is_builtin()
                    || head.get_namespace() == "symbolica"
                    || head.get_evaluation_info().is_some()
                    || crate::symbols::SYMBOL_REGISTRY.contains(&head)
                    || head == vk_symbol!("NumericalSource")
                    || (head.get_namespace() == crate::NAMESPACE
                        && head.get_stripped_name().starts_with("RustRed"))
            });
            if !matches!(symbol.as_view(), AtomView::Var(_) | AtomView::Fun(_))
                || Self::finite_precision(&symbol).is_some()
                || reserved_head
                || Self::is_recorded_period(symbol.as_view())
                || crate::inspection::constants::builtin_definition(symbol.as_view())
                    .map_err(failure)?
                    .is_some()
            {
                return Err(failure(format!(
                    "invalid finite constant identity {}",
                    constant.symbol
                )));
            }
            if constant.significant_digits == 0
                || constant.definition.trim().is_empty()
                || constant
                    .provenance
                    .as_object()
                    .is_none_or(|value| value.is_empty())
            {
                return Err(failure(
                    "finite constant lacks precision, definition or provenance",
                ));
            }
            let value = symbolica::domains::float::Float::parse(&constant.value, Some(64))
                .map_err(failure)?;
            let bound = symbolica::domains::float::Float::parse(&constant.absolute_error, Some(64))
                .map_err(failure)?;
            if !value.as_raw().is_finite()
                || !bound.as_raw().is_finite()
                || (bound.as_raw().is_sign_negative() && !bound.as_raw().is_zero())
            {
                return Err(failure(
                    "finite constant value/error must be finite and error nonnegative",
                ));
            }
            let bounded = match constant.error_quality.as_str() {
                "bounded" => true,
                "estimated" => false,
                _ => return Err(failure("unknown finite constant error quality")),
            };
            let source = crate::inspection::FiniteSource {
                value: constant.value,
                significant_digits: constant.significant_digits,
                absolute_error: Some(constant.absolute_error),
                error_is_bounded: bounded,
                identity: symbol.to_canonical_string(),
                master_id: None,
                epsilon_power: None,
                data_version: data.manifest.version.clone(),
                source_sha256: format!(
                    "{:x}",
                    Sha256::digest(serde_json::to_vec(&entry).map_err(failure)?)
                ),
            };
            if self.constants.insert(symbol, source).is_some() {
                return Err(failure(
                    "duplicate finite constant identity across loop packs",
                ));
            }
        }
        Ok(())
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

    pub(crate) fn coefficient_symbol(&self, master: &Atom, power: i64) -> Atom {
        vk_symbol!("RustRedCoefficient").call((&self.coefficient_tag, master, power))
    }

    pub(crate) fn coefficient_snapshot_tag(&self) -> &Atom {
        &self.coefficient_tag
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

    fn combine_master_coefficients(raw: Atom) -> Atom {
        let raw = raw.expand();
        // Keep the same whole-expression Float policy as combine_exact.
        if Self::has_floating_coefficients(&raw) {
            return raw;
        }
        let symbols: Vec<_> = Self::referenced_masters(&raw).into_iter().collect();
        let terms = raw.coefficient_list::<i32>(&symbols);
        // Collection leaves opaque functions and composite denominators in the
        // coefficients. Master-dependent coefficients are not independent, so
        // their cancellations still require general rational combination.
        if terms
            .iter()
            .any(|(_, coefficient)| !Self::referenced_masters(coefficient).is_empty())
        {
            return combine_exact(raw);
        }
        // Terminal products are independent indeterminates. Combine their
        // scalar coefficients before classifying any surviving unknowns,
        // without rationalizing the whole polynomial in those terminals.
        Atom::add_many(
            terms
                .into_iter()
                .map(|(monomial, coefficient)| monomial * combine_exact(coefficient)),
        )
    }

    fn unresolved(
        &self,
        expression: &Atom,
        require_coefficients: bool,
        announced: &BareMasterDeclarations,
    ) -> Vec<String> {
        let mut issues: Vec<_> = Self::referenced_masters(expression)
            .iter()
            .filter_map(|symbol| match self.masters.get(symbol) {
                Some(master) if !require_coefficients || !master.coefficients.is_empty() => None,
                Some(master) => Some(format!(
                    "{}: declared terminal has no computed coefficients; parent graph {}",
                    master.identity, master.graph
                )),
                None if announced.contains_key(symbol) => require_coefficients.then(|| format!(
                    "{}: declared terminal has no computed coefficients; native application {}; load its matching external-data pack for master coefficients and graphs",
                    symbol.to_canonical_string(), announced[symbol].application_path.display()
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
                        "{}: undeclared residual master; candidate parent graph DOT references (no unique parent inferred): {}",
                        symbol.to_canonical_string(),
                        graph
                            .map(String::as_str)
                            .unwrap_or("not present in this data snapshot")
                    ))
                }
            })
            .collect();
        issues.sort();
        issues
    }

    pub(super) fn retain(
        &self,
        raw: Atom,
        announced: &BareMasterDeclarations,
        mut observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
    ) -> Atom {
        let started = std::time::Instant::now();
        MonitorContext::native_phase(
            observation.as_deref_mut(),
            "retained_catalogue_exact_combination",
            Some(raw.as_view()),
        );
        let raw = Self::combine_master_coefficients(raw);
        log::debug!(target: "vakint::saved_application", "phase=retained_catalogue_exact_combination elapsed_seconds={:.6}", started.elapsed().as_secs_f64());
        let started = std::time::Instant::now();
        MonitorContext::native_phase(
            observation,
            "retained_catalogue_classification",
            Some(raw.as_view()),
        );
        let unresolved = self.unresolved(&raw, false, announced);
        if !unresolved.is_empty() {
            log::warn!(
                "Unresolved RustRed masters (data {}):\n  - {}",
                self.version,
                unresolved.join("\n  - ")
            );
        }
        log::debug!(target: "vakint::saved_application", "phase=retained_catalogue_classification elapsed_seconds={:.6}", started.elapsed().as_secs_f64());
        raw
    }

    pub(super) fn require_coefficients(
        &self,
        expression: &Atom,
        announced: &BareMasterDeclarations,
    ) -> Result<(), VakintError> {
        let unresolved = self.unresolved(expression, true, announced);
        if unresolved.is_empty() {
            Ok(())
        } else {
            Err(RustRedEvaluationError::MasterCapability {
                detail: format!(
                    "Unresolved RustRed masters (data {}):\n  - {}",
                    self.version,
                    unresolved.join("\n  - ")
                ),
            }
            .into())
        }
    }

    pub(super) fn finalize(
        &self,
        raw: Atom,
        settings: &VakintSettings,
        loop_count: usize,
        mode: RustRedMasterMode,
    ) -> Result<Atom, VakintError> {
        // User normalization factors can introduce finite literals after the
        // frontend staged the original input. Bind those remaining literals
        // before rational algebra; an exact zero prefactor must not launch an
        // unbounded search for its first nonzero Laurent coefficient.
        let raw = if mode == RustRedMasterMode::Numerical {
            crate::precision::legacy::stage_table(
                vk_symbol!("rustred_input")
                    .call(loop_count as i64)
                    .as_view(),
                raw.as_view(),
                vk_symbol!(&settings.epsilon_symbol),
                "Input",
            )?
        } else {
            raw
        };
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
        if mode == RustRedMasterMode::Retain {
            return Ok(self.retain(raw, &BTreeMap::new(), None));
        }
        let raw = combine_exact(raw);
        self.require_coefficients(&raw, &BTreeMap::new())?;
        let epsilon = vk_symbol!(&settings.epsilon_symbol);
        let target = settings
            .number_of_terms_in_epsilon_expansion
            .checked_sub(i64::try_from(loop_count).map_err(failure)?)
            .and_then(|value| value.checked_sub(1))
            .ok_or_else(|| failure("epsilon expansion depth overflow"))?;
        if !self.basis_aliases.contains_key(&raw)
            && let Some(master) = self.masters.get(&raw)
        {
            // Aliases retain their original validation order, so fallback
            // cannot hide malformed lower-order coefficient formulas.
            // A bare master needs its highest requested coefficient directly;
            // no other term or finite coefficient can cancel this dependency.
            // General expressions retain the exact combination below.
            let first = -i64::try_from(master.loop_count).map_err(failure)?;
            if target >= first {
                self.coefficient(&self.coefficient_symbol(&raw, target))?;
            }
        }
        let symbols: Vec<_> = Self::referenced_masters(&raw).into_iter().collect();
        let terms = raw.coefficient_list::<i32>(&symbols);
        let mut required = BTreeMap::<Atom, i64>::new();
        let mut known = Atom::Zero;
        for (monomial, coefficient) in terms {
            let coefficient = combine_exact(coefficient);
            if coefficient.is_zero() {
                continue;
            }
            let mut multiplicities = BTreeMap::new();
            Self::master_multiplicities(monomial.as_view(), &mut multiplicities)?;
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
                        let exact = crate::hpl::normalize_constants(exact.as_view())?;
                        // Numerical evaluation needs the original coefficient
                        // identity for precision diagnostics. Exact zeros must
                        // still disappear before epsilon availability checks.
                        exact_values.insert(placeholder.clone(), exact);
                    }
                }
                expansion += placeholder * Atom::var(epsilon).pow(Atom::num(power));
            }
            known = known.replace(symbol.to_pattern()).with(expansion);
        }
        // Exact zero coefficients and exact relations must cancel before an
        // absent/high-order or low-precision coefficient can cause an error.
        for (source, value) in &exact_values {
            if mode == RustRedMasterMode::Symbolic || value.is_zero() {
                known = known.replace(source.to_pattern()).with(value.clone());
            }
        }
        let mut result = combine_exact(
            known
                .series(epsilon, Atom::Zero.as_view(), Rational::from(target))
                .map_err(failure)?
                .to_atom(),
        );
        // Keeping exact coefficient identities improves diagnostics, but an
        // unavailable order can multiply a provably zero combination of known
        // lower orders. Resolve that algebraic question before rejecting it.
        if mode == RustRedMasterMode::Numerical
            && matches!(
                self.check_coefficients(&result, false),
                Err(VakintError::RustRedEvaluation(
                    RustRedEvaluationError::MasterCapability { .. }
                ))
            )
        {
            let replacements = exact_values
                .into_iter()
                .map(|(source, value)| Replacement::new(source, value))
                .collect::<Vec<_>>();
            result = combine_exact(result.replace_multiple(&replacements));
        }
        self.check_coefficients(&result, mode == RustRedMasterMode::Numerical)?;
        if mode == RustRedMasterMode::Numerical {
            result = self.resolve_numerical(result, settings)?;
        }
        if mode == RustRedMasterMode::Symbolic {
            // An exact formula can still depend on snapshot-local finite
            // periods. Bind those leaves before returning a symbolic result,
            // so loading another release cannot silently change its value.
            return self.bind_period_references(result);
        }
        Ok(result)
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

    pub(crate) fn coefficient(
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

    pub(crate) fn finite_precision(expression: &Atom) -> Option<u32> {
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

    fn check_coefficients(&self, expression: &Atom, numerical: bool) -> Result<(), VakintError> {
        let mut issues = Vec::new();
        for symbol in Self::coefficient_references(expression) {
            match self.coefficient(&symbol) {
                Err(RustRedEvaluationError::MasterCapability { detail }) => issues.push(detail),
                Err(error) => return Err(error.into()),
                Ok((master, coefficient, power)) if numerical => {
                    if let Some(exact) = &coefficient.exact {
                        if Self::finite_precision(&vk_parse!(exact).map_err(failure)?).is_some()
                            && coefficient
                                .significant_digits
                                .is_none_or(|digits| digits == 0)
                        {
                            return Err(failure(format!("master {} coefficient epsilon^{power} embeds a finite constant without validated significant digits; data {}", master.identity, self.version)).into());
                        }
                    } else {
                        Self::validate_decimal(
                            coefficient.numerical.as_deref().ok_or_else(|| {
                                failure("coefficient has no numerical representation")
                            })?,
                            &master.identity,
                            power,
                        )?;
                    }
                    // Source limits are reported only by tracked evaluation.
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

    /// Validate numerical intent without rounding any source before final evaluation.
    pub(crate) fn resolve_numerical(
        &self,
        expression: Atom,
        _settings: &VakintSettings,
    ) -> Result<Atom, VakintError> {
        let result = combine_exact(expression);
        let unresolved = self.unresolved(&result, true, &BTreeMap::new());
        if !unresolved.is_empty() {
            return Err(RustRedEvaluationError::MasterCapability {
                detail: format!(
                    "Cannot numerically evaluate unresolved masters (data {}):\n  - {}",
                    self.version,
                    unresolved.join("\n  - ")
                ),
            }
            .into());
        }
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
        self.check_coefficients(&result, true)?;
        // Values are resolved by prepare_evaluation into independently tracked
        // leaves when the complete Laurent coefficient is evaluated.
        Ok(result)
    }
}

#[cfg(test)]
#[path = "masters_acceptance_tests.rs"]
mod acceptance_tests;

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use symbolica::domains::float::{Complex, Float, Real};

    pub(super) fn catalog(exact: Option<&str>, digits: u32) -> (MasterCatalog, Atom) {
        crate::Vakint::initialize_vakint_symbols();
        let mut catalog = MasterCatalog {
            version: "test-v1".into(),
            coefficient_tag: Atom::var(vk_symbol!("test_master_data")),
            masters: BTreeMap::new(),
            constants: BTreeMap::new(),
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
                metadata: serde_json::Value::Null,
                data: None,
                coefficients: BTreeMap::from([
                    (
                        -1,
                        Coefficient {
                            exact: Some("1".into()),
                            numerical: Some("1".into()),
                            significant_digits: Some(digits),
                            absolute_error: None,
                            error_is_bounded: false,
                        },
                    ),
                    (
                        0,
                        Coefficient {
                            exact: exact.map(str::to_owned),
                            numerical: Some("2.5".into()),
                            significant_digits: Some(digits),
                            absolute_error: None,
                            error_is_bounded: false,
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
    fn retained_products_combine_coefficients_before_unknown_classification() {
        let (catalog, master) = catalog(Some("5/2"), 30);
        let unknown = catalog.symbol("unannounced-product", &[2]);
        let announced = catalog.symbol("announced-product", &[1, 1]);
        let declarations = BTreeMap::from([(
            announced.clone(),
            BareMasterDeclaration {
                application_path: "native/application.json".into(),
            },
        )]);
        let product = master.pow(Atom::num(2)) * &announced;
        let coefficient = vk_parse!("(1-ep)^-1-1-ep").unwrap();
        let rational_zero =
            &coefficient * vk_parse!("ep^-2").unwrap() - vk_parse!("(1-ep)^-1").unwrap();
        let spectator = vk_parse!("mUV^-6*mursq^5*coupling^10*𝑒^(-5*𝛾*ep)*dot(p,p)").unwrap();
        let raw = &product * &spectator * coefficient
            + &product * &spectator * vk_parse!("1+ep").unwrap()
            + rational_zero * unknown.pow(Atom::num(2)) * &master;
        let retained = catalog.retain(raw, &declarations, None);
        assert_eq!(
            catalog.retain(retained.clone(), &declarations, None),
            retained
        );
        let symbols = vec![master, announced];
        let terms = retained.coefficient_list::<i32>(&symbols);
        assert_eq!(terms.len(), 1);
        assert_eq!(terms[0].0, product);
        assert!(combine_exact(&terms[0].1 - spectator * vk_parse!("(1-ep)^-1").unwrap()).is_zero());
        assert!(
            catalog
                .unresolved(&retained, false, &declarations)
                .is_empty()
        );
        assert!(!catalog.masters.contains_key(&symbols[1]));
    }

    #[test]
    fn retained_nonpolynomial_master_dependence_keeps_general_cancellation() {
        let (catalog, master) = catalog(Some("5/2"), 30);
        let unknown = catalog.symbol("nonpolynomial-unknown", &[2]);
        let raw = &master
            + (unknown.pow(Atom::num(2)) + Atom::num(2) * &unknown + Atom::num(1))
                / (&unknown + Atom::num(1))
            - (&unknown + Atom::num(1));
        assert_eq!(catalog.retain(raw, &BTreeMap::new(), None), master);

        let opaque = vk_symbol!("opaque_master_function").call(&unknown);
        assert_eq!(
            catalog.retain(&opaque + &master, &BTreeMap::new(), None),
            combine_exact(opaque + master)
        );
    }

    #[test]
    fn retained_float_input_keeps_the_complete_expanded_representation() {
        let (catalog, master) = catalog(Some("5/2"), 30);
        let unknown = catalog.symbol("float-policy-unknown", &[2]);
        let zero = vk_parse!("ep^-2*((1-ep)^-1-1-ep)-(1-ep)^-1").unwrap();
        let raw = vk_parse!("(0.125+ep)/(1-ep)").unwrap() * master + zero * unknown;
        let expected = raw.expand();
        assert!(MasterCatalog::has_floating_coefficients(&expected));
        let retained = catalog.retain(raw, &BTreeMap::new(), None);
        assert_eq!(retained.as_view().get_data(), expected.as_view().get_data());
        assert_eq!(
            MasterCatalog::finite_precision(&retained),
            MasterCatalog::finite_precision(&expected)
        );
    }

    #[test]
    fn retained_coefficients_preserve_snapshot_and_independent_source_identities() {
        let (catalog, master) = catalog(Some("5/2"), 30);
        let coefficient = catalog.coefficient_symbol(&master, 0);
        let first = vk_parse!("NumericalSource(input_a,0,0,3/10,128)").unwrap();
        let second = vk_parse!("NumericalSource(input_b,0,0,3/10,128)").unwrap();
        let zero = vk_parse!("(1-ep)^-1-1-ep*(1-ep)^-1").unwrap();
        let raw =
            &master * (&coefficient + &first - &second) + master.pow(Atom::num(2)) * zero * &first;
        let retained = catalog.retain(raw, &BTreeMap::new(), None);
        let terms = retained.coefficient_list::<i32>(std::slice::from_ref(&master));
        assert_eq!(terms.len(), 1);
        assert_eq!(terms[0].0, master);
        assert!(!terms[0].1.is_zero());
        assert!(combine_exact(&terms[0].1 - coefficient - first + second).is_zero());
        assert!(MasterCatalog::referenced_masters(&terms[0].1).is_empty());
    }

    #[test]
    fn retained_scalar_and_zero_expressions_keep_exact_combination() {
        let (catalog, _) = catalog(Some("5/2"), 30);
        for expression in [
            Atom::Zero,
            vk_parse!("(1-ep)^-1-1-ep").unwrap(),
            vk_parse!("ep^-2*((1-ep)^-1-1-ep)-(1-ep)^-1").unwrap(),
        ] {
            assert_eq!(
                catalog.retain(expression.clone(), &BTreeMap::new(), None),
                combine_exact(expression)
            );
        }
    }

    #[test]
    #[ignore = "complete saved q0 exact-coefficient fixture; run explicitly without reductions"]
    fn retained_saved_q0_preserves_every_complete_master_product_coefficient() {
        crate::Vakint::initialize_vakint_symbols();
        let provenance = include_bytes!("../../examples/data/five_loop_q0/provenance.json");
        checked_bytes(
            provenance,
            "941d5a2e888557e8cabd8e8ab118935e95631be0a16110bd28a41e49436db36a",
        )
        .unwrap();
        let provenance: serde_json::Value = serde_json::from_slice(provenance).unwrap();
        let bytes = include_bytes!("../../examples/data/five_loop_q0/q0-retained.atom");
        let identity = &provenance["artifacts"]["q0-retained.atom"];
        assert_eq!(bytes.len() as u64, identity["bytes"].as_u64().unwrap());
        checked_bytes(bytes, identity["sha256"].as_str().unwrap()).unwrap();
        let mut stream = bytes.as_slice();
        let original = Atom::import(&mut stream, None).unwrap();
        assert!(stream.is_empty());
        assert!(!MasterCatalog::has_floating_coefficients(&original));
        let symbols: Vec<_> = MasterCatalog::referenced_masters(&original)
            .into_iter()
            .collect();
        assert_eq!(symbols.len(), 55);
        let before: BTreeMap<_, _> = original
            .coefficient_list::<i32>(&symbols)
            .into_iter()
            .collect();
        assert_eq!(before.len(), 58);
        assert!(
            before
                .values()
                .all(|coefficient| { MasterCatalog::referenced_masters(coefficient).is_empty() })
        );
        // This fixture is already globally combined. Do not run the previous
        // whole-expression rationalization or infer pre-merge speed from it.
        let combined = MasterCatalog::combine_master_coefficients(original);
        assert_eq!(
            MasterCatalog::referenced_masters(&combined),
            symbols.iter().cloned().collect()
        );
        let after: BTreeMap<_, _> = combined
            .coefficient_list::<i32>(&symbols)
            .into_iter()
            .collect();
        assert!(
            after
                .values()
                .all(|coefficient| { MasterCatalog::referenced_masters(coefficient).is_empty() })
        );
        let products: BTreeSet<_> = before.keys().chain(after.keys()).cloned().collect();
        assert_eq!(products.len(), 58);
        for product in products {
            let before = before.get(&product).cloned().unwrap_or(Atom::Zero);
            let after = after.get(&product).cloned().unwrap_or(Atom::Zero);
            assert!(
                combine_exact(before - after).is_zero(),
                "complete scalar coefficient differs for {product}"
            );
        }
    }

    #[test]
    fn retained_catalogue_terminals_are_not_unknown_but_substitution_needs_coefficients() {
        let (mut catalog, terminal) = catalog(None, 30);
        catalog
            .masters
            .get_mut(&terminal)
            .unwrap()
            .coefficients
            .clear();
        assert!(
            catalog
                .unresolved(&terminal, false, &BTreeMap::new())
                .is_empty()
        );
        assert_eq!(
            catalog
                .finalize(
                    terminal.clone(),
                    &settings(16, 2),
                    1,
                    RustRedMasterMode::Retain
                )
                .unwrap(),
            terminal
        );
        for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
            let error = catalog
                .finalize(terminal.clone(), &settings(16, 2), 1, mode)
                .unwrap_err()
                .to_string();
            assert!(error.contains("test-master"), "{error}");
            assert!(
                error.contains("declared terminal has no computed coefficients"),
                "{error}"
            );
            assert!(
                error.contains("test-parent.dot") && error.contains("test-v1"),
                "{error}"
            );
            assert!(!error.contains("undeclared residual master"), "{error}");
        }
    }

    #[test]
    fn native_bare_declarations_supply_identity_without_fabricating_master_data() {
        let (catalog, _) = catalog(Some("5/2"), 30);
        let declared = catalog.symbol("bare-native-family", &[1, 1]);
        let unknown = catalog.symbol("bare-native-family", &[2, 1]);
        let announced = BTreeMap::from([(
            declared.clone(),
            BareMasterDeclaration {
                application_path: "native/application.json".into(),
            },
        )]);
        assert!(catalog.unresolved(&declared, false, &announced).is_empty());
        assert_eq!(catalog.retain(declared.clone(), &announced, None), declared);
        let error = catalog
            .require_coefficients(&declared, &announced)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("declared terminal has no computed coefficients"),
            "{error}"
        );
        assert!(
            error.contains("native/application.json")
                && error.contains("matching external-data pack"),
            "{error}"
        );
        assert!(!error.contains("undeclared residual master"), "{error}");
        let unknown_issues = catalog.unresolved(&unknown, false, &announced);
        assert_eq!(unknown_issues.len(), 1);
        assert!(unknown_issues[0].contains("undeclared residual master"));
        // The native view is read-only; no data record, precision or graph was
        // added to the numerical registry by retaining an announced terminal.
        assert!(!catalog.masters.contains_key(&declared));
        let zero = vk_parse!("ep^-2*((1-ep)^-1-1-ep)-(1-ep)^-1").unwrap();
        let combined = combine_exact(zero * declared);
        assert!(catalog.require_coefficients(&combined, &announced).is_ok());
    }

    #[test]
    fn retained_postprocessing_phases_preserve_exact_cancellation_and_detachment() {
        use rustred::solver::{
            CandidateRoutedApplicationEnvelope, CandidateRoutedApplicationObserver,
            CandidateRoutedObservationOptions, CandidateRoutedObserverControl,
            CandidateRoutedObserverStatus,
        };
        let (catalog, master) = catalog(Some("5/2"), 30);
        let unknown = catalog.symbol("unannounced-phase-control", &[2]);
        let zero = vk_parse!("ep^-2*((1-ep)^-1-1-ep)-(1-ep)^-1").unwrap();
        let raw = &master + zero * unknown;
        for detach in [false, true] {
            let mut phases = Vec::new();
            {
                let mut callback = |event: &CandidateRoutedApplicationEnvelope<'_>| {
                    phases.push(event.caller_phase.unwrap().to_owned());
                    Ok(if detach {
                        CandidateRoutedObserverControl::Detach
                    } else {
                        CandidateRoutedObserverControl::Continue
                    })
                };
                let mut observer = CandidateRoutedApplicationObserver::new(
                    CandidateRoutedObservationOptions {
                        interval: std::time::Duration::ZERO,
                        ..Default::default()
                    },
                    &mut callback,
                );
                assert_eq!(
                    catalog.retain(raw.clone(), &BTreeMap::new(), Some(&mut observer)),
                    master,
                );
                assert_eq!(
                    observer.status(),
                    if detach {
                        CandidateRoutedObserverStatus::Detached
                    } else {
                        CandidateRoutedObserverStatus::Active
                    },
                );
            }
            let expected = if detach {
                vec!["retained_catalogue_exact_combination"]
            } else {
                vec![
                    "retained_catalogue_exact_combination",
                    "retained_catalogue_classification",
                ]
            };
            assert_eq!(phases, expected);
        }
    }

    #[test]
    fn residual_diagnostics_are_sorted_deduplicated_and_follow_exact_combination() {
        let (mut catalog, _) = catalog(Some("5/2"), 30);
        catalog.family_graphs.insert(
            MasterCatalog::family_symbol("test-family"),
            "test-parent.dot".into(),
        );
        let first = catalog.symbol("test-family", &[2]);
        let second = catalog.symbol("test-family", &[3]);
        let input = first.clone() + first.clone() + second.clone();
        let issues = catalog.unresolved(&input, false, &BTreeMap::new());
        assert_eq!(issues.len(), 2);
        assert!(issues.windows(2).all(|pair| pair[0] < pair[1]));
        assert!(issues.iter().all(|issue| issue.contains("test-parent.dot")));
        for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
            let error = catalog
                .finalize(input.clone(), &settings(16, 2), 1, mode)
                .unwrap_err()
                .to_string();
            assert_eq!(error.matches("undeclared residual master").count(), 2);
            assert!(
                error.contains("test-v1") && error.contains("\n  - "),
                "{error}"
            );
            // A nontrivial rational zero must be combined before classifying
            // terminals, in both substitution modes and in retained output.
            let zero = vk_parse!("ep^-2*((1-ep)^-1-1-ep)-(1-ep)^-1").unwrap();
            let cancelled = zero * input.clone();
            assert!(!cancelled.is_zero());
            assert!(
                catalog
                    .finalize(cancelled, &settings(16, 2), 1, mode)
                    .unwrap()
                    .is_zero()
            );
        }
        assert_eq!(
            catalog
                .finalize(
                    input.clone(),
                    &settings(16, 2),
                    1,
                    RustRedMasterMode::Retain
                )
                .unwrap(),
            combine_exact(input)
        );
    }

    fn evaluate(
        catalog: &MasterCatalog,
        expression: Atom,
        settings: &VakintSettings,
    ) -> crate::precision::NumericalEvaluationResult {
        crate::precision::evaluate_series(
            expression.as_view(),
            vk_symbol!("ep"),
            &HashMap::new(),
            settings,
            Some(catalog),
        )
        .unwrap()
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
        let prepared = catalog.prepare_evaluation(pole).unwrap();
        assert_eq!(prepared.expression, vk_parse!("ep^-1").unwrap());
        assert!(prepared.sources.is_empty());
        let expression = catalog
            .finalize(master, &settings(31, 2), 1, RustRedMasterMode::Numerical)
            .unwrap();
        let result = evaluate(&catalog, expression, &settings(31, 2));
        assert!(result.warnings().iter().any(|warning| {
            warning.requested_decimal_precision == 31
                && warning.sources.iter().any(|source| {
                    source.master_id.as_deref() == Some("test-master")
                        && source.epsilon_power == Some(0)
                        && source.available_decimal_precision == Some(30)
                })
        }));
    }

    #[test]
    fn symbolic_numerical_coefficients_resolve_with_active_data() {
        let (catalog, master) = catalog(None, 30);
        let symbolic = catalog
            .finalize(master, &settings(200, 2), 1, RustRedMasterMode::Symbolic)
            .unwrap();
        assert!(symbolic.contains_symbol(vk_symbol!("RustRedCoefficient")));
        let numerical = catalog
            .resolve_numerical(symbolic, &settings(31, 2))
            .unwrap();
        assert!(numerical.contains_symbol(vk_symbol!("RustRedCoefficient")));
        let prepared = catalog.prepare_evaluation(numerical.clone()).unwrap();
        assert_eq!(prepared.sources.len(), 1);
        assert_eq!(
            prepared.sources.values().next().unwrap().significant_digits,
            30
        );
        assert!(
            !evaluate(&catalog, numerical, &settings(31, 2))
                .warnings()
                .is_empty()
        );
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
    fn exact_lower_orders_cancel_a_missing_coefficient_in_numerical_mode() {
        let (catalog, master) = catalog(Some("5/2"), 80);
        let input = master.clone().pow(2) - Atom::num(2) * master / vk_parse!("ep").unwrap();
        let result = catalog
            .finalize(input, &settings(80, 3), 2, RustRedMasterMode::Numerical)
            .unwrap();
        let prepared = catalog.prepare_evaluation(result).unwrap();
        assert!(prepared.sources.is_empty());
        assert!(combine_exact(prepared.expression - vk_parse!("-ep^-2+25/4").unwrap()).is_zero());
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
    fn symbolic_formulas_bind_finite_periods_to_their_data_snapshot() {
        let (mut catalog, master) = catalog(Some("2*K5Moment(1)"), 30);
        let period = vk_parse!("K5Moment(1)").unwrap();
        catalog.constants.insert(
            period.clone(),
            crate::inspection::FiniteSource {
                value: "1.25".into(),
                significant_digits: 30,
                absolute_error: Some("1e-31".into()),
                error_is_bounded: true,
                identity: period.to_canonical_string(),
                master_id: None,
                epsilon_power: None,
                data_version: catalog.version.clone(),
                source_sha256: "certified-test-moment".into(),
            },
        );
        let symbolic = catalog
            .finalize(master, &settings(20, 2), 1, RustRedMasterMode::Symbolic)
            .unwrap();
        assert!(symbolic.contains_symbol(vk_symbol!("RustRedConstant")));
        let prepared = catalog.prepare_evaluation(symbolic.clone()).unwrap();
        assert_eq!(prepared.sources.len(), 1);
        assert_eq!(
            evaluate(&catalog, symbolic.clone(), &settings(20, 2)).coefficients()[&0]
                .value()
                .re,
            Float::parse("2.5", Some(100)).unwrap()
        );
        catalog.coefficient_tag = Atom::var(vk_symbol!("replacement_constant_data"));
        let error = catalog
            .prepare_evaluation(symbolic)
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
        assert!(numerical.contains_symbol(vk_symbol!("RustRedCoefficient")));
        let numerical = evaluate(&catalog, numerical, &settings(80, 2))
            .to_atom(vk_symbol!("ep"))
            .unwrap();
        assert!(!numerical.contains_symbol(vk_symbol!("RustRedCoefficient")));
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
        let expression = catalog
            .finalize(
                master.clone(),
                &settings(20, 2),
                1,
                RustRedMasterMode::Numerical,
            )
            .unwrap();
        let prepared = catalog.prepare_evaluation(expression.clone()).unwrap();
        assert_eq!(
            prepared.available_precision(),
            crate::precision::AvailablePrecision::FixedPrecision(5)
        );
        let result = evaluate(&catalog, expression, &settings(20, 2));
        assert!(result.warnings().iter().any(|warning| {
            warning
                .sources
                .iter()
                .any(|source| source.available_decimal_precision == Some(5))
        }));
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
            let numerical = evaluate(&catalog, numerical, &settings)
                .to_atom(vk_symbol!("ep"))
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
            let prepared = catalog.prepare_evaluation(coefficient).unwrap();
            let (actual, _) = crate::precision::evaluate_coefficient(
                prepared.expression.as_view(),
                &prepared.sources,
                &parameters,
                &settings,
            )
            .unwrap();
            assert!((actual.value().re.clone() - expected).norm() < tolerance);
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
    fn bare_alias_missing_higher_data_does_not_hide_corrupt_lower_formula() {
        let (mut catalog, alias) = catalog(Some("2"), 80);
        let (representative_catalog, _) = self::catalog(Some("2"), 80);
        let representative = catalog.symbol("test-family", &[2]);
        let mut rich = representative_catalog.masters.into_values().next().unwrap();
        rich.identity = "valid-high-order-representative".into();
        rich.coefficients.insert(
            1,
            Coefficient {
                exact: Some("7".into()),
                numerical: Some("7".into()),
                significant_digits: Some(80),
                absolute_error: None,
                error_is_bounded: false,
            },
        );
        catalog.masters.insert(representative.clone(), rich);
        catalog
            .masters
            .get_mut(&alias)
            .unwrap()
            .coefficients
            .get_mut(&0)
            .unwrap()
            .exact = Some("(".into());
        catalog.load_basis_aliases(vec![
            (
                1,
                representative.clone(),
                vk_parse!("corruption_test_basis").unwrap(),
            ),
            (
                1,
                alias.clone(),
                vk_parse!("corruption_test_basis").unwrap(),
            ),
        ]);
        assert!(catalog.basis_aliases.contains_key(&alias));
        assert!(
            catalog
                .finalize(
                    representative,
                    &settings(16, 3),
                    1,
                    RustRedMasterMode::Symbolic
                )
                .is_ok()
        );
        let error = catalog
            .finalize(alias, &settings(16, 3), 1, RustRedMasterMode::Symbolic)
            .unwrap_err();
        assert!(matches!(
            error,
            VakintError::RustRedEvaluation(RustRedEvaluationError::ExternalData { .. })
        ));
    }

    #[test]
    fn large_laurent_request_reports_missing_data_without_a_policy_work_limit() {
        let (catalog, master) = catalog(Some("2"), 80);
        let request = settings(16, 4099);
        assert_eq!(
            catalog
                .finalize_independent(master.clone(), &request, 1, RustRedMasterMode::Retain,)
                .unwrap(),
            master
        );
        assert!(
            catalog
                .finalize_independent(&master - &master, &request, 1, RustRedMasterMode::Symbolic,)
                .unwrap()
                .is_zero()
        );
        let error = catalog
            .finalize_independent(master, &request, 1, RustRedMasterMode::Symbolic)
            .unwrap_err();
        assert!(matches!(&error, VakintError::RustRedEvaluation(
            RustRedEvaluationError::MasterCapability { detail }
        ) if detail.contains("test-master")
            && detail.contains("requires unavailable epsilon coefficient 4097")
            && detail.contains("available powers [-1, 0]")
            && detail.contains("test-v1")));
        assert!(!error.to_string().contains("work limit"));
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
                metadata: serde_json::Value::Null,
                data: None,
                coefficients: BTreeMap::from([
                    (
                        -1,
                        Coefficient {
                            exact: Some("3".into()),
                            numerical: Some("3".into()),
                            significant_digits: Some(80),
                            absolute_error: None,
                            error_is_bounded: false,
                        },
                    ),
                    (
                        0,
                        Coefficient {
                            exact: Some("2".into()),
                            numerical: Some("2".into()),
                            significant_digits: Some(80),
                            absolute_error: None,
                            error_is_bounded: false,
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
            (alias.clone() - vk_parse!("3-4*ep").unwrap() * &master) / vk_parse!("ep").unwrap();
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

        // Distinct closed dummy names must not manufacture a pole requiring
        // unavailable orders before the declared basis relation can cancel.
        let left = vk_parse!("p(1,11)*p(2,11)").unwrap();
        let right = vk_parse!("p(1,12)*p(2,12)").unwrap();
        let raw = (&alias * &left - vk_parse!("3-4*ep").unwrap() * &master * &right)
            / vk_parse!("ep").unwrap();
        assert!(
            catalog
                .finalize(
                    raw.clone(),
                    &settings(80, 2),
                    1,
                    RustRedMasterMode::Numerical
                )
                .is_err()
        );
        let combined =
            crate::VakintTerm::scalar_product_form(&settings(80, 2), raw.as_view()).unwrap();
        assert!(
            catalog
                .finalize(combined, &settings(80, 2), 1, RustRedMasterMode::Numerical)
                .unwrap()
                .is_zero()
        );
        let uncancelled = &alias * &left / vk_parse!("ep").unwrap();
        let uncancelled =
            crate::VakintTerm::scalar_product_form(&settings(80, 2), uncancelled.as_view())
                .unwrap();
        assert!(
            catalog
                .finalize(
                    uncancelled,
                    &settings(80, 2),
                    1,
                    RustRedMasterMode::Numerical
                )
                .is_err()
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
                    absolute_error: None,
                    error_is_bounded: false,
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
    fn family_graph_remains_the_fallback_without_declared_parent_metadata() {
        let vendored = super::super::data::RustRedData::vendored().unwrap();
        let mut package = (*vendored.loops[&1]).clone();
        package.manifest.masters.clear();
        let package = Arc::new(package);
        let family = &package.manifest.families[0];
        let fingerprint = family["fingerprint"].as_str().unwrap();
        let expected =
            MasterCatalog::graph_location(&package, family["graph"].as_str().unwrap()).unwrap();
        let catalog = MasterCatalog::load(
            &BTreeMap::from([(1, Arc::clone(&package))]),
            "local-mixed-release",
        )
        .unwrap();
        assert!(catalog.masters.is_empty());
        assert_eq!(
            catalog.family_graphs[&MasterCatalog::family_symbol(fingerprint)],
            expected
        );
        let unknown = catalog.symbol(fingerprint, &[2]);
        let issues = catalog.unresolved(&unknown, false, &BTreeMap::new());
        assert_eq!(issues.len(), 1);
        assert!(issues[0].contains(&expected));
        assert!(issues[0].contains("no unique parent inferred"));
        assert!(expected.contains(&format!("/blob/{}/1L/", package.manifest.version)));
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
    #[test]
    fn finite_period_aliases_share_identity_and_best_authenticated_precision() {
        let (mut catalog, first) = catalog(None, 5);
        let (duplicate, _) = self::catalog(None, 50);
        let second = catalog.symbol("test-family", &[2]);
        let mut second_data = duplicate.masters.into_values().next().unwrap();
        second_data.identity = "higher-precision-equivalent".into();
        second_data.loop_count = 4;
        second_data.metadata = serde_json::json!({"exact_basis_expression":"vakint::PR4"});
        catalog.masters.insert(second.clone(), second_data);
        let first_data = catalog.masters.get_mut(&first).unwrap();
        first_data.loop_count = 4;
        first_data.metadata = serde_json::json!({"exact_basis_expression":"vakint::PR4"});
        let first_coefficient = catalog.coefficient_symbol(&first, 0);
        let second_coefficient = catalog.coefficient_symbol(&second, 0);
        let prepared = catalog
            .prepare_evaluation(first_coefficient.clone())
            .unwrap();
        assert_eq!(
            prepared.available_precision(),
            crate::precision::AvailablePrecision::FixedPrecision(50)
        );
        assert_eq!(
            prepared
                .sources
                .values()
                .next()
                .unwrap()
                .master_id
                .as_deref(),
            Some("higher-precision-equivalent")
        );
        let zero = catalog
            .prepare_evaluation(
                first_coefficient + second_coefficient
                    - Atom::num(2) * vk_parse!("PR4fin").unwrap(),
            )
            .unwrap();
        assert!(zero.expression.is_zero());
        assert!(zero.sources.is_empty());
        assert!(zero.origins.is_empty());
    }

    #[test]
    fn finite_formula_leaves_preserve_their_cap_without_metadata_recursion() {
        let (catalog, master) = catalog(Some("2.5`100+symbolica::zeta(3)"), 5);
        let prepared = catalog
            .prepare_evaluation(catalog.coefficient_symbol(&master, 0))
            .unwrap();
        assert_eq!(
            prepared.available_precision(),
            crate::precision::AvailablePrecision::FixedPrecision(5)
        );
        assert_eq!(prepared.sources.len(), 1);
        assert!(!matches!(prepared.expression.as_view(), AtomView::Num(_)));
        assert!(
            prepared
                .origins
                .iter()
                .any(|origin| origin.master_id.as_deref() == Some("test-master")
                    && origin.epsilon_power == Some(0)
                    && origin.available_decimal_precision == Some(5))
        );
    }

    #[test]
    fn lookup_names_are_typed_and_family_ambiguity_is_explicit() {
        use crate::inspection::{ConstantSymbol, InspectionError, MasterSymbol};
        let (mut catalog, first) = catalog(Some("1"), 20);
        catalog.masters.get_mut(&first).unwrap().metadata = serde_json::json!({"name":"P1"});
        let (second_catalog, _) = self::catalog(Some("2"), 20);
        let mut other = second_catalog.masters.into_values().next().unwrap();
        other.metadata = serde_json::json!({"name":"P1"});
        let second = catalog.symbol("another-family", &[1]);
        catalog.masters.insert(second, other);
        let catalog = Arc::new(catalog);
        assert!(matches!(
            catalog.get_constant("foreign::z3".parse::<ConstantSymbol>().unwrap(), 16),
            Err(VakintError::Inspection(
                InspectionError::UnknownConstant { .. }
            ))
        ));
        assert!(matches!(
            catalog.get_master("missing".parse::<MasterSymbol>().unwrap(), 16),
            Err(VakintError::Inspection(
                InspectionError::UnknownMaster { .. }
            ))
        ));
        assert!(
            matches!(catalog.get_master("P1".parse::<MasterSymbol>().unwrap(), 16), Err(VakintError::Inspection(InspectionError::AmbiguousMaster { candidates, .. })) if candidates.len() == 2)
        );
        assert!(
            catalog
                .master_symbols()
                .iter()
                .any(|symbol| symbol.as_atom() == &first)
        );
    }
    #[test]
    fn bounded_period_lookup_caps_storage_and_keeps_its_snapshot() {
        let (mut catalog, symbol) = catalog(None, 30);
        let master = catalog.masters.get_mut(&symbol).unwrap();
        master.loop_count = 4;
        master.metadata = serde_json::json!({"exact_basis_expression":"vakint::PR4"});
        let coefficient = master.coefficients.get_mut(&0).unwrap();
        coefficient.absolute_error = Some("1e-38".into());
        coefficient.error_is_bounded = true;
        let catalog = Arc::new(catalog);
        let constant = catalog.get_constant("PR4fin".parse().unwrap(), 40).unwrap();
        assert_eq!(
            constant.available_precision(),
            crate::precision::AvailablePrecision::FixedPrecision(30)
        );
        assert!(constant.numerical_value().prec() <= crate::precision::decimal_bits(30));
        assert_eq!(
            constant.precision().real().quality(),
            crate::precision::ErrorQuality::Bounded
        );
        assert!(constant.warnings().iter().any(|warning| {
            warning.requested_decimal_precision == 40
                && warning
                    .sources
                    .iter()
                    .any(|source| source.available_decimal_precision == Some(30))
        }));
        assert!(Arc::strong_count(&catalog) > 1);
        drop(catalog);
        assert_eq!(constant.numerical_value(), &Float::with_val(128, 2.5));
    }
    #[test]
    fn an_uncertain_real_zero_remains_inspectable_but_cannot_be_an_atom() {
        let (mut catalog, symbol) = catalog(None, 30);
        let master = catalog.masters.get_mut(&symbol).unwrap();
        master.loop_count = 4;
        master.metadata = serde_json::json!({"exact_basis_expression":"vakint::PR4"});
        let coefficient = master.coefficients.get_mut(&0).unwrap();
        coefficient.numerical = Some("0".into());
        coefficient.absolute_error = Some("1e-20".into());
        coefficient.error_is_bounded = true;
        let constant = Arc::new(catalog)
            .get_constant("PR4fin".parse().unwrap(), 16)
            .unwrap();
        assert!(constant.numerical_value().as_raw().is_zero());
        assert!(!constant.precision().real().exact_zero());
        assert!(constant.precision().real().absolute_error() > &Float::with_val(64, 0));
        assert!(!constant.warnings().is_empty());
        assert!(constant.to_atom().is_err());
    }
    #[test]
    fn instance_numerical_evaluation_reports_terminal_status_before_tag_evaluation() {
        let (mut catalog, master) = catalog(Some("5/2"), 80);
        let settings = settings(40, 2);
        let evaluate = |catalog: &MasterCatalog, expression: &Atom| {
            crate::precision::evaluate_series(
                expression.as_view(),
                vk_symbol!("ep"),
                &HashMap::new(),
                &settings,
                Some(catalog),
            )
        };
        let error = evaluate(&catalog, &master).unwrap_err().to_string();
        assert!(error.contains("Retained masters require") && error.contains("RustRedMaster"));
        let unknown = catalog.symbol("undeclared-family", &[2]);
        let error = evaluate(&catalog, &unknown).unwrap_err().to_string();
        assert!(
            error.contains("undeclared residual master") && error.contains("parent graph"),
            "{error}"
        );
        catalog
            .masters
            .get_mut(&master)
            .unwrap()
            .coefficients
            .clear();
        let error = evaluate(&catalog, &master).unwrap_err().to_string();
        assert!(
            error.contains("declared terminal has no computed coefficients")
                && error.contains("test-parent.dot"),
            "{error}"
        );
        assert!(
            evaluate(
                &catalog,
                &(unknown.clone() - unknown + master.clone() - master)
            )
            .is_ok()
        );
    }
}
