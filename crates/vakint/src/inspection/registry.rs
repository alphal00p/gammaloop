use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::sync::Arc;

use sha2::{Digest, Sha256};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    coefficient::CoefficientView,
    domains::float::Float,
};

use super::{
    ConstantSymbol, DataProvenance, FiniteSource, InspectionError, MasterCoefficient, MasterSymbol,
    MasterSymbolicExpression, PreparedExpression, VakintConstant, VakintMaster, constants, graph,
};
use crate::{
    VakintError, VakintSettings, precision,
    rustred_evaluation::{
        RustRedEvaluationError,
        masters::{Master, MasterCatalog, combine_exact},
    },
    utils::vakint_macros::{vk_parse, vk_symbol},
};

fn invalid(detail: impl std::fmt::Display) -> VakintError {
    RustRedEvaluationError::ExternalData {
        detail: detail.to_string(),
    }
    .into()
}

// These names identify recorded coefficients, not independent decimal tables.
const PERIOD_BINDINGS: [(&str, &str, i64); 5] = [
    ("PR4fin", "PR4", 0),
    ("PR4dfin", "PR4d", 0),
    ("PR11dep0", "PR11d", 0),
    ("PR4ep1", "PR4", 1),
    ("PR4dep1", "PR4d", 1),
];

#[derive(Clone)]
struct RecordedOrigin {
    source: precision::PrecisionSource,
    dependencies: BTreeSet<Atom>,
}

// Only semantic leaves participate in attribution. In particular, the tags
// inside a finite source are not mathematical dependencies of the result.
fn dependencies(expression: &Atom) -> BTreeSet<Atom> {
    let mut result = BTreeSet::new();
    expression.visitor(&mut |part| {
        if matches!(part, AtomView::Var(_) | AtomView::Fun(_)) {
            result.insert(part.to_owned());
            false
        } else {
            true
        }
    });
    result
}

impl MasterCatalog {
    fn provenance(&self, master: &Master, power: Option<i64>) -> DataProvenance {
        DataProvenance {
            version: master
                .data
                .as_ref()
                .map_or(&self.version, |data| &data.manifest.version)
                .clone(),
            source_sha256: format!(
                "{:x}",
                Sha256::digest(master.metadata.to_string().as_bytes())
            ),
            master_id: Some(master.identity.clone()),
            epsilon_power: power,
        }
    }

    /// Enumerate canonical identities, including declared but uncomputed terminals.
    pub(crate) fn master_symbols(&self) -> Vec<MasterSymbol> {
        self.masters.keys().cloned().map(MasterSymbol).collect()
    }

    fn lookup_master(&self, requested: &MasterSymbol) -> Result<(Atom, &Master), VakintError> {
        if let Some(master) = self.masters.get(requested.as_atom()) {
            return Ok((requested.as_atom().clone(), master));
        }
        let matches: Vec<_> = self
            .masters
            .iter()
            .filter(|(_, master)| {
                [
                    master.identity.as_str(),
                    master.metadata["name"].as_str().unwrap_or(""),
                ]
                .iter()
                .any(|name| {
                    !name.is_empty()
                        && vk_parse!(name).is_ok_and(|atom| atom == requested.as_atom())
                }) || master.metadata["exact_basis_expression"]
                    .as_str()
                    .is_some_and(|expression| {
                        vk_parse!(expression).is_ok_and(|atom| atom == requested.as_atom())
                    })
            })
            .collect();
        match matches.as_slice() {
            [] => Err(InspectionError::UnknownMaster {
                symbol: requested.to_string(),
            }
            .into()),
            [(symbol, master)] => Ok(((*symbol).clone(), master)),
            _ => Err(InspectionError::AmbiguousMaster {
                symbol: requested.to_string(),
                candidates: matches
                    .iter()
                    .map(|(symbol, _)| symbol.to_canonical_string())
                    .collect(),
            }
            .into()),
        }
    }

    /// Finite period aliases are exact references to a recorded master coefficient.
    pub(crate) fn is_recorded_period(symbol: AtomView<'_>) -> bool {
        PERIOD_BINDINGS
            .iter()
            .any(|(name, _, _)| symbol == Atom::var(vk_symbol!(name)).as_view())
    }

    pub(super) fn period_binding(&self, symbol: AtomView<'_>) -> Result<Option<Atom>, VakintError> {
        if let Some((symbol, _)) = self.constant_source(symbol)? {
            return Ok(Some(self.constant_symbol(symbol)));
        }
        let binding = PERIOD_BINDINGS
            .into_iter()
            .find_map(|(name, basis, power)| {
                (symbol == Atom::var(vk_symbol!(name)).as_view()).then_some((basis, power))
            });
        let Some((basis, power)) = binding else {
            return Ok(None);
        };
        let basis = Atom::var(vk_symbol!(basis));
        let candidate = self.masters.iter().find(|(_, master)| {
            master.loop_count == 4
                && master.metadata["exact_basis_expression"]
                    .as_str()
                    .is_some_and(|value| vk_parse!(value).is_ok_and(|atom| atom == basis))
                && master.coefficients.get(&power).is_some_and(|coefficient| {
                    coefficient.numerical.is_some() || coefficient.exact.is_some()
                })
        });
        Ok(candidate.map(|(master, _)| self.coefficient_symbol(master, power)))
    }

    /// Bind finite leaves without expanding definitions or selecting numerical payloads.
    pub(crate) fn bind_period_references(&self, expression: Atom) -> Result<Atom, VakintError> {
        let mut replacements = BTreeMap::new();
        let mut error = None;
        expression.visitor(&mut |part| {
            if error.is_some()
                || precision::snapshot_reference(part)
                || precision::legacy::is_source(part)
            {
                return false;
            }
            match self.period_binding(part) {
                Ok(Some(binding)) => {
                    replacements.insert(part.to_owned(), binding);
                    false
                }
                Ok(None) => true,
                Err(value) => {
                    error = Some(value);
                    false
                }
            }
        });
        if let Some(error) = error {
            return Err(error);
        }
        let replacements = replacements
            .into_iter()
            .map(|(from, to)| symbolica::id::Replacement::new(from, to))
            .collect::<Vec<_>>();
        Ok(expression.replace_multiple(&replacements))
    }

    fn constant_symbol(&self, symbol: &Atom) -> Atom {
        // The same name in another snapshot may have different validated data.
        vk_symbol!("RustRedConstant").call((self.coefficient_snapshot_tag(), symbol))
    }

    fn constant_source(
        &self,
        symbol: AtomView<'_>,
    ) -> Result<Option<(&Atom, &FiniteSource)>, VakintError> {
        if let AtomView::Fun(function) = symbol
            && function.get_symbol() == vk_symbol!("RustRedConstant")
        {
            let arguments = function.iter().collect::<Vec<_>>();
            let [tag, identity] = arguments.as_slice() else {
                return Err(invalid("malformed finite constant reference"));
            };
            if *tag != self.coefficient_snapshot_tag().as_view() {
                return Err(invalid(format!(
                    "finite constant reference belongs to another data snapshot: {symbol}"
                )));
            }
            return self
                .constants
                .get_key_value(&(*identity).to_owned())
                .map(Some)
                .ok_or_else(|| invalid(format!("unknown finite constant reference {symbol}")));
        }
        Ok(self.constants.get_key_value(&symbol.to_owned()))
    }

    // An exact basis identity fixes a coefficient independently of how many
    // separately rounded files happen to contain it. Choose its identity and
    // its best available payload independently, so precision upgrades cannot
    // break algebraic cancellation between an alias and a named period.
    fn coefficient_group<'a>(
        &'a self,
        requested: &Atom,
        master: &'a Master,
        power: i64,
    ) -> Vec<(Atom, &'a Master)> {
        let Some(basis) = master.metadata["exact_basis_expression"]
            .as_str()
            .and_then(|value| vk_parse!(value).ok())
            .filter(|value| Self::finite_precision(value).is_none())
        else {
            return vec![(requested.clone(), master)];
        };
        let candidates: Vec<_> = self
            .masters
            .iter()
            .filter(|(_, other)| {
                other.loop_count == master.loop_count
                    && other.coefficients.contains_key(&power)
                    && other.metadata["exact_basis_expression"]
                        .as_str()
                        .is_some_and(|value| vk_parse!(value).is_ok_and(|value| value == basis))
            })
            .map(|(symbol, master)| (self.coefficient_symbol(symbol, power), master))
            .collect();
        if candidates.is_empty() {
            vec![(requested.clone(), master)]
        } else {
            candidates
        }
    }

    // Select an authenticated representation for the entire alias group before
    // exposing a leaf. This is capability analysis, with no numerical evaluation
    // or retry. Exact formulas win ties; a finite formula must not replace a
    // more accurate independent numerical payload.
    fn prepare_coefficient(
        &self,
        requested: &Atom,
        master: &Master,
        power: i64,
        sources: &mut BTreeMap<Atom, FiniteSource>,
        visiting: &mut BTreeSet<Atom>,
        origins: &mut BTreeMap<Atom, RecordedOrigin>,
    ) -> Result<Atom, VakintError> {
        let candidates = self.coefficient_group(requested, master, power);
        let canonical = &candidates[0].0;
        if !visiting.insert(canonical.clone()) {
            return Err(invalid(format!(
                "cyclic master coefficient definition {canonical}"
            )));
        }
        let mut best = None;
        for (identity, candidate) in &candidates {
            let coefficient = &candidate.coefficients[&power];
            if let Some(formula) = &coefficient.exact {
                let mut trial_sources = sources.clone();
                let mut trial_origins = origins.clone();
                let formula = self.bind_finite_literals(
                    vk_parse!(formula).map_err(invalid)?,
                    identity,
                    coefficient.significant_digits,
                    &self.provenance(candidate, Some(power)),
                    &mut trial_sources,
                )?;
                let value =
                    self.prepare_inner(formula, &mut trial_sources, visiting, &mut trial_origins)?;
                let value = combine_exact(crate::hpl::normalize_constants(value.as_view())?);
                let capability = Self::dependency_precision(&value, &trial_sources);
                // None denotes arbitrary precision and ranks above every cap.
                let rank = (capability.is_none(), capability.unwrap_or(0), true);
                if best
                    .as_ref()
                    .is_none_or(|(old_rank, _, _, _)| rank > *old_rank)
                {
                    best = Some((rank, value, trial_sources, trial_origins));
                }
            }
            if let Some(value) = &coefficient.numerical {
                let significant_digits = coefficient
                    .significant_digits
                    .filter(|digits| *digits > 0)
                    .ok_or_else(|| invalid("finite coefficient has no validated precision"))?;
                Self::validate_decimal(value, &candidate.identity, power)?;
                let rank = (false, significant_digits, false);
                if best
                    .as_ref()
                    .is_none_or(|(old_rank, _, _, _)| rank > *old_rank)
                {
                    let mut trial_sources = sources.clone();
                    let provenance = self.provenance(candidate, Some(power));
                    let identity = sources
                        .get(canonical)
                        .filter(|source| {
                            source.value == *value
                                && source.significant_digits == significant_digits
                                && source.source_sha256 == provenance.source_sha256
                        })
                        .map_or_else(
                            || format!("{} epsilon^{power}", candidate.identity),
                            |source| source.identity.clone(),
                        );
                    trial_sources.insert(
                        canonical.clone(),
                        FiniteSource {
                            value: value.clone(),
                            significant_digits,
                            absolute_error: coefficient.absolute_error.clone(),
                            error_is_bounded: coefficient.error_is_bounded,
                            identity,
                            master_id: Some(candidate.identity.clone()),
                            epsilon_power: Some(power),
                            data_version: provenance.version,
                            source_sha256: provenance.source_sha256,
                        },
                    );
                    best = Some((rank, canonical.clone(), trial_sources, origins.clone()));
                }
            }
        }
        visiting.remove(canonical);
        let (_, value, selected_sources, selected_origins) =
            best.ok_or_else(|| invalid("coefficient has no usable representation"))?;
        *sources = selected_sources;
        *origins = selected_origins;
        Ok(value)
    }

    pub(crate) fn constant_symbols(&self) -> Result<Vec<ConstantSymbol>, VakintError> {
        let mut symbols: BTreeSet<_> = constants::builtin_symbols().into_iter().collect();
        symbols.extend(self.constants.keys().cloned());
        for (name, _, _) in PERIOD_BINDINGS {
            let atom = Atom::var(vk_symbol!(name));
            if self.period_binding(atom.as_view())?.is_some() {
                symbols.insert(atom);
            }
        }
        for master in self.masters.values() {
            for coefficient in master.coefficients.values() {
                if let Some(formula) = &coefficient.exact {
                    let expression = vk_parse!(formula).map_err(invalid)?;
                    let expression = crate::hpl::normalize_constants(expression.as_view())?;
                    let mut error = None;
                    expression.visitor(&mut |part| {
                        match constants::builtin_definition(part) {
                            Ok(Some(_)) => {
                                symbols.insert(part.to_owned());
                                return false;
                            }
                            Ok(None) => {}
                            Err(value) => error = Some(value),
                        }
                        true
                    });
                    if let Some(error) = error {
                        return Err(error);
                    }
                }
            }
        }
        Ok(symbols.into_iter().map(ConstantSymbol).collect())
    }

    /// Expand exact definitions while retaining the identity and accuracy of finite leaves.
    pub(crate) fn prepare_evaluation(
        &self,
        expression: Atom,
    ) -> Result<PreparedExpression, VakintError> {
        let mut sources = BTreeMap::new();
        let mut origins = BTreeMap::new();
        let expression = self.prepare_inner(
            combine_exact(expression),
            &mut sources,
            &mut BTreeSet::new(),
            &mut origins,
        )?;
        let expression = combine_exact(crate::hpl::normalize_constants(expression.as_view())?);
        // A source that cancels exactly is no longer a limitation of this value.
        let mut used = BTreeSet::new();
        expression.visitor(&mut |part| {
            if sources.contains_key(&part.to_owned()) {
                used.insert(part.to_owned());
                return false;
            }
            true
        });
        sources.retain(|key, _| used.contains(key));
        let surviving = dependencies(&expression);
        let origins = origins
            .into_values()
            .filter(|origin| !origin.dependencies.is_disjoint(&surviving))
            .map(|origin| origin.source)
            .collect();
        Ok(PreparedExpression {
            expression,
            sources,
            origins,
        })
    }

    fn prepare_inner(
        &self,
        expression: Atom,
        sources: &mut BTreeMap<Atom, FiniteSource>,
        visiting: &mut BTreeSet<Atom>,
        origins: &mut BTreeMap<Atom, RecordedOrigin>,
    ) -> Result<Atom, VakintError> {
        if visiting.len() > 128 {
            return Err(invalid("constant definition dependency depth exceeds 128"));
        }
        let mut replacements = BTreeMap::new();
        let mut error = None;
        expression.visitor(&mut |part| {
            if error.is_some() { return false }
            let symbol = part.to_owned();
            // The arguments of tracked leaves are identity metadata, not
            // expression children or recursive coefficient definitions.
            if sources.contains_key(&symbol)
                && !matches!(part, AtomView::Fun(fun) if fun.get_symbol() == vk_symbol!("RustRedCoefficient"))
            { return false }
            let replacement = (|| -> Result<Option<Atom>, VakintError> {
                if let Some((identity, source)) = self.constant_source(part)? {
                    let key = self.constant_symbol(identity);
                    origins.insert(key.clone(), RecordedOrigin {
                        source: precision::PrecisionSource::from(source),
                        dependencies: BTreeSet::from([key.clone()]),
                    });
                    sources.insert(key.clone(), source.clone());
                    return Ok(Some(key));
                }
                if let AtomView::Fun(function) = part
                    && function.get_symbol() == vk_symbol!("RustRedFiniteLiteral")
                {
                    let arguments = function.iter().collect::<Vec<_>>();
                    if arguments.len() != 3 { return Err(invalid("malformed finite coefficient literal")); }
                    let identity = arguments[0].to_owned();
                    let (master, coefficient, power) = self.coefficient(&identity)?;
                    let formula = coefficient.exact.as_ref().ok_or_else(|| invalid("finite literal has no owning formula"))?;
                    let mut recovered = BTreeMap::new();
                    let _ = self.bind_finite_literals(vk_parse!(formula).map_err(invalid)?, &identity,
                        coefficient.significant_digits, &self.provenance(master, Some(power)), &mut recovered)?;
                    let source = recovered.remove(&symbol).ok_or_else(|| invalid("finite literal identity disagrees with its owning formula"))?;
                    origins.insert(symbol.clone(), RecordedOrigin {
                        source: precision::PrecisionSource::from(&source), dependencies: BTreeSet::from([symbol.clone()]),
                    });
                    sources.insert(symbol.clone(), source);
                    return Ok(Some(symbol.clone()));
                }
                if precision::legacy::is_source(part) {
                    if let Some(original) = precision::legacy::original_source(part)?
                        && self.period_binding(original.as_view())?.is_some()
                    {
                        return self.prepare_inner(original, sources, visiting, origins).map(Some);
                    }
                    // Preserve historical payload and source identity together;
                    // none of the encoded arguments are expression children.
                    sources.extend(precision::legacy::collect_sources(part)?);
                    return Ok(Some(symbol.clone()));
                }
                if matches!(part, AtomView::Fun(fun) if fun.get_symbol() == vk_symbol!("RustRedCoefficient")) {
                    let (master, _, power) = self.coefficient(&symbol)?;
                    let value = self.prepare_coefficient(&symbol, master, power, sources, visiting, origins)?;
                    origins.insert(symbol.clone(), RecordedOrigin {
                        source: precision::PrecisionSource {
                        identity: format!("{} epsilon^{power}", master.identity),
                        master_id: Some(master.identity.clone()), epsilon_power: Some(power),
                        available_decimal_precision: Self::dependency_precision(&value, sources),
                        source_sha256: self.provenance(master, Some(power)).source_sha256,
                        },
                        dependencies: dependencies(&value),
                    });
                    return Ok(Some(value));
                }
                if let Some(binding) = self.period_binding(part)? {
                    let (master, _, power) = self.coefficient(&binding)?;
                    let provenance = self.provenance(master, Some(power));
                    let value = self.prepare_inner(binding, sources, visiting, origins)?;
                    if let Some(source) = sources.get_mut(&value) {
                        let name = symbol.to_canonical_string();
                        if !source.identity.contains(&name) { source.identity = format!("{name} ({})", source.identity); }
                    }
                    origins.insert(symbol.clone(), RecordedOrigin {
                        source: precision::PrecisionSource {
                            identity: symbol.to_canonical_string(), master_id: provenance.master_id,
                            epsilon_power: Some(power), available_decimal_precision: Self::dependency_precision(&value, sources),
                            source_sha256: provenance.source_sha256,
                        }, dependencies: dependencies(&value),
                    });
                    return Ok(Some(value));
                }
                if let Some(definition) = constants::builtin_definition(part)? {
                    let source = precision::PrecisionSource {
                        identity: symbol.to_canonical_string(), master_id: None, epsilon_power: None,
                        available_decimal_precision: None,
                        source_sha256: format!("{:x}", Sha256::digest(definition.to_canonical_string().as_bytes())),
                    };
                    if definition == symbol {
                        origins.entry(symbol.clone()).or_insert_with(|| RecordedOrigin {
                            source, dependencies: dependencies(&symbol),
                        });
                        return Ok(Some(symbol.clone()));
                    }
                    if !visiting.insert(symbol.clone()) { return Err(invalid(format!("cyclic constant definition {symbol}"))) }
                    let value = self.prepare_inner(definition, sources, visiting, origins)?;
                    visiting.remove(&symbol);
                    origins.entry(symbol.clone()).or_insert_with(|| RecordedOrigin {
                        source, dependencies: dependencies(&value),
                    });
                    return Ok(Some(value));
                }
                Ok(None)
            })();
            match replacement {
                Ok(Some(value)) => { if value != symbol { replacements.insert(symbol, value); } false }
                Ok(None) => true,
                Err(value) => { error = Some(value); false }
            }
        });
        if let Some(error) = error {
            return Err(error);
        }
        let replacements: Vec<_> = replacements
            .into_iter()
            .map(|(from, to)| symbolica::id::Replacement::new(from, to))
            .collect();
        Ok(expression.replace_multiple(&replacements))
    }

    fn dependency_precision(
        expression: &Atom,
        sources: &BTreeMap<Atom, FiniteSource>,
    ) -> Option<u32> {
        let mut digits: Option<u32> = None;
        expression.visitor(&mut |part| {
            if let Some(source) = sources.get(&part.to_owned()) {
                digits = Some(digits.map_or(source.significant_digits, |current| {
                    current.min(source.significant_digits)
                }));
                return false;
            }
            true
        });
        digits
    }

    pub(crate) fn validate_decimal(
        value: &str,
        master: &str,
        power: i64,
    ) -> Result<(), VakintError> {
        let number = Float::parse(value, Some(64)).map_err(|error| invalid(format!("master {master} coefficient epsilon^{power} has an invalid numerical value: {error}")))?;
        if !number.as_raw().is_finite() {
            return Err(invalid(format!(
                "master {master} coefficient epsilon^{power} has a non-finite numerical value"
            )));
        }
        Ok(())
    }

    fn bind_finite_literals(
        &self,
        expression: Atom,
        coefficient: &Atom,
        digits: Option<u32>,
        provenance: &DataProvenance,
        sources: &mut BTreeMap<Atom, FiniteSource>,
    ) -> Result<Atom, VakintError> {
        let mut replacements = BTreeMap::new();
        let mut error = None;
        expression.visitor(&mut |part| {
            let AtomView::Num(number) = part else { return true };
            let CoefficientView::Float(real, imaginary) = number.get_coeff_view() else { return true };
            let Some(digits) = digits.filter(|digits| *digits > 0) else {
                error = Some(invalid(format!("master {} coefficient epsilon^{} embeds a finite constant without validated significant digits; data {}", provenance.master_id.as_deref().unwrap_or("unknown"), provenance.epsilon_power.unwrap_or(0), provenance.version)));
                return false;
            };
            let mut value = Atom::Zero;
            for (index, component) in [real.to_float(), imaginary.to_float()].into_iter().enumerate() {
                if !component.as_raw().is_finite() { error = Some(invalid("non-finite literal in exact coefficient formula")); return false }
                if component.as_raw().is_zero() && index == 1 { continue }
                let digest = format!("rustred_literal_{:x}", Sha256::digest(part.to_canonical_string().as_bytes()));
                let key = vk_symbol!("RustRedFiniteLiteral").call((coefficient, Atom::var(vk_symbol!(&digest)), index as i64));
                let significant_digits = digits.min((f64::from(component.prec()) * std::f64::consts::LOG10_2).floor() as u32);
                sources.insert(key.clone(), FiniteSource {
                    value: component.as_raw().to_string_radix(10, None), significant_digits,
                    absolute_error: None, error_is_bounded: false,
                    identity: format!("{coefficient}: finite literal {part}"),
                    master_id: provenance.master_id.clone(), epsilon_power: provenance.epsilon_power,
                    data_version: provenance.version.clone(), source_sha256: provenance.source_sha256.clone(),
                });
                value += if index == 0 { key } else { key * vk_parse!("𝑖").unwrap() };
            }
            replacements.insert(part.to_owned(), value);
            false
        });
        if let Some(error) = error {
            return Err(error);
        }
        Ok(expression.replace_multiple(
            replacements
                .into_iter()
                .map(|(from, to)| symbolica::id::Replacement::new(from, to))
                .collect::<Vec<_>>(),
        ))
    }

    pub(crate) fn get_constant(
        self: &Arc<Self>,
        symbol: ConstantSymbol,
        prec: u32,
    ) -> Result<VakintConstant, VakintError> {
        if precision::legacy::is_source(symbol.as_atom().as_view()) {
            // These real leaves carry their source and precision themselves.
            // Validate the complete identity before admitting one as a constant;
            // a matching function name alone does not authenticate its payload.
            let _ = precision::legacy::collect_sources(symbol.as_atom().as_view())?;
        } else if constants::builtin_definition(symbol.as_atom().as_view())?.is_none()
            && self.period_binding(symbol.as_atom().as_view())?.is_none()
        {
            return Err(InspectionError::UnknownConstant {
                symbol: symbol.to_string(),
            }
            .into());
        }
        let prepared = self.prepare_evaluation(symbol.as_atom().clone())?;
        let available_precision = prepared.available_precision();
        let settings = VakintSettings {
            run_time_decimal_precision: prec,
            ..Default::default()
        };
        let (number, mut warnings) = precision::evaluate_coefficient(
            prepared.expression.as_view(),
            &prepared.sources,
            &HashMap::new(),
            &settings,
        )?;
        let cap = match available_precision {
            precision::AvailablePrecision::FixedPrecision(cap) => cap.min(prec),
            precision::AvailablePrecision::ArbitraryPrecision => prec,
        };
        let number = number.capped_decimal_precision(cap);
        for warning in &mut warnings {
            if warning.achieved_decimal_precision.is_some() {
                warning.achieved_decimal_precision = number.precision().real().significant_digits();
            }
        }
        number.real_midpoint()?;
        let mut provenance = prepared.sources.values().next().map_or_else(
            || DataProvenance {
                version: "exact-definition".into(),
                source_sha256: format!(
                    "{:x}",
                    Sha256::digest(prepared.expression.to_canonical_string().as_bytes())
                ),
                master_id: None,
                epsilon_power: None,
            },
            |source| DataProvenance {
                version: source.data_version.clone(),
                source_sha256: source.source_sha256.clone(),
                master_id: source.master_id.clone(),
                epsilon_power: source.epsilon_power,
            },
        );
        if let Some(binding) = self.period_binding(symbol.as_atom().as_view())?
            && matches!(binding.as_view(), AtomView::Fun(function) if function.get_symbol() == vk_symbol!("RustRedCoefficient"))
        {
            let (master, _, power) = self.coefficient(&binding)?;
            provenance = self.provenance(master, Some(power));
        }
        precision::attach_origins(&mut warnings, &prepared.origins, None, &provenance.version);
        for warning in &warnings {
            log::warn!("{warning}");
        }

        Ok(VakintConstant {
            symbol,
            requested_decimal_precision: prec,
            available_precision,
            numerical: number,
            warnings,
            provenance,
            _snapshot: Arc::clone(self),
        })
    }

    pub(crate) fn get_master(
        self: &Arc<Self>,
        requested: MasterSymbol,
        prec: u32,
    ) -> Result<VakintMaster, VakintError> {
        let (symbol, master) = self.lookup_master(&requested)?;
        let data = master
            .data
            .as_ref()
            .ok_or_else(|| InspectionError::UnavailableGraph {
                master: master.identity.clone(),
                detail: "no authenticated graph source".into(),
            })?;
        let graph = graph::master_graph(data, &master.metadata)?;
        let mut coefficients = BTreeMap::new();
        let settings = VakintSettings {
            run_time_decimal_precision: prec,
            ..Default::default()
        };
        for (&power, coefficient) in &master.coefficients {
            let identity = self.coefficient_symbol(&symbol, power);
            let prepared = self.prepare_evaluation(identity.clone())?;
            let (number, mut warnings) = precision::evaluate_coefficient(
                prepared.expression.as_view(),
                &prepared.sources,
                &HashMap::new(),
                &settings,
            )?;
            let cap = match prepared.available_precision() {
                precision::AvailablePrecision::FixedPrecision(cap) => cap.min(prec),
                precision::AvailablePrecision::ArbitraryPrecision => prec,
            };
            let number = number.capped_decimal_precision(cap);
            for warning in &mut warnings {
                if warning.achieved_decimal_precision.is_some() {
                    warning.achieved_decimal_precision =
                        number.precision().real().significant_digits();
                }
            }
            precision::attach_origins(
                &mut warnings,
                &prepared.origins,
                Some(power),
                &self.provenance(master, Some(power)).version,
            );
            for warning in &warnings {
                log::warn!("{warning}");
            }

            let symbolic_expression = match &coefficient.exact {
                Some(expression) => {
                    let expression = self.bind_finite_literals(
                        vk_parse!(expression).map_err(invalid)?,
                        &identity,
                        coefficient.significant_digits,
                        &self.provenance(master, Some(power)),
                        &mut BTreeMap::new(),
                    )?;
                    let expression = self.prepare_evaluation(expression)?.expression;
                    MasterSymbolicExpression::Expression(crate::hpl::to_symbolica_expression(
                        expression.as_view(),
                    )?)
                }
                None => {
                    match prepared.available_precision() {
                        precision::AvailablePrecision::FixedPrecision(digits) => {
                            MasterSymbolicExpression::NumericalOnly(digits)
                        }
                        precision::AvailablePrecision::ArbitraryPrecision => {
                            // An authenticated alias can supply an exact
                            // definition absent from this particular row.
                            MasterSymbolicExpression::Expression(
                                crate::hpl::to_symbolica_expression(prepared.expression.as_view())?,
                            )
                        }
                    }
                }
            };
            let provenance = if coefficient.exact.is_none() {
                prepared
                    .sources
                    .values()
                    .next()
                    .map(|source| DataProvenance {
                        version: source.data_version.clone(),
                        source_sha256: source.source_sha256.clone(),
                        master_id: source.master_id.clone(),
                        epsilon_power: source.epsilon_power,
                    })
                    .unwrap_or_else(|| self.provenance(master, Some(power)))
            } else {
                self.provenance(master, Some(power))
            };
            number.real_midpoint()?;
            coefficients.insert(
                power,
                MasterCoefficient {
                    symbol: identity,
                    requested_decimal_precision: prec,
                    epsilon_power: power,
                    available_precision: prepared.available_precision(),
                    numerical: number,
                    warnings,
                    provenance,
                    symbolic_expression,
                    _snapshot: Arc::clone(self),
                },
            );
        }
        Ok(VakintMaster {
            symbol: MasterSymbol(symbol),
            requested_decimal_precision: prec,
            master_id: master.identity.clone(),
            family_sha256: master.metadata["family_sha256"]
                .as_str()
                .ok_or_else(|| invalid("master has no family digest"))?
                .into(),
            loop_count: master.loop_count,
            normalization: data.manifest.normalization.clone(),
            provenance: self.provenance(master, None),
            graph,
            coefficients,
            _snapshot: Arc::clone(self),
        })
    }
}

#[cfg(test)]
#[path = "registry_tests.rs"]
mod tests;
