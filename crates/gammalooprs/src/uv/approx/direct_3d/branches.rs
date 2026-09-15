use std::{
    collections::BTreeMap,
    ops::Neg,
    sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    },
};

use color_eyre::Result;
use eyre::{ensure, eyre};
use spenso::shadowing::{ANTISYM, CYCLIC, SYM};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Indeterminate, Symbol, SymbolAttribute, SymbolBuilder},
    domains::rational::Rational,
    id::Replacement,
    poly::series::SeriesDepth,
    symbol,
};

use crate::{
    cff::{
        expression::{OrientationID, energy_map_replacements_gs},
        surface::LinearEnergyExpr,
    },
    graph::Graph,
    integrands::process::param_builder::FnMapEntry,
    numerator::symbolica_ext::NumeratorAtomExt,
    utils::{GS, external_energy_atom_from_index, ose_atom_from_index},
    uv::{
        Integrands,
        approx::{OrientationProjection, local_3d::OrientationIntegrands},
    },
};

static NUMERATOR_SCOPE: AtomicUsize = AtomicUsize::new(0);

/// The energy substitution belonging to one complete generalized residue key.
///
/// Stored production branches recover their map from `selector_host`. Reduced
/// sources retain their explicit map because it is part of their residue key,
/// not temporary projection metadata. Every factor in a branch is mapped with
/// this one authority; the Taylor operator only changes its rescaling series.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) enum DirectEnergyMap {
    Production,
    Source(Vec<LinearEnergyExpr>),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct DirectResidueKey {
    pub(crate) selector_host: OrientationID,
    pub(crate) energy_map: DirectEnergyMap,
}

impl DirectResidueKey {
    pub(crate) fn production(selector_host: OrientationID) -> Self {
        Self {
            selector_host,
            energy_map: DirectEnergyMap::Production,
        }
    }

    pub(crate) fn source(
        selector_host: OrientationID,
        edge_energy_map: Vec<LinearEnergyExpr>,
    ) -> Self {
        Self {
            selector_host,
            energy_map: DirectEnergyMap::Source(edge_energy_map),
        }
    }

    pub(crate) fn map_numerator(
        &self,
        orientation: OrientationProjection<'_>,
        graph: &Graph,
        numerator: &Atom,
    ) -> Result<Atom> {
        let source_map = match &self.energy_map {
            DirectEnergyMap::Production => None,
            DirectEnergyMap::Source(map) => Some(map.as_slice()),
        };
        orientation.map_numerator(graph, self.selector_host, source_map, numerator)
    }

    pub(crate) fn source_edge_energy_map(&self) -> Option<&[LinearEnergyExpr]> {
        match &self.energy_map {
            DirectEnergyMap::Production => None,
            DirectEnergyMap::Source(map) => Some(map),
        }
    }

    fn selector(&self, materialize_key_selector: bool) -> Atom {
        if materialize_key_selector {
            self.selector_host.atom()
        } else {
            Atom::one()
        }
    }
}

/// Factorized direct-3D bodies grouped by their complete generalized residue
/// key. Distinct source maps hosted by the same production selector remain
/// distinct branches, while exact duplicate keys are coalesced additively.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct DirectResidueBranches(Vec<(DirectResidueKey, Integrands)>);

impl DirectResidueBranches {
    /// Reserve a persistent scope tag, including names imported from archives.
    pub(crate) fn numerator_scope() -> (usize, Atom) {
        loop {
            let scope = NUMERATOR_SCOPE.fetch_add(1, Ordering::Relaxed);
            let name = format!("gammalooprs::uv::numerator_scope_{scope}");
            if !symbolica::state::State::symbol_iter().any(|(symbol, existing)| {
                existing == name || symbol.get_aliases().iter().any(|alias| alias == &name)
            }) {
                return (scope, Atom::var(symbol!(name.as_str())));
            }
        }
    }

    #[cfg(test)]
    pub(crate) fn production(selector_host: OrientationID, integrands: Integrands) -> Result<Self> {
        Self::from_keyed([(DirectResidueKey::production(selector_host), integrands)])
    }

    pub(crate) fn from_transient(source: &OrientationIntegrands) -> Result<Self> {
        Self::from_keyed(source.iter_orientations().map(
            |(selector_host, source_map, integrands)| {
                let key = source_map.map_or_else(
                    || DirectResidueKey::production(selector_host),
                    |map| DirectResidueKey::source(selector_host, map.to_vec()),
                );
                (key, integrands.clone())
            },
        ))
    }

    pub(super) fn from_keyed(
        keyed: impl IntoIterator<Item = (DirectResidueKey, Integrands)>,
    ) -> Result<Self> {
        let mut branches: Vec<(DirectResidueKey, Integrands)> = Vec::new();
        for (key, integrands) in keyed {
            if let Some((_, existing)) = branches
                .iter_mut()
                .find(|(existing_key, _)| *existing_key == key)
            {
                *existing = existing.clone().zip_add([integrands])?;
            } else {
                branches.push((key, integrands));
            }
        }
        let fallback_zero = branches
            .iter()
            .find(|(_, integrands)| integrands.iter().all(|(_, atom)| atom.is_zero()))
            .cloned();
        branches.retain(|(_, integrands)| !integrands.iter().all(|(_, atom)| atom.is_zero()));
        if branches.is_empty() {
            branches.extend(fallback_zero);
        }
        if branches.is_empty() {
            return Err(eyre!("direct local-3D residue branches cannot be empty"));
        }
        Ok(Self(branches))
    }

    pub(crate) fn iter_keys(&self) -> impl Iterator<Item = (&DirectResidueKey, &Integrands)> {
        self.0.iter().map(|(key, integrands)| (key, integrands))
    }

    /// Prepare one numerator with Taylor-inert coefficients of the complete
    /// affine energy maps. The returned tuple is `(numerator, parameters, rows)`;
    /// each row retains its complete residue key and ordered numeric arguments.
    /// OSEs, external energies and the continuous sampling scale M stay visible
    /// to Taylor, while the integer sampling node is a separate coefficient.
    #[allow(clippy::type_complexity)]
    pub(crate) fn prepare_numerator(
        &self,
        orientation: OrientationProjection<'_>,
        graph: &Graph,
        numerator: &Atom,
        scope: usize,
    ) -> Result<(Atom, Vec<Atom>, Vec<(DirectResidueKey, Vec<Atom>)>)> {
        let orientations = orientation.exact_orientations()?;
        let maps = self
            .iter_keys()
            .map(|(key, _)| {
                let production = orientations.get(key.selector_host).ok_or_else(|| {
                    eyre!(
                        "missing production energy map for orientation {}",
                        key.selector_host.0
                    )
                })?;
                Ok(key
                    .source_edge_energy_map()
                    .unwrap_or(&production.edge_energy_map))
            })
            .collect::<Result<Vec<_>>>()?;
        let edge_count = maps[0].len();
        if maps.iter().any(|map| map.len() != edge_count) {
            return Err(eyre!(
                "a direct numerator family requires equal edge-map support; absent entries are not zero samples"
            ));
        }

        let coefficient = symbol!("gammalooprs::uv::numerator_coefficient"; Scalar);
        let mut captured = false;
        for atom in std::iter::once(numerator).chain(
            self.iter_keys()
                .flat_map(|(_, integrands)| integrands.iter().map(|(_, atom)| atom)),
        ) {
            let _ = atom.replace_map(|view, _, _| {
                if let AtomView::Fun(function) = view
                    && function.get_symbol() == coefficient
                    && function.get_nargs() == 2
                    && usize::try_from(function.get(0)).ok() == Some(scope)
                {
                    captured = true;
                }
            });
        }
        if captured {
            return Err(eyre!(
                "direct numerator coefficient scope {scope} is already in use"
            ));
        }

        // The support is the union of full affine rows, not their orientation
        // signs. Distinct occurrences and distinct maps on one host keep their
        // own coefficients, including explicit zero samples.
        let mut support = BTreeMap::new();
        for (row, map) in maps.iter().enumerate() {
            for (slot, energy) in map.iter().enumerate() {
                let energy = energy.clone().canonical();
                let terms = energy
                    .internal_terms
                    .iter()
                    .map(|(edge, value)| (0, usize::from(*edge), value, ose_atom_from_index(*edge)))
                    .chain(energy.external_terms.iter().map(|(edge, value)| {
                        (
                            1,
                            usize::from(*edge),
                            value,
                            external_energy_atom_from_index(*edge),
                        )
                    }))
                    .chain([
                        (
                            2,
                            0,
                            &energy.uniform_scale_coeff,
                            Atom::var(GS.numerator_sampling_scale),
                        ),
                        (3, 0, &energy.constant, Atom::one()),
                    ]);
                for (kind, edge, value, basis) in terms {
                    let value = Atom::num(value.clone());
                    if !value.is_zero() {
                        support
                            .entry((slot, kind, edge))
                            .or_insert_with(|| (basis, vec![Atom::Zero; maps.len()]))
                            .1[row] = value;
                    }
                }
            }
        }
        let mut energies = vec![Atom::Zero; edge_count];
        let mut parameters = Vec::with_capacity(support.len());
        let mut rows = self
            .iter_keys()
            .map(|(key, _)| (key.clone(), Vec::with_capacity(support.len())))
            .collect::<Vec<_>>();
        for ((slot, _, _), (basis, values)) in support {
            let parameter = coefficient.call_args([Atom::num(scope), Atom::num(parameters.len())]);
            energies[slot] += &parameter * basis;
            parameters.push(parameter);
            for ((_, arguments), value) in rows.iter_mut().zip(values) {
                arguments.push(value);
            }
        }
        let numerator = numerator.replace_multiple(energy_map_replacements_gs(energies, graph));
        let mut used = vec![false; parameters.len()];
        let _ = numerator.replace_map(|view, _, _| {
            if let AtomView::Fun(function) = view
                && function.get_symbol() == coefficient
                && function.get_nargs() == 2
                && usize::try_from(function.get(0)).ok() == Some(scope)
                && let Ok(index) = usize::try_from(function.get(1))
                && let Some(used) = used.get_mut(index)
            {
                *used = true;
            }
        });
        parameters = parameters
            .into_iter()
            .zip(&used)
            .filter_map(|(parameter, used)| used.then_some(parameter))
            .collect();
        for (_, arguments) in &mut rows {
            *arguments = std::mem::take(arguments)
                .into_iter()
                .zip(&used)
                .filter_map(|(argument, used)| used.then_some(argument))
                .collect();
        }
        Ok((numerator, parameters, rows))
    }

    /// Build the multiplicative identity on this branch family's exact cut
    /// support. A selected raised residue need not contain the uncut root key.
    pub(crate) fn identity_integrands(&self) -> Integrands {
        self.0
            .first()
            .expect("direct residue branches are never empty")
            .1
            .iter()
            .map(|(index, _)| (*index, Atom::one()))
            .collect()
    }

    #[cfg(test)]
    pub(crate) fn factorized_sum(&self) -> Atom {
        self.0
            .iter()
            .map(|(_, integrands)| {
                integrands
                    .resolved()
                    .expect("diagnostic branch numerator definitions must resolve")
                    .iter()
                    .map(|(_, atom)| atom)
                    .sum::<Atom>()
            })
            .sum()
    }

    pub(crate) fn zip_add(&self, other: &Self) -> Result<Self> {
        Self::from_keyed(
            self.iter_keys()
                .map(|(key, integrands)| (key.clone(), integrands.clone()))
                .chain(
                    other
                        .iter_keys()
                        .map(|(key, integrands)| (key.clone(), integrands.clone())),
                ),
        )
    }

    pub(crate) fn zip_mul_unmapped(&self, other: &Integrands) -> Result<Self> {
        Self::from_keyed(
            self.iter_keys()
                .map(|(key, integrands)| Ok((key.clone(), integrands.zip_mul(other)?)))
                .collect::<Result<Vec<_>>>()?,
        )
    }

    pub(crate) fn map(&self, mut map: impl FnMut(&Atom) -> Atom) -> Self {
        Self(
            self.0
                .iter()
                .map(|(key, integrands)| (key.clone(), integrands.map(&mut map)))
                .collect(),
        )
    }

    pub(crate) fn fallible_map(
        &self,
        mut map: impl FnMut(&DirectResidueKey, &Atom) -> Result<Atom>,
    ) -> Result<Self> {
        Self::from_keyed(
            self.iter_keys()
                .map(|(key, integrands)| {
                    Ok((key.clone(), integrands.fallible_map(|atom| map(key, atom))?))
                })
                .collect::<Result<Vec<_>>>()?,
        )
    }

    /// Transform each shared definition once, preserving its formal binding.
    pub(crate) fn map_numerators(
        &self,
        mut map: impl FnMut(&FnMapEntry) -> Result<Atom>,
    ) -> Result<Self> {
        let mut definitions: BTreeMap<Atom, Arc<FnMapEntry>> = BTreeMap::new();
        let mut transformed = BTreeMap::new();
        for (_, integrands) in self.iter_keys() {
            for entry in integrands.numerators() {
                if let Some(existing) = definitions.get(&entry.lhs) {
                    ensure!(
                        existing == entry,
                        "conflicting retained numerator {}",
                        entry.lhs
                    );
                    continue;
                }
                let rhs = map(entry)?;
                transformed.insert(
                    entry.lhs.clone(),
                    if rhs == entry.rhs {
                        Arc::clone(entry)
                    } else {
                        Arc::new(FnMapEntry {
                            rhs,
                            ..entry.as_ref().clone()
                        })
                    },
                );
                definitions.insert(entry.lhs.clone(), Arc::clone(entry));
            }
        }
        Self::from_keyed(
            self.iter_keys()
                .map(|(key, integrands)| {
                    Ok((
                        key.clone(),
                        integrands.clone().with_numerators(
                            integrands
                                .numerators()
                                .iter()
                                .map(|entry| Arc::clone(&transformed[&entry.lhs])),
                        )?,
                    ))
                })
                .collect::<Result<Vec<_>>>()?,
        )
    }

    /// Apply a substitution to both scalar roots and their shared tensor bodies.
    pub(crate) fn map_expressions(
        &self,
        mut map: impl FnMut(&Atom) -> Result<Atom>,
    ) -> Result<Self> {
        self.map_numerators(|entry| map(&entry.rhs))?
            .fallible_map(|_, atom| map(atom))
    }

    pub(crate) fn multiply_key_mapped(
        &self,
        orientation: OrientationProjection<'_>,
        graph: &Graph,
        factor: &Atom,
        scope: (usize, Atom),
    ) -> Result<Self> {
        if matches!(factor.as_view(), AtomView::Num(_)) || factor.is_zero() {
            return Ok(self.map(|atom| atom * factor));
        }
        let (rhs, parameters, rows) =
            self.prepare_numerator(orientation, graph, factor, scope.0)?;
        let lhs = symbol!("gammalooprs::uv::numerator_family")
            .call_args(std::iter::once(scope.1.clone()).chain(parameters.iter().cloned()));
        let entry = Arc::new(FnMapEntry {
            lhs,
            rhs,
            args: parameters
                .iter()
                .cloned()
                .map(Indeterminate::try_from)
                .collect::<std::result::Result<Vec<_>, _>>()
                .map_err(|error| eyre!(error))?,
            tags: vec![scope.1],
        });
        Self::from_keyed(
            self.iter_keys()
                .zip(rows)
                .map(|((key, integrands), (row_key, arguments))| {
                    // The complete branch key owns one energy substitution;
                    // every selected cut order reuses that same mapped factor.
                    ensure!(*key == row_key, "prepared numerator residue order changed");
                    let call =
                        entry
                            .lhs
                            .replace_multiple(parameters.iter().zip(arguments).map(
                                |(parameter, value)| {
                                    Replacement::new(parameter.to_pattern(), value)
                                },
                            ));
                    Ok((
                        key.clone(),
                        integrands.map(|atom| atom * &call).with_numerators(
                            integrands
                                .numerators()
                                .iter()
                                .cloned()
                                .chain([Arc::clone(&entry)]),
                        )?,
                    ))
                })
                .collect::<Result<Vec<_>>>()?,
        )
    }

    /// Expand shared bodies once, retaining formal coefficient families until
    /// branch specialization. Native finite series supply conservative lower
    /// bounds; a vanishing probe is not a claim of an identically zero body.
    pub(crate) fn series_preserving_numerators(
        &self,
        variable: Symbol,
        center: AtomView<'_>,
        depth: i64,
        scope: Atom,
    ) -> Result<Self> {
        let all_definitions = self.0[0].1.clone().with_numerators(
            self.iter_keys()
                .flat_map(|(_, integrands)| integrands.numerators().iter().cloned()),
        )?;
        let definitions = all_definitions.numerators();
        if definitions.is_empty() {
            return self.fallible_map(|_, atom| {
                atom.series_preserving_factors(variable, center, depth, &[])
            });
        }
        let family = symbol!("gammalooprs::uv::numerator_family");
        let delta = Atom::var(variable) - center;
        // Match the polynomial owners already recognized by EnergyPowerAnalyzer.
        // Inert projectors keep their original attributes outside this operation.
        let multilinear = |head: Symbol, arity: usize| {
            head.is_linear()
                || (head == GS.dot && arity == 2)
                || [*CYCLIC, *SYM, *ANTISYM].contains(&head)
        };
        let mut linear_heads = BTreeMap::<Symbol, Symbol>::new();
        let mut probes = Vec::with_capacity(definitions.len());
        let mut bodies = Vec::with_capacity(definitions.len());
        for (index, entry) in definitions.iter().enumerate() {
            let parameters = entry
                .args
                .iter()
                .cloned()
                .map(Atom::from)
                .collect::<Vec<_>>();
            ensure!(
                !entry.tags.is_empty()
                    && entry.lhs
                        == family.call_args(
                            entry.tags.iter().cloned().chain(parameters.iter().cloned())
                        )
                    && !entry.lhs.contains_symbol(variable)
                    && definitions[..index]
                        .iter()
                        .all(|other| other.tags != entry.tags)
                    && parameters
                        .iter()
                        .enumerate()
                        .all(|(i, parameter)| !parameters[..i].contains(parameter)),
                "invalid or ambiguous retained numerator binding {}",
                entry.lhs
            );
            // Specialization may remove leading coefficients, but must never
            // introduce a pole in the Taylor variable through a map parameter.
            let contains_parameter =
                |atom: AtomView<'_>| parameters.iter().any(|parameter| atom.contains(parameter));
            let mut invalid = None;
            let mut inert_heads = std::collections::BTreeSet::new();
            entry.rhs.visitor(&mut |atom| {
                if invalid.is_some()
                    || parameters
                        .iter()
                        .any(|parameter| parameter.as_view() == atom)
                {
                    return false;
                }
                let regular = match atom {
                    AtomView::Fun(call) => {
                        let known = multilinear(call.get_symbol(), call.get_nargs());
                        if known && !call.get_symbol().is_linear() && atom.contains_symbol(variable)
                        {
                            inert_heads.insert(call.get_symbol());
                        }
                        !contains_parameter(atom) || known
                    }
                    AtomView::Pow(power) if contains_parameter(atom) => {
                        Rational::try_from(power.get_exp())
                            .is_ok_and(|power| power.is_integer() && !power.is_negative())
                    }
                    _ => true,
                };
                if !regular {
                    invalid = Some(atom.to_owned());
                }
                regular
            });
            ensure!(
                invalid.is_none(),
                "retained numerator must be polynomial in map parameters: {:?}",
                invalid
            );
            for head in inert_heads {
                if linear_heads.contains_key(&head) {
                    continue;
                }
                let linear = loop {
                    let id = NUMERATOR_SCOPE.fetch_add(1, Ordering::Relaxed);
                    let name = format!("gammalooprs::uv::numerator_series_linear_{id}");
                    if symbolica::state::State::symbol_iter().any(|(symbol, existing)| {
                        existing == name || symbol.get_aliases().iter().any(|alias| alias == &name)
                    }) {
                        continue;
                    }
                    let mut attributes = head.get_attributes();
                    attributes.push(SymbolAttribute::Linear);
                    break SymbolBuilder::new(symbolica::wrap_symbol!(name.as_str()))
                        .with_attributes(attributes)
                        .build()
                        .map_err(|error| eyre!("{error}"))?;
                };
                linear_heads.insert(head, linear);
            }
            // Normalize only multilinear argument slots, never distribute the
            // surrounding numerator product or expand projector permutations.
            let normalized = entry.rhs.replace_map_bottom_up(|atom, _, out| {
                if let AtomView::Fun(call) = atom
                    && atom.contains_symbol(variable)
                    && let Some(linear) = linear_heads.get(&call.get_symbol())
                {
                    **out = linear.call_args(call.iter());
                }
            });
            let normalized = normalized.replace_map_bottom_up(|atom, _, out| {
                if let AtomView::Fun(call) = atom
                    && let Some((original, _)) = linear_heads
                        .iter()
                        .find(|(_, linear)| **linear == call.get_symbol())
                {
                    **out = original.call_args(call.iter());
                }
            });
            let mut hidden_taylor = None;
            normalized.visitor(&mut |atom| {
                if let AtomView::Fun(call) = atom
                    && multilinear(call.get_symbol(), call.get_nargs())
                    && atom.contains_symbol(variable)
                {
                    hidden_taylor = Some(atom.to_owned());
                    return false;
                }
                hidden_taylor.is_none()
            });
            ensure!(
                hidden_taylor.is_none(),
                "Taylor dependence remains inside a multilinear numerator call: {:?}",
                hidden_taylor
            );
            ensure!(
                linear_heads
                    .values()
                    .all(|head| !normalized.contains_symbol(*head)),
                "temporary multilinear numerator head escaped normalization"
            );
            // Keep independent tensor/vertex factors outside native coefficient
            // collection, just as series_preserving_factors does for scalar roots.
            // The entire dependent product still owns native Laurent precision.
            let mut independent = Atom::one();
            let mut dependent = normalized.clone();
            if let AtomView::Mul(product) = normalized.as_view() {
                dependent = Atom::one();
                for factor in product.iter() {
                    if factor.contains_symbol(variable) {
                        dependent *= factor;
                    } else {
                        independent *= factor;
                    }
                }
            } else if !normalized.contains_symbol(variable) {
                independent = normalized.clone();
                dependent = Atom::one();
            }
            let probe = dependent.series(variable, center, SeriesDepth::absolute(0))?;
            let bound = if probe.is_zero() {
                probe.absolute_order()
            } else {
                probe.absolute_order() - probe.relative_order()
            };
            probes.push(
                FnMapEntry {
                    rhs: if normalized.is_zero() {
                        Atom::Zero
                    } else {
                        &entry.lhs * delta.pow(Atom::num(bound.clone()))
                    },
                    ..entry.as_ref().clone()
                }
                .replacement(),
            );
            bodies.push((independent, dependent, bound));
        }

        let mut extra = Rational::from(0);
        let mut root_bounds = BTreeMap::new();
        for (branch_index, (_, integrands)) in self.iter_keys().enumerate() {
            for (cut, atom) in integrands.iter() {
                let mut invalid = None;
                atom.visitor(&mut |view| {
                    if invalid.is_some() {
                        return false;
                    }
                    let regular = match view {
                        AtomView::Fun(call) if call.get_symbol() == family => {
                            let entry = integrands.numerators().iter().find(|entry| {
                                call.get_nargs() == entry.tags.len() + entry.args.len()
                                    && entry
                                        .tags
                                        .iter()
                                        .enumerate()
                                        .all(|(i, tag)| call.get(i) == tag.as_view())
                            });
                            let valid = entry.is_some_and(|entry| {
                                (entry.tags.len()..call.get_nargs())
                                    .all(|i| Rational::try_from(call.get(i)).is_ok())
                                    && !view.contains_symbol(variable)
                            });
                            if !valid {
                                invalid = Some(view.to_owned());
                            }
                            return false;
                        }
                        AtomView::Fun(_) => !view.contains_symbol(family),
                        AtomView::Pow(power) if view.contains_symbol(family) => {
                            Rational::try_from(power.get_exp())
                                .is_ok_and(|power| power.is_integer() && !power.is_negative())
                        }
                        _ => true,
                    };
                    if !regular {
                        invalid = Some(view.to_owned());
                    }
                    regular
                });
                ensure!(
                    invalid.is_none(),
                    "invalid retained numerator occurrence (expected numeric arguments and polynomial use): {:?}",
                    invalid
                );
                let proxy = atom.replace_multiple(&probes);
                let probe = proxy.series(variable, center, SeriesDepth::absolute(0))?;
                let bound = if probe.is_zero() {
                    probe.absolute_order()
                } else {
                    probe.absolute_order() - probe.relative_order()
                };
                extra = extra.max(Rational::from(depth) - &bound);
                root_bounds.insert((branch_index, *cut), bound);
            }
        }

        let mut replacements = Vec::with_capacity(definitions.len());
        let mut coefficients = BTreeMap::new();
        for (index, (entry, (independent, dependent, bound))) in
            definitions.iter().zip(bodies).enumerate()
        {
            let series =
                dependent.series(variable, center, SeriesDepth::absolute(bound + &extra))?;
            let mut retained = Vec::new();
            let mut rhs = Atom::Zero;
            for (exponent, coefficient) in series.terms() {
                let coefficient = &independent * coefficient;
                if coefficient.is_zero() {
                    continue;
                }
                let tag = symbol!("gammalooprs::uv::numerator_series_tag").call_args([
                    scope.clone(),
                    Atom::num(index),
                    Atom::num(exponent.clone()),
                ]);
                let lhs = family.call_args(
                    std::iter::once(tag.clone()).chain(entry.args.iter().cloned().map(Atom::from)),
                );
                rhs += &lhs * delta.pow(Atom::num(exponent));
                retained.push(Arc::new(FnMapEntry {
                    lhs,
                    rhs: coefficient,
                    args: entry.args.clone(),
                    tags: vec![tag],
                }));
            }
            replacements.push(
                FnMapEntry {
                    rhs,
                    ..entry.as_ref().clone()
                }
                .replacement(),
            );
            coefficients.insert(entry.lhs.clone(), retained);
        }
        crate::debug_tags!(#generation, #profile, #uv, #numerator;
            stage = "shared_numerator_taylor_coefficients",
            source_families = definitions.len(),
            residue_rows = self.0.len(),
            coefficient_families = coefficients.values().map(Vec::len).sum::<usize>(),
            source_body_bytes = definitions.iter().map(|entry| entry.rhs.as_view().get_byte_size()).sum::<usize>(),
            coefficient_body_bytes = coefficients.values().flatten().map(|entry| entry.rhs.as_view().get_byte_size()).sum::<usize>(),
            "Retained shared Taylor coefficients before residue specialization"
        );
        Self::from_keyed(
            self.iter_keys()
                .enumerate()
                .map(|(branch_index, (key, integrands))| {
                    let roots = integrands
                        .iter()
                        .map(|(cut, atom)| {
                            let series = if root_bounds[&(branch_index, *cut)] > depth {
                                Atom::Zero
                            } else {
                                atom.replace_multiple(&replacements)
                                    .series_preserving_factors(variable, center, depth, &[])?
                            };
                            Ok((*cut, series))
                        })
                        .collect::<Result<Integrands>>()?;
                    Ok((
                        key.clone(),
                        roots.with_numerators(
                            integrands
                                .numerators()
                                .iter()
                                .flat_map(|entry| coefficients[&entry.lhs].iter().cloned()),
                        )?,
                    ))
                })
                .collect::<Result<Vec<_>>>()?,
        )
    }

    /// Materialize residue selectors only at the evaluator boundary.
    pub(crate) fn materialize(&self, materialize_key_selector: bool) -> Result<Integrands> {
        let mut branches = self.iter_keys();
        let (first_key, first) = branches
            .next()
            .ok_or_else(|| eyre!("direct local-3D residue branches cannot be empty"))?;
        let first_selector = first_key.selector(materialize_key_selector);
        first
            .map(|atom| atom * &first_selector)
            .zip_add(branches.map(|(key, integrands)| {
                let selector = key.selector(materialize_key_selector);
                integrands.map(|atom| atom * &selector)
            }))
    }
}

impl Neg for DirectResidueBranches {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(
            self.0
                .into_iter()
                .map(|(key, integrands)| (key, -integrands))
                .collect(),
        )
    }
}

#[cfg(test)]
mod tests {
    use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
    use spenso::structure::representation::{LibraryRep, Minkowski, RepName};
    use symbolica::{domains::rational::Rational, function, id::Replacement};
    use typed_index_collections::TiVec;

    use crate::{
        cff::{
            CutCFFIndex,
            expression::{OrientationData, OrientationExpression},
        },
        dot,
        graph::{Graph, parse::IntoGraph},
        initialisation::test_initialise,
        settings::global::OrientationPattern,
        utils::GS,
    };

    use super::*;

    #[test]
    fn shared_numerator_series_normalizes_known_multilinear_owners() -> Result<()> {
        test_initialise()?;
        let variable = symbol!("shared_multilinear_test::t"; Scalar);
        let s = symbol!("shared_multilinear_test::s"; Scalar);
        let z = symbol!("shared_multilinear_test::z"; Scalar);
        let family = symbol!("gammalooprs::uv::numerator_family");
        let cut = CutCFFIndex::new_all_none();
        let vector_rep = function!(LibraryRep::from(Minkowski {}).symbol(), 4);
        let vectors = [0, 1, 2, 3].map(|edge| GS.emr_mom(EdgeIndex(edge), &vector_rep));
        let start = idenso::bis!(4, Atom::var(symbol!("shared_multilinear_test::start")));
        let end = idenso::bis!(4, Atom::var(symbol!("shared_multilinear_test::end")));
        let bis = function!(
            LibraryRep::from(idenso::representations::Bispinor {}).symbol(),
            4
        );
        let linear = symbol!("shared_multilinear_test::linear"; Linear);
        for (head, closed) in [
            (GS.dot, false),
            (*CYCLIC, false),
            (*SYM, false),
            (*ANTISYM, false),
            (linear, false),
            (*SYM, true),
        ] {
            let arguments = if head == GS.dot {
                vectors.clone()
            } else {
                vectors.each_ref().map(|vector| idenso::gamma!(vector))
            };
            let [a, b, c, d] = &arguments;
            let apply = |left: &Atom, right: &Atom| {
                let operator = function!(head, left, right);
                if head == GS.dot {
                    operator
                } else if closed {
                    spenso::trace!(&bis, operator)
                } else {
                    spenso::chain!(&start, &end, operator)
                }
            };
            let original_attributes = head.get_attributes();
            for center in [Atom::Zero, Atom::num(2)] {
                let delta = Atom::var(variable) - &center;
                let body = apply(
                    &(Atom::var(s) * a + &delta * b),
                    &(Atom::var(z) * c + &delta * d),
                );
                let tag = DirectResidueBranches::numerator_scope().1;
                let definition = Arc::new(FnMapEntry {
                    lhs: function!(family, &tag, s, z),
                    rhs: body,
                    args: vec![s.into(), z.into()],
                    tags: vec![tag.clone()],
                });
                let rows = [[1, 2], [-1, 0], [0, 0]];
                let branches = DirectResidueBranches::from_keyed(
                    rows.into_iter()
                        .enumerate()
                        .map(|(row, [sign, node])| {
                            let call = function!(family, &tag, sign, node);
                            let root = call / (delta.pow(2) * (Atom::one() - &delta));
                            Ok((
                                DirectResidueKey::source(
                                    OrientationID(0),
                                    vec![LinearEnergyExpr::uniform_scale(row as i64)],
                                ),
                                Integrands::from_iter([(cut, root)])
                                    .with_numerators([Arc::clone(&definition)])?,
                            ))
                        })
                        .collect::<Result<Vec<_>>>()?,
                )?;
                let first = branches.series_preserving_numerators(
                    variable,
                    center.as_view(),
                    0,
                    DirectResidueBranches::numerator_scope().1,
                )?;
                let repeated = first.series_preserving_numerators(
                    variable,
                    center.as_view(),
                    0,
                    DirectResidueBranches::numerator_scope().1,
                )?;
                assert_eq!(first.0[0].1.numerators().len(), 3);
                assert_eq!(repeated.0[0].1.numerators().len(), 3);
                for expanded in [&first, &repeated] {
                    assert_eq!(expanded.0.len(), rows.len());
                    for ((_, integrands), [sign, node]) in expanded.iter_keys().zip(rows) {
                        for entry in integrands.numerators() {
                            assert!(!entry.rhs.contains_symbol(Symbol::DERIVATIVE));
                            assert!(entry.rhs.contains_symbol(head));
                            entry.rhs.visitor(&mut |atom| {
                                if let AtomView::Fun(call) = atom {
                                    assert!(
                                        !call
                                            .get_symbol()
                                            .get_stripped_name()
                                            .starts_with("numerator_series_linear_")
                                    );
                                }
                                true
                            });
                        }
                        // This oracle uses the owner's explicit bilinear identity,
                        // not the generic native derivative of an inert head.
                        let polynomial = Atom::num(sign * node) * apply(a, c)
                            + &delta
                                * (Atom::num(sign) * apply(a, d) + Atom::num(node) * apply(b, c))
                            + delta.pow(2) * apply(b, d);
                        let expected = (polynomial / (delta.pow(2) * (Atom::one() - &delta)))
                            .series(variable, &center, SeriesDepth::absolute(0))?
                            .to_atom();
                        let actual = integrands.resolved()?;
                        assert!(
                            (actual.iter().next().unwrap().1 - &expected)
                                .expand()
                                .is_zero(),
                            "multilinear owner {head}, center {center}, row [{sign}, {node}]"
                        );
                    }
                    for (_, integrands) in expanded.iter_keys().skip(1) {
                        for (first, next) in expanded.0[0]
                            .1
                            .numerators()
                            .iter()
                            .zip(integrands.numerators())
                        {
                            assert!(Arc::ptr_eq(first, next));
                        }
                    }
                }
            }
            assert_eq!(
                head.get_attributes(),
                original_attributes,
                "normalization must not mutate physical owner attributes"
            );
        }
        Ok(())
    }

    #[test]
    fn shared_numerator_series_does_not_capture_multilinear_scope_names() -> Result<()> {
        test_initialise()?;
        let tag = DirectResidueBranches::numerator_scope().1;
        let series_scope = DirectResidueBranches::numerator_scope().1;
        let next = NUMERATOR_SCOPE.load(Ordering::Relaxed);
        let occupied_name = format!("gammalooprs::uv::numerator_series_linear_{next}");
        let occupied = symbol!(occupied_name.as_str(); Scalar);
        let alias = format!("gammalooprs::uv::numerator_series_linear_{}", next + 1);
        let owner_name = format!("gammalooprs::uv::multilinear_collision_owner_{next}");
        let alias_owner = SymbolBuilder::new(symbolica::wrap_symbol!(owner_name.as_str()))
            .with_aliases([alias])
            .build()
            .map_err(|error| eyre!("{error}"))?;
        let variable = symbol!("shared_multilinear_collision_test::t"; Scalar);
        let parameter = symbol!("shared_multilinear_collision_test::p"; Scalar);
        let a = Atom::var(symbol!("shared_multilinear_collision_test::a"));
        let b = Atom::var(symbol!("shared_multilinear_collision_test::b"));
        let c = Atom::var(symbol!("shared_multilinear_collision_test::c"));
        let family = symbol!("gammalooprs::uv::numerator_family");
        let spectator = Atom::var(occupied) * Atom::var(alias_owner);
        let body = &spectator
            * function!(
                *SYM,
                Atom::var(parameter) * &a + Atom::var(variable) * &b,
                &c
            );
        let definition = Arc::new(FnMapEntry {
            lhs: function!(family, &tag, parameter),
            rhs: body,
            args: vec![parameter.into()],
            tags: vec![tag.clone()],
        });
        let cut = CutCFFIndex::new_all_none();
        let branches = DirectResidueBranches::production(
            OrientationID(0),
            Integrands::from_iter([(cut, function!(family, &tag, 2) / Atom::var(variable))])
                .with_numerators([definition])?,
        )?;
        let expanded = branches.series_preserving_numerators(
            variable,
            Atom::Zero.as_view(),
            0,
            series_scope,
        )?;
        let actual = expanded.0[0].1.resolved()?;
        let expected = &spectator
            * (Atom::num(2) * function!(*SYM, &a, &c) / Atom::var(variable)
                + function!(*SYM, &b, &c));
        assert!(
            (actual.iter().next().unwrap().1 - expected)
                .expand()
                .is_zero()
        );
        assert!(
            expanded.0[0]
                .1
                .numerators()
                .iter()
                .all(|entry| entry.rhs.contains_symbol(occupied)
                    && entry.rhs.contains_symbol(alias_owner))
        );
        Ok(())
    }

    #[test]
    fn shared_numerator_series_preserves_independent_tensor_factors() -> Result<()> {
        test_initialise()?;
        let variable = symbol!("shared_series_factor_test::t");
        let [a, b, c, d] = ["a", "b", "c", "d"].map(|name| {
            Atom::var(symbol!(
                format!("shared_series_factor_test::{name}").as_str()
            ))
        });
        let index = LibraryRep::from(Minkowski {})
            .to_symbolic([Atom::var(symbol!("shared_series_factor_test::mu"))]);
        let g = (GS.emr_vec_index(EdgeIndex(0), &index) + GS.emr_vec_index(EdgeIndex(1), &index))
            * (Atom::var(symbol!("shared_series_factor_test::u"))
                + Atom::var(symbol!("shared_series_factor_test::v")));
        let family = symbol!("gammalooprs::uv::numerator_family");
        let cut = CutCFFIndex::new_all_none();
        let expected_coefficients = [&g * &a * &c, &g * (&a * &d + &b * &c), &g * &b * &d]
            .into_iter()
            .collect::<std::collections::BTreeSet<_>>();
        for center in [Atom::Zero, Atom::num(2)] {
            let delta = Atom::var(variable) - &center;
            let tag = DirectResidueBranches::numerator_scope().1;
            let call = family.call_args([tag.clone()]);
            let body = &g * (&a + &b * &delta) * (&c + &d * &delta);
            let definition = Arc::new(FnMapEntry {
                lhs: call.clone(),
                rhs: body.clone(),
                args: vec![],
                tags: vec![tag],
            });
            let source = DirectResidueBranches::production(
                OrientationID(0),
                Integrands::from_iter([(cut, &call / delta.pow(2))])
                    .with_numerators([Arc::clone(&definition)])?,
            )?;
            let mut expanded = source.series_preserving_numerators(
                variable,
                center.as_view(),
                0,
                DirectResidueBranches::numerator_scope().1,
            )?;
            for round in 0..2 {
                let integrands = &expanded.0[0].1;
                assert_eq!(
                    integrands
                        .numerators()
                        .iter()
                        .map(|entry| entry.rhs.clone())
                        .collect::<std::collections::BTreeSet<_>>(),
                    expected_coefficients,
                    "independent tensor sums must remain factored in every shared coefficient",
                );
                let expected = (&body / delta.pow(2))
                    .series(variable, &center, SeriesDepth::absolute(0))?
                    .to_atom();
                let actual = integrands.resolved()?;
                assert!(
                    (actual.iter().next().unwrap().1 - &expected)
                        .expand()
                        .is_zero()
                );
                if round == 0 {
                    expanded = expanded.series_preserving_numerators(
                        variable,
                        center.as_view(),
                        0,
                        DirectResidueBranches::numerator_scope().1,
                    )?;
                }
            }

            let zero_body = &g * ((&a + &b * &delta) / &delta - &a / &delta - &b);
            let cancelled = DirectResidueBranches::production(
                OrientationID(0),
                Integrands::from_iter([(cut, &call / delta.pow(2))]).with_numerators([
                    Arc::new(FnMapEntry {
                        rhs: zero_body,
                        ..definition.as_ref().clone()
                    }),
                ])?,
            )?
            .series_preserving_numerators(
                variable,
                center.as_view(),
                0,
                DirectResidueBranches::numerator_scope().1,
            )?;
            assert!(cancelled.0[0].1.numerators().is_empty());
            assert!(cancelled.0[0].1.iter().all(|(_, atom)| atom.is_zero()));
        }
        Ok(())
    }

    #[test]
    fn shared_numerator_series_matches_specialized_laurent_products() -> Result<()> {
        test_initialise()?;
        let variable = symbol!("direct_series_test::t");
        let parameters = ["sigma", "node", "tau", "other_node"]
            .map(|name| Atom::var(symbol!(format!("direct_series_test::{name}").as_str())));
        let family = symbol!("gammalooprs::uv::numerator_family");
        let cut = CutCFFIndex::new_all_none();
        let raised = CutCFFIndex {
            lu_cut_order: Some(1),
            ..cut
        };
        for center in [Atom::Zero, Atom::num(2)] {
            let delta = Atom::var(variable) - &center;
            // Independent sign/node pairs include a zero leading coefficient.
            // Both OSE square roots must contribute their Taylor derivatives.
            let body = (&parameters[0] * (Atom::one() + &delta).pow(Atom::num((1, 2)))
                + &parameters[1])
                * (&parameters[2] * (Atom::num(4) + Atom::num(2) * &delta).pow(Atom::num((1, 2)))
                    + &parameters[3]);
            let tags = [
                DirectResidueBranches::numerator_scope().1,
                DirectResidueBranches::numerator_scope().1,
            ];
            let entries = [body, delta.pow(3) * (Atom::one() + &delta)]
                .into_iter()
                .enumerate()
                .map(|(index, rhs)| {
                    Arc::new(FnMapEntry {
                        lhs: family.call_args(
                            std::iter::once(tags[index].clone()).chain(parameters.iter().cloned()),
                        ),
                        rhs,
                        args: parameters
                            .iter()
                            .cloned()
                            .map(Indeterminate::try_from)
                            .collect::<std::result::Result<_, _>>()
                            .unwrap(),
                        tags: vec![tags[index].clone()],
                    })
                })
                .collect::<Vec<_>>();
            let branches = DirectResidueBranches::from_keyed(
                [[1, -1, -1, 2], [-1, 2, 1, 0], [1, 3, -1, -2]]
                    .into_iter()
                    .enumerate()
                    .map(|(row, values)| {
                        let calls = tags.each_ref().map(|tag| {
                            family.call_args(
                                std::iter::once(tag.clone()).chain(values.map(Atom::num)),
                            )
                        });
                        let roots = Integrands::from_iter([
                            (cut, &calls[0] / (delta.pow(3) * (Atom::one() - &delta))),
                            (
                                raised,
                                &calls[0] * &calls[1] / delta.pow(4)
                                    + &calls[0]
                                        * ((Atom::one() - &delta).pow(-1) - Atom::one() - &delta),
                            ),
                        ])
                        .with_numerators(entries.iter().cloned())?;
                        Ok((
                            DirectResidueKey::source(
                                OrientationID(0),
                                vec![LinearEnergyExpr::uniform_scale(row as i64)],
                            ),
                            roots,
                        ))
                    })
                    .collect::<Result<Vec<_>>>()?,
            )?;
            for depth in [-1, 0, 2] {
                let expanded = branches.series_preserving_numerators(
                    variable,
                    center.as_view(),
                    depth,
                    DirectResidueBranches::numerator_scope().1,
                )?;
                assert_eq!(
                    expanded.0.len(),
                    branches.0.len(),
                    "same-host source maps must stay distinct"
                );
                for ((key, actual), (expected_key, original)) in
                    expanded.iter_keys().zip(branches.iter_keys())
                {
                    assert_eq!(key, expected_key);
                    let oracle = original.resolved()?.fallible_map(|atom| {
                        Ok(atom
                            .series(variable, &center, SeriesDepth::absolute(depth))?
                            .to_atom())
                    })?;
                    let difference = actual
                        .resolved()?
                        .checked_zip(&oracle, |_, actual, expected| {
                            Ok((actual - expected).expand())
                        })?;
                    assert!(
                        difference.iter().all(|(_, atom)| atom.is_zero()),
                        "center {center}, depth {depth}: {difference:?}"
                    );
                }
                let single = DirectResidueBranches::from_keyed([branches.0[0].clone()])?
                    .series_preserving_numerators(
                        variable,
                        center.as_view(),
                        depth,
                        DirectResidueBranches::numerator_scope().1,
                    )?;
                assert_eq!(
                    expanded.0[0].1.numerators().len(),
                    single.0[0].1.numerators().len(),
                    "definition count must not grow with branch count"
                );
                for (_, integrands) in expanded.iter_keys().skip(1) {
                    for (first, shared) in expanded.0[0]
                        .1
                        .numerators()
                        .iter()
                        .zip(integrands.numerators())
                    {
                        assert!(Arc::ptr_eq(first, shared));
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn shared_numerator_series_handles_finite_zero_probes_and_cancellation() -> Result<()> {
        test_initialise()?;
        let variable = symbol!("direct_zero_series_test::t");
        let t = Atom::var(variable);
        let a = Atom::var(symbol!("direct_zero_series_test::a"));
        let b = Atom::var(symbol!("direct_zero_series_test::b"));
        let family = symbol!("gammalooprs::uv::numerator_family");
        let cut = CutCFFIndex::new_all_none();
        for body in [
            Atom::Zero,
            (&a + &b * &t) / &t - &a / &t - &b,
            t.pow(5) * (&a + &t),
        ] {
            let tag = DirectResidueBranches::numerator_scope().1;
            let call = family.call_args([tag.clone()]);
            let entry = Arc::new(FnMapEntry {
                lhs: call.clone(),
                rhs: body,
                args: vec![],
                tags: vec![tag],
            });
            for root in [
                &call / t.pow(7),
                &call * ((&a + &b * &t) / &t - &a / &t - &b),
            ] {
                let original =
                    Integrands::from_iter([(cut, root)]).with_numerators([Arc::clone(&entry)])?;
                let branches =
                    DirectResidueBranches::production(OrientationID(0), original.clone())?;
                let expanded = branches.series_preserving_numerators(
                    variable,
                    Atom::Zero.as_view(),
                    0,
                    DirectResidueBranches::numerator_scope().1,
                )?;
                let expected = original.resolved()?.fallible_map(|atom| {
                    Ok(atom
                        .series(variable, Atom::Zero, SeriesDepth::absolute(0))?
                        .to_atom())
                })?;
                let difference = expanded.0[0]
                    .1
                    .resolved()?
                    .checked_zip(&expected, |_, a, b| Ok((a - b).expand()))?;
                assert!(difference.iter().all(|(_, atom)| atom.is_zero()));
            }
        }
        Ok(())
    }

    #[test]
    fn shared_numerator_series_rejects_singular_specialization_and_hidden_calls() -> Result<()> {
        test_initialise()?;
        let variable = symbol!("direct_invalid_series_test::t");
        let parameter = symbol!("direct_invalid_series_test::z");
        let t = Atom::var(variable);
        let z = Atom::var(parameter);
        let family = symbol!("gammalooprs::uv::numerator_family");
        let outer = symbol!("direct_invalid_series_test::outer");
        let tag = DirectResidueBranches::numerator_scope().1;
        let lhs = family.call_args([tag.clone(), z.clone()]);
        let call = family.call_args([tag.clone(), Atom::Zero]);
        let cut = CutCFFIndex::new_all_none();
        let cases = [
            ((&z + &t).pow(-1), call.clone()),
            ((&z + &t).pow(Atom::num((1, 2))), call.clone()),
            (outer.call_args([z.clone()]), call.clone()),
            (&z + &t, call.pow(-1)),
            (&z + &t, outer.call_args([call.clone()])),
            (&z + &t, family.call_args([tag.clone(), t.clone()])),
        ];
        for (rhs, root) in cases {
            let entry = Arc::new(FnMapEntry {
                lhs: lhs.clone(),
                rhs,
                args: vec![parameter.into()],
                tags: vec![tag.clone()],
            });
            let branches = DirectResidueBranches::production(
                OrientationID(0),
                Integrands::from_iter([(cut, root)]).with_numerators([entry])?,
            )?;
            assert!(
                branches
                    .series_preserving_numerators(
                        variable,
                        Atom::Zero.as_view(),
                        0,
                        DirectResidueBranches::numerator_scope().1
                    )
                    .is_err()
            );
        }
        Ok(())
    }

    #[test]
    fn prepared_numerator_keeps_full_affine_rows_and_independent_sampling_nodes() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(
            digraph G {
                edge [particle="scalar_1"];
                node [num=1];
                a -> b [id=0];
                a -> b [id=1];
            },
            "scalars"
        )?;
        let production = vec![OrientationExpression {
            data: OrientationData::new(EdgeVec::from_iter([
                Orientation::Default,
                Orientation::Reversed,
            ])),
            loop_energy_map: Vec::new(),
            edge_energy_map: vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::ose(EdgeIndex(1), -1),
            ],
            variants: Vec::new(),
        }]
        .into_iter()
        .collect::<TiVec<OrientationID, _>>();
        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::default();
        let orientation = OrientationProjection::exact(&production, &options, &pattern, false);
        let host = OrientationID(0);
        let keys = [
            DirectResidueKey::production(host),
            DirectResidueKey::source(
                host,
                vec![
                    LinearEnergyExpr {
                        internal_terms: vec![(EdgeIndex(0), 1.into())],
                        uniform_scale_coeff: 2.into(),
                        ..LinearEnergyExpr::zero()
                    },
                    LinearEnergyExpr::zero(),
                ],
            ),
            DirectResidueKey::source(
                host,
                vec![
                    LinearEnergyExpr {
                        internal_terms: vec![
                            (EdgeIndex(0), (-3).into()),
                            (EdgeIndex(1), 3.into()),
                            (EdgeIndex(0), 1.into()),
                        ],
                        external_terms: vec![(EdgeIndex(2), (5, 2).into())],
                        uniform_scale_coeff: (-3).into(),
                        constant: 7.into(),
                    },
                    LinearEnergyExpr::uniform_scale(4),
                ],
            ),
            DirectResidueKey::source(host, vec![LinearEnergyExpr::zero(); 2]),
        ];
        let branches = DirectResidueBranches::from_keyed(keys.iter().cloned().map(|key| {
            (
                key,
                [(CutCFFIndex::new_all_none(), Atom::one())]
                    .into_iter()
                    .collect(),
            )
        }))?;
        let q0 = GS.emr_mom(EdgeIndex(0), GS.cind(0));
        let numerator = (q0.clone() + Atom::var(symbol!("direct_3d_test::a")))
            * (q0 + Atom::var(symbol!("direct_3d_test::b")));
        let (template, parameters, rows) =
            branches.prepare_numerator(orientation, &graph, &numerator, 21)?;

        assert_eq!(
            parameters.len(),
            5,
            "unused second-edge coefficients must be pruned"
        );
        assert_eq!(
            rows.len(),
            keys.len(),
            "equal hosts must not merge different maps"
        );
        assert!(
            matches!(template.as_view(), AtomView::Mul(product) if product.iter().count() == 2)
        );
        assert!(template.contains_symbol(GS.numerator_sampling_scale));
        for ((key, arguments), expected_key) in rows.iter().zip(&keys) {
            assert_eq!(key, expected_key);
            assert_eq!(arguments.len(), parameters.len());
            assert!(
                arguments
                    .iter()
                    .all(|argument| matches!(argument.as_view(), AtomView::Num(_)))
            );
            let specialized = template.replace_multiple(parameters.iter().zip(arguments).map(
                |(parameter, argument)| {
                    Replacement::new(parameter.to_pattern(), argument.to_pattern())
                },
            ));
            assert_eq!(
                specialized,
                key.map_numerator(orientation, &graph, &numerator)?
            );
        }
        assert_eq!(rows[1].1, [1, 0, 0, 2, 0].map(Atom::num));
        assert_eq!(
            rows[2].1,
            [
                Atom::num(-2),
                Atom::num(3),
                Atom::num(Rational::from((5, 2))),
                Atom::num(-3),
                Atom::num(7),
            ]
        );
        assert!(rows[3].1.iter().all(|argument| argument.is_zero()));
        assert!(
            branches
                .prepare_numerator(orientation, &graph, &template, 21)
                .is_err()
        );

        let mismatched = DirectResidueBranches::from_keyed([
            (keys[0].clone(), branches.identity_integrands()),
            (
                DirectResidueKey::source(host, vec![LinearEnergyExpr::zero()]),
                branches.identity_integrands(),
            ),
        ])?;
        assert!(
            mismatched
                .prepare_numerator(orientation, &graph, &numerator, 22)
                .is_err()
        );
        Ok(())
    }

    #[test]
    fn prepared_numerator_reuses_vector_loop_and_external_mapping_rules() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph G {
            edge [num=1 mass=0]
            node [num=1]
            ext [style=invis]
            ext -> a [id=0]
            a -> b [id=1 lmb_id=0]
            a -> b [id=2]
            b -> ext [id=3]
        })?;
        let production = vec![OrientationExpression {
            data: OrientationData::new(EdgeVec::from_iter([
                Orientation::Undirected,
                Orientation::Default,
                Orientation::Reversed,
                Orientation::Undirected,
            ])),
            loop_energy_map: Vec::new(),
            edge_energy_map: vec![
                LinearEnergyExpr::zero(),
                LinearEnergyExpr::ose(EdgeIndex(1), 1),
                LinearEnergyExpr::ose(EdgeIndex(2), -1),
                LinearEnergyExpr::zero(),
            ],
            variants: Vec::new(),
        }]
        .into_iter()
        .collect::<TiVec<OrientationID, _>>();
        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::default();
        let orientation = OrientationProjection::exact(&production, &options, &pattern, false);
        let branches = DirectResidueBranches::production(
            OrientationID(0),
            [(CutCFFIndex::new_all_none(), Atom::one())]
                .into_iter()
                .collect(),
        )?;
        let index =
            LibraryRep::from(Minkowski {}).to_symbolic([Atom::var(symbol!("direct_3d_test::mu"))]);
        let external = GS.emr_mom(EdgeIndex(0), &index);
        let numerator = &external
            * (GS.emr_mom(EdgeIndex(1), &index) + function!(GS.loop_mom, 0, &index))
            * function!(GS.loop_mom, 0, GS.cind(1));
        let (template, parameters, rows) =
            branches.prepare_numerator(orientation, &graph, &numerator, 23)?;
        assert_eq!(parameters.len(), 1);
        let specialized = template.replace_multiple(parameters.iter().zip(&rows[0].1).map(
            |(parameter, argument)| Replacement::new(parameter.to_pattern(), argument.to_pattern()),
        ));
        assert_eq!(
            specialized,
            rows[0].0.map_numerator(orientation, &graph, &numerator)?
        );
        assert!(
            specialized
                .replace(external.to_pattern())
                .with(Atom::Zero)
                .is_zero()
        );
        Ok(())
    }

    #[test]
    fn one_source_residue_map_is_a_homomorphism_for_all_factorized_factors() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(
            digraph G {
                edge [particle="scalar_1"];
                node [num=1];
                a -> b [id=0];
                a -> b [id=1];
            },
            "scalars"
        )?;
        let production = vec![OrientationExpression {
            data: OrientationData::new(EdgeVec::from_iter([
                Orientation::Default,
                Orientation::Reversed,
            ])),
            loop_energy_map: Vec::new(),
            edge_energy_map: vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::ose(EdgeIndex(1), -1),
            ],
            variants: Vec::new(),
        }]
        .into_iter()
        .collect::<TiVec<OrientationID, _>>();
        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::default();
        let orientation = OrientationProjection::exact(&production, &options, &pattern, false);

        let a = Atom::var(symbol!("direct_3d_test::A"));
        let b = Atom::var(symbol!("direct_3d_test::B"));
        let host = OrientationID(0);
        let index = CutCFFIndex::new_all_none();
        let key = DirectResidueKey::source(
            host,
            vec![LinearEnergyExpr::uniform_scale(2), LinearEnergyExpr::zero()],
        );
        let q0 = GS.emr_mom(EdgeIndex(0), GS.cind(0));

        // Map both factors independently through the same complete residue
        // key. Keeping them factorized must not change that common Q0 -> 2M
        // substitution authority.
        let first = key.map_numerator(orientation, &graph, &(q0.clone() + &a))?;
        let branches = DirectResidueBranches::from_keyed([(
            key,
            [(index, first.clone())].into_iter().collect(),
        )])?;
        let completed = branches.fallible_map(|key, body| {
            let later = key.map_numerator(orientation, &graph, &(q0.clone() + &b))?;
            Ok(body * later)
        })?;
        let sampled_energy = Atom::num(2) * Atom::var(GS.numerator_sampling_scale);
        let later = &sampled_energy + &b;
        assert_eq!(first, &sampled_energy + &a);
        let expected = first * later;

        let explicit = completed.materialize(false)?;
        assert_eq!(explicit.iter().next().unwrap().1, &expected);
        let localized = completed.materialize(true)?;
        let localized_body = localized.iter().next().unwrap().1;
        assert_eq!(host.select(localized_body.as_view()), expected);
        assert_eq!(
            OrientationID(1).select(localized_body.as_view()),
            Atom::Zero
        );
        // Several cut orders share one map, while two maps on the same host
        // remain different branches. Attach an untouched factorized numerator
        // through the public production operation and compare each cut body
        // with the explicit substitution oracle.
        let raised = CutCFFIndex {
            lu_cut_order: Some(1),
            ..index
        };
        let factor = (q0.clone() + &a) * (&q0 + &b);
        let branches = DirectResidueBranches::from_keyed([2, 3].map(|scale| {
            (
                DirectResidueKey::source(
                    host,
                    vec![
                        LinearEnergyExpr::uniform_scale(scale),
                        LinearEnergyExpr::zero(),
                    ],
                ),
                [(index, Atom::num(5)), (raised, Atom::num(7))]
                    .into_iter()
                    .collect(),
            )
        }))?;
        let mapped = branches.multiply_key_mapped(
            orientation,
            &graph,
            &factor,
            DirectResidueBranches::numerator_scope(),
        )?;
        let expected_sum = [2, 3].into_iter().fold(Atom::Zero, |sum, scale| {
            let energy = Atom::num(scale) * Atom::var(GS.numerator_sampling_scale);
            sum + (&energy + &a) * (&energy + &b)
        });
        let expected = [
            (index, Atom::num(5) * &expected_sum),
            (raised, Atom::num(7) * &expected_sum),
        ]
        .into_iter()
        .collect::<Integrands>();
        let localized = mapped.materialize(true)?.resolved()?;
        for actual in [
            mapped.materialize(false)?.resolved()?,
            localized.map(|atom| host.select(atom.as_view())),
        ] {
            let difference = actual.checked_zip(&expected, |_, actual, expected| {
                Ok((actual - expected).expand())
            })?;
            assert!(
                difference.iter().all(|(_, atom)| atom.is_zero()),
                "every cut must preserve the complete mapped numerator value: {difference:?}"
            );
        }
        assert_eq!(
            localized.map(|atom| OrientationID(1).select(atom.as_view())),
            expected.map(|_| Atom::Zero)
        );
        Ok(())
    }

    #[test]
    fn selected_cut_support_survives_root_to_first_sector_identity() -> Result<()> {
        let selected = CutCFFIndex {
            left_threshold_order: None,
            right_threshold_order: None,
            lu_cut_order: Some(1),
        };
        let body = Atom::var(symbol!("direct_3d_test::selected_cut_body"));
        let root = DirectResidueBranches::production(
            OrientationID(0),
            [(selected, body.clone())].into_iter().collect(),
        )?;

        let identity = root.identity_integrands();
        assert_eq!(
            identity.iter().collect::<Vec<_>>(),
            vec![(&selected, &Atom::one())]
        );
        let first_sector = root.zip_mul_unmapped(&identity)?;
        assert_eq!(first_sector.factorized_sum(), body);
        Ok(())
    }

    #[test]
    fn captured_gl1_numerator_shares_one_body_across_all_residue_rows() -> Result<()> {
        use crate::cff::orientations::GraphOrientation;

        test_initialise()?;
        // Replay the historical source140 Taylor-entry tensor without distributing
        // its vertex sums. This tests mapping and selector support, not a fresh
        // physical generation or the subsequent Taylor operation.
        let fixture: serde_json::Value = serde_json::from_str(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/resources/uv_parametric_numerator/gl1-source140-pre-taylor.json"
        )))?;
        let graph: Graph = fixture["graph"]
            .as_str()
            .unwrap()
            .into_graph(&crate::utils::load_generic_model("sm"))?;
        let mut internal_edges = graph
            .underlying
            .iter_edges()
            .filter_map(|(pair, edge, _)| pair.is_paired().then_some(usize::from(edge)))
            .collect::<Vec<_>>();
        internal_edges.sort_unstable();
        assert_eq!(internal_edges, vec![2, 3, 4, 5, 6]);
        let numerator = symbolica::parse!(fixture["numerator"].as_str().unwrap());
        let rows = fixture["rows"]
            .as_array()
            .unwrap()
            .iter()
            .map(|row| {
                Ok((
                    OrientationID(row["selector"].as_u64().unwrap() as usize),
                    row["signs"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .map(|sign| sign.as_i64().unwrap())
                        .collect::<Vec<_>>(),
                    symbolica::parse!(row["weight"].as_str().unwrap()),
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        assert_eq!(rows.len(), 18);
        for (index, (id, signs, _)) in rows.iter().enumerate() {
            assert_eq!(*id, OrientationID(index));
            assert_eq!(signs.len(), 5);
            assert!(signs.iter().all(|sign| [-1, 1].contains(sign)));
        }
        let orientation_for_signs = |signs: &[i64]| {
            OrientationData::new(EdgeVec::from_iter(
                [Orientation::Undirected; 2]
                    .into_iter()
                    .chain(signs.iter().map(|sign| {
                        if *sign == 1 {
                            Orientation::Default
                        } else {
                            Orientation::Reversed
                        }
                    })),
            ))
        };
        let production =
            rows.iter()
                .map(|(_, signs, _)| OrientationExpression {
                    data: orientation_for_signs(signs),
                    loop_energy_map: Vec::new(),
                    edge_energy_map: [LinearEnergyExpr::zero(), LinearEnergyExpr::zero()]
                        .into_iter()
                        .chain(
                            signs.iter().enumerate().map(|(edge, sign)| {
                                LinearEnergyExpr::ose(EdgeIndex(edge + 2), *sign)
                            }),
                        )
                        .collect(),
                    variants: Vec::new(),
                })
                .collect::<TiVec<OrientationID, _>>();
        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::default();
        let orientation = OrientationProjection::exact(&production, &options, &pattern, false);
        let cut = CutCFFIndex::new_all_none();
        let source = DirectResidueBranches::from_keyed(rows.iter().map(|(id, _, weight)| {
            (
                DirectResidueKey::production(*id),
                Integrands::from_iter([(cut, weight.clone())]),
            )
        }))?;
        let shared = source.multiply_key_mapped(
            orientation,
            &graph,
            &numerator,
            DirectResidueBranches::numerator_scope(),
        )?;
        let definition = &shared.0[0].1.numerators()[0];
        assert_eq!(definition.args.len(), 5);
        assert!(matches!(definition.rhs.as_view(), AtomView::Mul(_)));
        let mut expected = Vec::new();
        for ((key, integrands), (id, _, weight)) in shared.iter_keys().zip(&rows) {
            assert_eq!(key.selector_host, *id);
            assert_eq!(integrands.numerators().len(), 1);
            assert!(Arc::ptr_eq(definition, &integrands.numerators()[0]));
            let body = key.map_numerator(orientation, &graph, &numerator)? * weight;
            assert_eq!(integrands.resolved()?.iter().next().unwrap().1, &body);
            expected.push(body);
        }
        let localized = shared.materialize(true)?;
        assert_eq!(localized.numerators().len(), 1);
        let localized = localized.resolved()?;
        let root = localized.iter().next().unwrap().1;
        // Evaluate the existing physical selector owner on every captured sign
        // row. The expected selector IDs come from the historical truth table,
        // independently of the newly assembled production orientation list.
        let physical_root = root.replace_multiple(
            production
                .iter_enumerated()
                .map(|(id, row)| Replacement::new(id.atom(), row.data.orientation_delta())),
        );
        let truth_table = fixture["truth_table"].as_array().unwrap();
        assert_eq!(truth_table.len(), 32);
        let mut visited = std::collections::BTreeSet::new();
        let mut active = 0;
        for case in truth_table {
            let signs = case["signs"]
                .as_array()
                .unwrap()
                .iter()
                .map(|sign| sign.as_i64().unwrap())
                .collect::<Vec<_>>();
            assert_eq!(signs.len(), 5);
            assert!(signs.iter().all(|sign| [-1, 1].contains(sign)));
            assert!(
                visited.insert(signs.clone()),
                "truth-table rows must be distinct"
            );
            let selectors = case["selectors"].as_array().unwrap();
            assert!(selectors.len() <= 1);
            active += selectors.len();
            let oracle = Atom::add_many(
                selectors
                    .iter()
                    .map(|id| &expected[id.as_u64().unwrap() as usize]),
            );
            assert_eq!(orientation_for_signs(&signs).select(&physical_root), oracle);
        }
        assert_eq!(active, 18);
        assert_eq!(
            shared
                .materialize(false)?
                .resolved()?
                .iter()
                .next()
                .unwrap()
                .1,
            &Atom::add_many(expected)
        );
        Ok(())
    }
}
