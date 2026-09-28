use spenso::structure::partial::{PartialStructure, PartialStructureExt};
use std::{
    collections::{HashMap, HashSet},
    sync::Arc,
};
use symbolica::atom::{AliasedAtom, Atom, AtomCore, AtomView};

use super::{SymbolicTensor, aliases::AliasInterfaces, inference::TensorInferenceError};
use crate::shorthands::schoonschip::{
    ContractionStatus, DotNormalizer, SchoonschipSettings, SchoonschipWithSettings,
    SimplificationCandidates, SlotContraction,
};

/// Contraction scheduling and whether rank-one tensors may form compact products.
/// Metric-only contraction retains explicit vector endpoints, including Vakint inputs.
#[derive(Clone, Copy, Debug)]
pub struct ContractionSettings<'a> {
    order: Option<&'a [usize]>,
    rank_one: bool,
}

impl Default for ContractionSettings<'_> {
    fn default() -> Self {
        Self {
            order: None,
            rank_one: true,
        }
    }
}

impl<'a> ContractionSettings<'a> {
    /// Visit normalized top-level factors in this explicit permutation.
    pub fn with_order(mut self, order: &'a [usize]) -> Self {
        self.order = Some(order);
        self
    }

    /// Contract metrics while retaining explicit rank-one tensor endpoints.
    pub fn without_rank_one_tensors(mut self) -> Self {
        self.rank_one = false;
        self
    }
}

impl SymbolicTensor<PartialStructure> {
    fn contraction_is_noop(&self) -> bool {
        let observed = SimplificationCandidates::scan(self.expression.as_view(), [], || true);
        observed.complete
            && observed.intrinsic
            && observed.traversal_complete
            && !observed.repeated_indices
            && !observed.brackets
            // A no-op leaf needs its exact validated interface, not the
            // stronger self-dual proof required for arbitrary algebraic rewrites.
            && (self.rewrites_preserve_interface()
                || (matches!(self.expression.as_view(), AtomView::Fun(_))
                    && self.structure.open_positions().is_empty()
                    && !self.expression.as_view().needs_normalization()
                    && self.established_interface_is_valid()))
            // The observer conservatively marks explicit vector leaves as dot
            // candidates. Only the existing normalizer can discharge that flag;
            // intrinsic observation above excludes user callbacks from this check.
            && (!observed.dots || DotNormalizer::run(self.expression.as_view()) == self.expression)
    }

    /// Contract the metric/vector graph while keeping generated sums in aliases.
    /// An explicit order permutes the top-level factors after dot normalization.
    /// Materializing the result is a separate operation on the returned tensor.
    pub fn contract(
        &self,
        settings: ContractionSettings<'_>,
    ) -> Result<SymbolicTensor<AliasInterfaces, AliasedAtom>, TensorInferenceError> {
        let contracted = self.contract_parts(settings)?;
        let aliases = contracted
            .aliases
            .into_iter()
            .map(|(handle, body)| {
                let scalar = PartialStructure::from_logical_slots([]);
                Ok((
                    Self::checked_parts(handle, scalar.clone())?,
                    Self::checked_parts(body, scalar)?,
                ))
            })
            .collect::<Result<Vec<_>, TensorInferenceError>>()?;
        let mut result = contracted.root.with_aliases(aliases)?;
        // The completion bit certifies the full metric/vector policy only.
        // A metric-only result can retain explicit vector pairs for a later call.
        result.proofs.contracted =
            settings.rank_one && contracted.status == ContractionStatus::Complete;
        Ok(result)
    }

    pub(crate) fn contract_parts(
        &self,
        settings: ContractionSettings<'_>,
    ) -> Result<crate::shorthands::schoonschip::FactorizedContraction<Self>, TensorInferenceError>
    {
        let order = settings.order;
        let source = DotNormalizer::run(self.expression.as_view());
        if let Some(order) = order {
            let count = match source.as_view() {
                AtomView::Mul(product) => product.iter().len(),
                _ => 1,
            };
            let mut positions = order.to_vec();
            positions.sort_unstable();
            if positions != (0..count).collect::<Vec<_>>() {
                return Err(TensorInferenceError::Invalid(
                    "contraction order must list every normalized top-level factor exactly once"
                        .into(),
                ));
            }
        }
        if order.is_none() && source == self.expression && self.contraction_is_noop() {
            // Complete intrinsic observation and the carried interface proof
            // exclude callbacks, opaque payloads and AUTO. The alias owner
            // still visits every original reachable definition.
            return Ok(crate::shorthands::schoonschip::FactorizedContraction {
                root: self.clone(),
                aliases: Vec::new(),
                status: ContractionStatus::Complete,
                literal_relabellings: Vec::new(),
            });
        }
        if let Some(contracted) =
            SlotContraction::new().contract_factorized(source.as_view(), order, settings.rank_one)
        {
            return Ok(crate::shorthands::schoonschip::FactorizedContraction {
                root: self.with_rewritten_expression(contracted.root)?,
                aliases: contracted.aliases,
                status: contracted.status,
                literal_relabellings: contracted.literal_relabellings,
            });
        }
        // Callback-sensitive leaves retain the existing ordered rewrite and
        // checked finisher. A changed result can expose intrinsic work (for
        // example alpha-equivalent dummy labels after metadata cleanup); retry
        // the same graph owner only on that changed input.
        // Preserve the established network sum schedule: each already-present
        // term is rewritten once, then collected. A nested whole-sum discovery
        // probe would normalize its first callback-sensitive term twice.
        // The typed interface includes ports inside chain/trace words. Their
        // builtin projectors can decline intrinsic graph admission; the ordered
        // reducer must still apply substitutions to those established ports.
        let rewrite = if settings.rank_one {
            SchoonschipSettings::default().with_chain_like_functions()
        } else {
            SchoonschipSettings::default().without_rank1_tensors()
        };
        let mut literal_relabellings = Vec::new();
        let rewritten = match source.as_view() {
            AtomView::Add(sum) => Atom::add_many(sum.iter().map(|term| {
                SchoonschipWithSettings { settings: &rewrite }.run(term, &mut literal_relabellings)
            })),
            _ => SchoonschipWithSettings { settings: &rewrite }
                .run(source.as_view(), &mut literal_relabellings),
        };
        if rewritten != source
            && let Some(contracted) = SlotContraction::new().contract_factorized(
                rewritten.as_view(),
                None,
                settings.rank_one,
            )
        {
            literal_relabellings.extend(contracted.literal_relabellings);
            return Ok(crate::shorthands::schoonschip::FactorizedContraction {
                root: self.with_rewritten_expression(contracted.root)?,
                aliases: contracted.aliases,
                status: contracted.status,
                literal_relabellings,
            });
        }
        let root = self.with_rewritten_expression(rewritten)?;
        // The ordered fallback supports spaces declined by the coefficient
        // collector (including dual representations). Certify only its checked,
        // changed result when the same no-work observation proves it finished.
        // An incomplete collector frontier returned above retains its budget.
        let complete = root.expression != self.expression && root.contraction_is_noop();
        Ok(crate::shorthands::schoonschip::FactorizedContraction {
            root,
            aliases: Vec::new(),
            status: if complete {
                ContractionStatus::Complete
            } else {
                ContractionStatus::Deferred
            },
            literal_relabellings,
        })
    }
}

impl SymbolicTensor<AliasInterfaces, AliasedAtom> {
    /// Contract the root and its reachable definition templates without resolving
    /// their DAG. A completed default result reuses the same sealed allocation.
    /// Disconnected definitions are preserved without invoking their callbacks.
    pub fn contract(
        self: &Arc<Self>,
        settings: ContractionSettings<'_>,
    ) -> Result<Arc<Self>, TensorInferenceError> {
        if let Some(completed) = self.completed_contraction(settings.order) {
            return Ok(completed);
        }
        let mut definitions = self
            .aliases()?
            .into_iter()
            .map(|pair| (pair.0.expression.clone(), pair))
            .collect::<HashMap<_, _>>();
        // Literal relabellings always derive from this pass's uncontracted
        // templates. Keep completed bodies separate until every reachable use
        // has been registered, so traversal order cannot mix rewrite stages.
        let mut contracted_bodies = HashMap::new();
        let root = self.root().contract_parts(settings)?;
        let mut status = root.status;
        // Generated definitions already belong to the frontier just processed.
        // Re-entering them can endlessly wrap a finished weight in fresh aliases,
        // or bypass an incomplete frontier's explicit work budget. Their existing
        // references remain reachable below, but only earlier templates and new
        // literal relabellings require their own contraction in this call.
        let mut processed = root
            .aliases
            .iter()
            .map(|(handle, _)| handle.clone())
            .collect::<HashSet<_>>();
        Self::retain_contracted_definitions(&root, &mut definitions)?;
        let root = root.root;
        loop {
            // Alias lookup is literal and stops at a registered use. Template
            // traversal never opens tensor metadata or guesses a handle family.
            let mut reachable = HashSet::new();
            let mut pending = vec![root.expression.clone()];
            while let Some(expression) = pending.pop() {
                expression.visitor(&mut |node| {
                    if let Some((handle, body)) = definitions.get(node.get_data()) {
                        if reachable.insert(handle.expression.clone()) {
                            let body = contracted_bodies.get(&handle.expression).unwrap_or(body);
                            pending.push(body.expression.clone());
                        }
                        return false;
                    }
                    true
                });
            }
            let mut next = reachable
                .into_iter()
                .filter(|handle| !processed.contains(handle))
                .collect::<Vec<_>>();
            next.sort_unstable_by(|a, b| AtomView::cmp(&a.as_view(), &b.as_view()));
            if next.is_empty() {
                break;
            }
            for handle in next {
                let body = definitions[&handle].1.clone();
                let contracted = body.contract_parts(ContractionSettings {
                    order: None,
                    ..settings
                })?;
                status = status.max(contracted.status);
                processed.extend(contracted.aliases.iter().map(|(handle, _)| handle.clone()));
                Self::retain_contracted_definitions(&contracted, &mut definitions)?;
                contracted_bodies.insert(handle.clone(), contracted.root);
                processed.insert(handle);
            }
        }
        for (handle, body) in contracted_bodies {
            definitions.get_mut(&handle).unwrap().1 = body;
        }
        let mut result = Arc::new(root.with_aliases(definitions.into_values())?);
        if status != ContractionStatus::Capped {
            // Two open definitions can carry contractions across their literal
            // boundaries. Reuse the selected collector's occurrence-local
            // exposure; independently finished templates do not certify their
            // product. Deferred opaque work may expose its registered bodies,
            // but an exhausted frontier must never be restarted. A deferred
            // result stays incomplete until the shared fixed point checks the
            // newly exposed expression in its next contraction pass.
            let reachable = result.reachable_definitions();
            let registry = result
                .aliases()?
                .into_iter()
                .filter(|(handle, _)| reachable.contains(&handle.expression))
                .collect::<Vec<_>>();
            if registry.iter().any(|(_, body)| !body.is_scalar()) {
                use spenso::{
                    network::{library::symbolic::ETS, tags::SPENSO_TAG},
                    structure::slot::{SlotMatch, SlotMatcher},
                };
                let mut slots = SlotMatcher::default();
                result = result.map_domains(|domain, _| {
                    use spenso::network::parsing::AtomStructureExt;
                    // Every template is already contracted. Only repeated
                    // literal ports can connect their definitions in this
                    // domain; opening disjoint ones needlessly distributes
                    // their factored bodies through the collection tape.
                    if !domain.expression.has_repeated_explicit_indices() {
                        return Ok((domain, Vec::new()));
                    }
                    let mut selected_complete = false;
                    let mapped = domain.collect_with_map(
                        Some(&mut |_, _, _| {
                            selected_complete = true;
                            Ok(())
                        }),
                        &registry,
                        |value| {
                            matches!(value, AtomView::Fun(function)
                                if (function.get_symbol() == ETS.metric
                                    || function.get_symbol().has_tag(&SPENSO_TAG.rank1))
                                && function.iter().any(|argument|
                                    matches!(slots.classify(argument), SlotMatch::Explicit(_))))
                        },
                        |selected, available, _| {
                            let contracted = selected.contract_parts(ContractionSettings {
                                order: None,
                                ..settings
                            })?;
                            status = status.max(contracted.status);
                            let mut definitions = available
                                .iter()
                                .cloned()
                                .map(|pair| (pair.0.expression.clone(), pair))
                                .collect::<HashMap<_, _>>();
                            Self::retain_contracted_definitions(&contracted, &mut definitions)?;
                            Ok((contracted.root, definitions.into_values().collect()))
                        },
                    )?;
                    if !selected_complete {
                        status = status.max(ContractionStatus::Deferred);
                    }
                    Ok(mapped)
                })?;
            }
        }
        Arc::make_mut(&mut result).proofs.contracted =
            settings.rank_one && status == ContractionStatus::Complete;
        Ok(result)
    }

    pub(crate) fn retain_contracted_definitions(
        contracted: &crate::shorthands::schoonschip::FactorizedContraction<
            SymbolicTensor<PartialStructure>,
        >,
        definitions: &mut HashMap<
            Atom,
            (
                SymbolicTensor<PartialStructure>,
                SymbolicTensor<PartialStructure>,
            ),
        >,
    ) -> Result<(), TensorInferenceError> {
        for (source, target) in &contracted.literal_relabellings {
            Self::register_literal_use(source, target, definitions)?;
        }
        for (handle, body) in &contracted.aliases {
            let scalar = PartialStructure::from_logical_slots([]);
            let pair = (
                SymbolicTensor::checked_parts(handle.clone(), scalar.clone())?,
                SymbolicTensor::checked_parts(body.clone(), scalar)?,
            );
            if let Some(previous) = definitions.get(handle)
                && previous != &pair
            {
                return Err(TensorInferenceError::invalid(
                    "conflicting contraction alias definitions",
                ));
            }
            definitions.insert(handle.clone(), pair);
        }
        Ok(())
    }
}

#[cfg(test)]
mod alias_tests;
#[cfg(test)]
mod factorized_tests;
#[cfg(test)]
mod network_tests;

#[cfg(test)]
pub mod test {
    use super::super::AbstractIndex;
    use super::*;
    use crate::shorthands::schoonschip::Schoonschip;
    use crate::test_support::test_initialize;
    use spenso::{g, mink, p};

    #[test]
    fn schoonschip_simplifies_metrics() {
        test_initialize();
        let expr = g!(mink!(4, mu), mink!(4, nu)) * p!(mink!(4, nu));

        assert_eq!(
            crate::test_support::contracted_atom(expr.as_view()).unwrap(),
            p!(mink!(4, mu))
        );
    }

    #[test]
    fn metric_only_retains_vector_pairs_and_default_finishes_them() {
        test_initialize();
        let mu = mink!(4, 94401);
        let nu = mink!(4, 94402);
        let spectator = symbolica::parse_lit!((metric_policy_x + metric_policy_y) ^ 30);
        let source = &spectator * g!(&mu, &nu) * p!(&nu) * spenso::q!(&mu);
        let value = SymbolicTensor::infer(source.clone()).unwrap();
        let policy = ContractionSettings::default().without_rank_one_tensors();
        let metric_only = Arc::new(value.contract(policy).unwrap());
        let explicit = metric_only.resolved().unwrap();
        assert_eq!(explicit.expression, &spectator * p!(&nu) * spenso::q!(&nu));
        assert!(!metric_only.contraction_complete());
        assert_eq!(explicit.structure, value.structure);
        assert_eq!(
            metric_only.contract(policy).unwrap().resolved().unwrap(),
            explicit
        );

        let full = metric_only.contract(Default::default()).unwrap();
        assert_eq!(
            full.resolved().unwrap().expression,
            spectator * g!(p!(mink!(4)), spenso::q!(mink!(4)))
        );
        assert!(full.contraction_complete());
        assert!(Arc::ptr_eq(
            &full,
            &full.contract(Default::default()).unwrap()
        ));
    }

    #[test]
    fn metric_only_preserves_factored_branches_and_scoped_alias_powers() {
        test_initialize();
        let mu = mink!(4, 94411);
        let nu = mink!(4, 94412);
        let spectator = symbolica::parse_lit!((metric_alias_x + metric_alias_y) ^ 30);
        let body =
            SymbolicTensor::infer((g!(&mu, &nu) * p!(&nu) * spenso::q!(&mu) + Atom::one()).pow(2))
                .unwrap();
        let handle = body.alias_handle().unwrap();
        let root = SymbolicTensor::infer(spectator.clone())
            .unwrap()
            .multiply(&handle)
            .unwrap();
        let value = Arc::new(root.with_aliases([(handle, body)]).unwrap());
        let metric_only = value
            .contract(ContractionSettings::default().without_rank_one_tensors())
            .unwrap();
        assert!(!metric_only.contraction_complete());
        assert_eq!(
            metric_only.resolved().unwrap().expression,
            &spectator * (p!(&nu) * spenso::q!(&nu) + Atom::one()).pow(2)
        );
        let full = metric_only.contract(Default::default()).unwrap();
        assert_eq!(
            full.resolved().unwrap().expression,
            spectator * (g!(p!(mink!(4)), spenso::q!(mink!(4))) + Atom::one()).pow(2)
        );
        assert!(full.contraction_complete());
    }

    #[test]
    fn metric_only_relabels_final_vector_port_and_preserves_metadata() {
        use symbolica::atom::FunctionBuilder;
        test_initialize();
        let mu = mink!(4, 94421);
        let nu = mink!(4, 94422);
        let metadata = symbolica::function!(
            symbolica::symbol!("metric_policy_scalar_metadata"; Scalar),
            mink!(4, 94423)
        );
        let vector = spenso::vector_symbol!("metric_policy_metadata");
        let vector_at = |port: &Atom| {
            FunctionBuilder::new(vector)
                .add_arg(&metadata)
                .add_arg(port)
                .finish()
        };
        let invalid = FunctionBuilder::new(vector)
            .add_arg(mink!(4, 94423))
            .add_arg(&nu)
            .finish();
        assert!(SymbolicTensor::infer(g!(&mu, &nu) * invalid).is_err());
        let value = SymbolicTensor::infer(g!(&mu, &nu) * vector_at(&nu)).unwrap();
        assert!(
            SlotContraction::new()
                .contract_factorized(value.expression.as_view(), None, false)
                .is_some()
        );
        let result = value
            .contract(ContractionSettings::default().without_rank_one_tensors())
            .unwrap()
            .resolved()
            .unwrap();
        assert_eq!(result.expression, vector_at(&mu));
        assert_eq!(result.structure, value.structure);
    }

    #[test]
    fn contraction_returns_checked_aliases_and_preserves_typed_zero() {
        use spenso::structure::{
            partial::{PartialIndex, PartialStructureExt},
            representation::{LibraryRep, Minkowski, RepName},
            slot::IsAbstractSlot,
        };
        use symbolica::atom::{Atom, AtomCore, FunctionBuilder};

        test_initialize();
        let p = spenso::vector_symbol!("contract_alias_p");
        let q = spenso::vector_symbol!("contract_alias_q");
        let r = spenso::vector_symbol!("contract_alias_r");
        let rep = LibraryRep::from(Minkowski {}).new_rep(4);
        let slot = rep.slot::<AbstractIndex, _>(92101).to_atom();
        let vector = |head| FunctionBuilder::new(head).add_arg(&slot).finish();
        let expression = (vector(p) + vector(q)) * vector(r);
        let source = SymbolicTensor::infer(expression).unwrap();
        let contracted = source.contract(Default::default()).unwrap();
        assert!(!contracted.expression.get_aliases().is_empty());
        assert_eq!(
            contracted.resolved().unwrap().expression.expand(),
            source.expression.expand().schoonschip().expand()
        );
        let zero = SymbolicTensor::checked_parts(
            Atom::Zero,
            PartialStructure::from_logical_slots([
                rep.slot(PartialIndex::Explicit(AbstractIndex::from(92103)))
            ]),
        )
        .unwrap();
        assert_eq!(zero.contract(Default::default()).unwrap().root(), zero);
    }

    #[test]
    fn generated_contraction_weights_are_not_reentered_as_new_templates() {
        test_initialize();
        let left = p!(mink!(4, 92801));
        let right = spenso::q!(mink!(4, 92801));
        let spectator = p!(mink!(4, 92801));
        let body = SymbolicTensor::infer((&left + &right) * &spectator).unwrap();
        let handle = body.alias_handle().unwrap();
        let source = Arc::new(handle.clone().with_aliases([(handle, body)]).unwrap());
        let result = source.contract(Default::default()).unwrap();
        assert!(
            result.expression.get_aliases().len() <= 8,
            "one completed frontier must retain a finite definition registry"
        );
        let expected = g!(p!(mink!(4)), p!(mink!(4))) + g!(spenso::q!(mink!(4)), p!(mink!(4)));
        assert_eq!(
            result.resolved().unwrap().expression.expand(),
            expected.expand()
        );
        assert!(Arc::ptr_eq(
            &result,
            &result.contract(Default::default()).unwrap()
        ));
    }

    #[test]
    fn completed_noop_contraction_does_not_wrap_existing_coefficients() {
        test_initialize();
        let a = mink!(4, 92803);
        let b = mink!(4, 92804);
        let weight = SymbolicTensor::infer(
            Atom::var(symbolica::symbol!("noop_contraction_weight"))
                * g!(p!(mink!(4)), p!(mink!(4))),
        )
        .unwrap();
        let handle = weight.alias_handle().unwrap();
        let metric = SymbolicTensor::infer(g!(&a, &b)).unwrap();
        let root = metric.multiply(&handle).unwrap();
        let value = Arc::new(root.with_aliases([(handle, weight)]).unwrap());
        assert!(!value.contraction_complete());
        let contracted = value.contract(Default::default()).unwrap();
        assert!(contracted.contraction_complete());
        assert_eq!(contracted.root(), value.root());
        assert_eq!(contracted.aliases().unwrap(), value.aliases().unwrap());

        // Domain rewrites legitimately clear completion. Reconstructing the
        // same checked registry must still be stable without that shortcut.
        let unattested = Arc::new(
            contracted
                .root()
                .with_aliases(contracted.aliases().unwrap())
                .unwrap(),
        );
        assert!(!unattested.contraction_complete());
        let repeated = unattested.contract(Default::default()).unwrap();
        assert_eq!(repeated.root(), contracted.root());
        assert_eq!(repeated.aliases().unwrap(), contracted.aliases().unwrap());
        assert!(repeated.contraction_complete());
        assert!(Arc::ptr_eq(
            &repeated,
            &repeated.contract(Default::default()).unwrap()
        ));
    }

    #[test]
    fn contraction_noop_observation_keeps_free_metric_polynomials_factored() {
        use spenso::structure::{representation::RepName, slot::IsAbstractSlot};
        test_initialize();
        let [a, b, c, d] = [92811, 92812, 92813, 92814].map(|index| {
            spenso::structure::representation::LibraryRep::from(
                spenso::structure::representation::Minkowski {},
            )
            .new_rep(4)
            .slot::<AbstractIndex, _>(index)
            .to_atom()
        });
        let spectator = symbolica::parse_lit!((noop_x + noop_y) ^ 8);
        let expression = &spectator * (g!(&a, &b) * g!(&c, &d) + g!(&a, &c) * g!(&b, &d));
        let source = SymbolicTensor::infer(expression.clone()).unwrap();
        let result = source.contract(Default::default()).unwrap();
        assert!(result.contraction_complete());
        assert_eq!(result.root(), source);
        assert!(result.expression.get_aliases().is_empty());
        assert_eq!(result.resolved().unwrap().expression, expression);
    }

    #[test]
    fn contraction_noop_root_still_reduces_scoped_definition_powers() {
        use spenso::structure::{
            abstract_index::AbstractIndex,
            representation::{LibraryRep, Minkowski, RepName},
            slot::IsAbstractSlot,
        };
        use symbolica::atom::FunctionBuilder;
        test_initialize();
        let rep = LibraryRep::from(Minkowski {}).new_rep(4);
        let index = AbstractIndex::Normal(92821).scoped(symbolica::symbol!("noop_scope"));
        let slot = rep.slot::<AbstractIndex, _>(index).to_atom();
        let p = spenso::vector_symbol!("noop_scoped_p");
        let q = spenso::vector_symbol!("noop_scoped_q");
        let vector = |head, arg: &Atom| FunctionBuilder::new(head).add_arg(arg).finish();
        let x = Atom::var(symbolica::symbol!("noop_scoped_weight"));
        for power in [2, 3] {
            let body = SymbolicTensor::infer((vector(p, &slot) * vector(q, &slot) + &x).pow(power))
                .unwrap();
            let handle = body.alias_handle().unwrap();
            let source = Arc::new(handle.clone().with_aliases([(handle, body)]).unwrap());
            let result = source.contract(Default::default()).unwrap();
            let compact = rep.to_symbolic([]);
            let expected = (g!(vector(p, &compact), vector(q, &compact)) + &x).pow(power);
            assert_eq!(
                result.resolved().unwrap().expression.expand(),
                expected.expand()
            );
            assert!(result.contraction_complete());
            assert!(Arc::ptr_eq(
                &result,
                &result.contract(Default::default()).unwrap()
            ));
        }
    }

    #[test]
    fn contraction_checks_callback_induced_rank_loss() {
        use spenso::network::library::symbolic::ETS;
        use symbolica::atom::{Atom, AtomView, FunctionBuilder};

        test_initialize();
        let a = mink!(4, 92201);
        let b = mink!(4, 92203);
        let target = b.clone();
        let head = spenso::tensor_symbol!(
            "contract_alias_callback",
            norm = move |node, output| {
                if let AtomView::Fun(function) = node
                    && function.iter().any(|argument| argument == target.as_view())
                {
                    **output = Atom::one();
                }
            }
        );
        let metric = FunctionBuilder::new(ETS.metric)
            .add_arg(&a)
            .add_arg(&b)
            .finish();
        let tensor = FunctionBuilder::new(head).add_arg(&a).finish();
        let source = SymbolicTensor::infer(metric * tensor).unwrap();
        assert!(source.contract(Default::default()).is_err());
        assert!(
            source
                .contract(ContractionSettings::default().without_rank_one_tensors())
                .is_err()
        );
    }
}

#[cfg(test)]
mod completion_tests {
    use super::*;
    use spenso::{
        network::library::symbolic::ETS,
        structure::{
            abstract_index::AbstractIndex,
            representation::{LibraryRep, RepName},
            slot::IsAbstractSlot,
        },
    };
    use std::sync::atomic::{AtomicUsize, Ordering};
    use symbolica::atom::FunctionBuilder;

    #[test]
    fn dual_fallback_completion_keeps_callback_refusal() {
        crate::test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let vector = spenso::vector_symbol!(
            "completion_callback_q",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let rep = LibraryRep::new_dual("completion_callback::R")
            .unwrap()
            .new_rep(4);
        let i = rep
            .slot::<AbstractIndex, _>(AbstractIndex::Symbol(
                symbolica::symbol!("completion_callback::i").into(),
            ))
            .to_atom();
        let j = AbstractIndex::Symbol(symbolica::symbol!("completion_callback::j").into());
        let expected = FunctionBuilder::new(vector).add_arg(&i).finish();
        let expression = FunctionBuilder::new(ETS.metric)
            .add_arg(&i)
            .add_arg(rep.dual().slot::<AbstractIndex, _>(j).to_atom())
            .finish()
            * FunctionBuilder::new(vector)
                .add_arg(rep.slot::<AbstractIndex, _>(j).to_atom())
                .finish();
        let source = SymbolicTensor::infer(expression).unwrap();
        calls.store(0, Ordering::Relaxed);
        let result = source.contract(Default::default()).unwrap();
        assert!(!result.contraction_complete());
        assert!(calls.load(Ordering::Relaxed) > 0);
        assert_eq!(result.resolved().unwrap().expression(), &expected);
    }
}
