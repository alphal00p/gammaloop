//! Shared scheduling of tensor algebra passes. Materialization is an explicit
//! operation on the result, never a cleanup step of this scheduler.

use spenso::{
    network::tags::SPENSO_TAG,
    structure::{partial::PartialStructure, representation::LibraryRep},
};
use std::sync::{Arc, LazyLock};
use symbolica::atom::{AliasedAtom, AtomView, Symbol};

use crate::{
    color::ColorSimplifySettings,
    dirac::{AGS, GammaOutput, GammaSimplifySettings},
    epsilon::EpsilonSimplifierPass,
    representations::{Bispinor, ColorAdjoint, ColorFundamental, ColorSextet},
    shorthands::{chain::Chain, schoonschip::SimplificationCandidates},
};

use super::{SymbolicTensor, aliases::AliasInterfaces, inference::TensorInferenceError};

// Share the contractor's syntactic observations across domain passes. The
// representation heads select domains; the remaining heads admit open words.
static DOMAIN_HEADS: LazyLock<[Symbol; 17]> = LazyLock::new(|| {
    [
        LibraryRep::from(Bispinor {}).symbol(),
        LibraryRep::from(ColorFundamental {}).symbol(),
        LibraryRep::from(ColorAdjoint {}).symbol(),
        LibraryRep::from(ColorSextet {}).symbol(),
        *crate::epsilon::EPSILON_SYMBOL,
        SPENSO_TAG.chain,
        SPENSO_TAG.trace,
        SPENSO_TAG.bracket,
        AGS.gamma,
        AGS.gamma0,
        AGS.gamma5,
        AGS.projm,
        AGS.projp,
        AGS.sigma,
        AGS.charge_conjugation,
        AGS.gammaconj,
        AGS.gammaadj,
    ]
});

/// Select tensor identities without changing their dimension conventions.
#[derive(Clone, Copy, Debug)]
pub struct SimplifySettings {
    pub metrics: bool,
    pub gamma: Option<GammaSimplifySettings>,
    pub color: Option<ColorSimplifySettings>,
    pub epsilon: bool,
    pub max_passes: usize,
}

impl Default for SimplifySettings {
    fn default() -> Self {
        Self {
            metrics: true,
            gamma: None,
            color: None,
            epsilon: false,
            max_passes: 16,
        }
    }
}

impl SimplifySettings {
    pub fn hep() -> Self {
        Self {
            gamma: Some(GammaSimplifySettings::default()),
            color: Some(ColorSimplifySettings::default()),
            epsilon: true,
            ..Self::default()
        }
    }

    pub fn validate(&self) -> Result<(), TensorInferenceError> {
        if self.max_passes == 0 {
            return Err(TensorInferenceError::invalid("max_passes must be positive"));
        }
        Ok(())
    }
}

impl SymbolicTensor<PartialStructure> {
    /// Complete selected identities, retaining typed literal definitions.
    pub fn simplify(
        &self,
        settings: &SimplifySettings,
    ) -> Result<Arc<SymbolicTensor<AliasInterfaces, AliasedAtom>>, TensorInferenceError> {
        Arc::new(self.clone().with_aliases([])?).simplify(settings)
    }

    /// Simplify Dirac identities without materializing trace definitions.
    pub fn simplify_gamma(
        &self,
        settings: GammaSimplifySettings,
    ) -> Result<Arc<SymbolicTensor<AliasInterfaces, AliasedAtom>>, TensorInferenceError> {
        self.simplify(&SimplifySettings {
            metrics: settings.output == GammaOutput::Reduced,
            gamma: Some(settings),
            epsilon: settings.output == GammaOutput::Reduced,
            ..SimplifySettings::default()
        })
    }

    pub fn simplify_color(
        &self,
        settings: ColorSimplifySettings,
    ) -> Result<Arc<SymbolicTensor<AliasInterfaces, AliasedAtom>>, TensorInferenceError> {
        self.simplify(&SimplifySettings {
            metrics: settings.simplify_non_color,
            color: Some(settings),
            ..SimplifySettings::default()
        })
    }

    pub fn simplify_epsilon(
        &self,
    ) -> Result<Arc<SymbolicTensor<AliasInterfaces, AliasedAtom>>, TensorInferenceError> {
        self.simplify(&SimplifySettings {
            epsilon: true,
            ..SimplifySettings::default()
        })
    }
}

impl SymbolicTensor<AliasInterfaces, AliasedAtom> {
    /// Run one shared fixed point over the root and reachable definitions.
    /// Domain passes keep handles opaque; the contractor registers any literal
    /// port changes before a later pass can observe them.
    pub fn simplify(
        self: &Arc<Self>,
        settings: &SimplifySettings,
    ) -> Result<Arc<Self>, TensorInferenceError> {
        settings.validate()?;
        let mut current = Arc::clone(self);
        for _ in 0..settings.max_passes {
            let before = Arc::clone(&current);
            if settings.metrics {
                current = current.contract(Default::default())?;
            }
            if settings.metrics
                || settings.gamma.is_some()
                || settings.color.is_some()
                || settings.epsilon
            {
                let registry = current.aliases()?;
                // Epsilon has no identifying representation port: its body can
                // be hidden entirely behind a Minkowski tensor alias.
                let registry_epsilon = settings.epsilon
                    && registry.iter().any(|(_, body)| {
                        let observed = SimplificationCandidates::scan(
                            body.expression.as_view(),
                            [*crate::epsilon::EPSILON_SYMBOL],
                            || true,
                        );
                        observed.symbols[0]
                    });
                current = current.map_domains(|mut domain, _| {
                    let mut definitions = Vec::new();
                    let mut observed = SimplificationCandidates::scan(
                        domain.expression.as_view(),
                        *DOMAIN_HEADS,
                        || true,
                    );
                    if settings.metrics && (!observed.complete || observed.symbols[5]) {
                        let expression = domain.expression.normalize_chains();
                        if expression != domain.expression {
                            domain = domain.with_rewritten_expression(expression)?;
                            observed = SimplificationCandidates::scan(
                                domain.expression.as_view(),
                                *DOMAIN_HEADS,
                                || true,
                            );
                        }
                    }
                    if let Some(gamma) = settings.gamma
                        && (!observed.complete
                            || observed.symbols[0]
                            || observed.symbols[5..].iter().any(|&seen| seen))
                    {
                        let previous = domain.expression.clone();
                        (domain, definitions) = domain.simplify_gamma_parts(&registry, gamma)?;
                        if domain.expression != previous {
                            observed = SimplificationCandidates::scan(
                                domain.expression.as_view(),
                                *DOMAIN_HEADS,
                                || true,
                            );
                        }
                    }
                    if let Some(color) = settings.color
                        && (!observed.complete || observed.symbols[1..4].iter().any(|&seen| seen))
                    {
                        let previous = domain.expression.clone();
                        // The next domain can use handles emitted by this root's
                        // earlier stage; their bodies still wait for the next pass.
                        let available = registry
                            .iter()
                            .chain(&definitions)
                            .cloned()
                            .collect::<Vec<_>>();
                        let (rewritten, emitted) =
                            domain.simplify_color_parts(&available, color)?;
                        domain = rewritten;
                        definitions.extend(emitted);
                        if domain.expression != previous {
                            observed = SimplificationCandidates::scan(
                                domain.expression.as_view(),
                                *DOMAIN_HEADS,
                                || true,
                            );
                        }
                    }
                    if settings.epsilon && (observed.symbols[4] || registry_epsilon)
                    {
                        use spenso::{
                            network::library::symbolic::ETS,
                            structure::slot::{SlotMatch, SlotMatcher},
                        };
                        let available = registry.iter().chain(&definitions).cloned().collect::<Vec<_>>();
                        let mut slots = SlotMatcher::default();
                        let (rewritten, emitted) = domain.collect_with_map(
                            None,
                            &available,
                            |value| matches!(value, AtomView::Fun(function)
                                if function.get_symbol() == *crate::epsilon::EPSILON_SYMBOL
                                    || ((function.get_symbol() == ETS.metric
                                        || function.get_symbol().has_tag(&SPENSO_TAG.rank1))
                                        && function.iter().any(|argument|
                                            matches!(slots.classify(argument), SlotMatch::Explicit(_))))),
                            |selected, _, _| {
                                let expression = EpsilonSimplifierPass::step(selected.expression.as_view());
                                Ok((selected.with_rewritten_expression(expression)?, Vec::new()))
                            },
                        )?;
                        domain = rewritten;
                        definitions.extend(emitted);
                    } else if settings.epsilon && !observed.complete {
                        // An opaque label is not evidence of an epsilon domain.
                        // Retain the conservative primitive visit without opening
                        // unrelated metric/vector definitions into a selected tape.
                        let expression = EpsilonSimplifierPass::step(domain.expression.as_view());
                        domain = domain.with_rewritten_expression(expression)?;
                    }
                    Ok((domain, definitions))
                })?;
            }
            if settings.metrics {
                current = current.contract(Default::default())?;
            }
            if Arc::ptr_eq(&current, &before)
                || (current.expression.get_root() == before.expression.get_root()
                    && current.expression.get_aliases() == before.expression.get_aliases()
                    && current.structure.root() == before.structure.root()
                    && current.structure.definitions() == before.structure.definitions()
                    && current.is_metric == before.is_metric
                    && current.is_composite == before.is_composite)
            {
                return Ok(current);
            }
        }
        Err(TensorInferenceError::invalid(format!(
            "simplification did not stabilize within {} passes",
            settings.max_passes
        )))
    }

    pub fn simplify_gamma(
        self: &Arc<Self>,
        settings: GammaSimplifySettings,
    ) -> Result<Arc<Self>, TensorInferenceError> {
        self.simplify(&SimplifySettings {
            metrics: settings.output == GammaOutput::Reduced,
            gamma: Some(settings),
            epsilon: settings.output == GammaOutput::Reduced,
            ..SimplifySettings::default()
        })
    }

    pub fn simplify_color(
        self: &Arc<Self>,
        settings: ColorSimplifySettings,
    ) -> Result<Arc<Self>, TensorInferenceError> {
        self.simplify(&SimplifySettings {
            metrics: settings.simplify_non_color,
            color: Some(settings),
            ..SimplifySettings::default()
        })
    }

    pub fn simplify_epsilon(self: &Arc<Self>) -> Result<Arc<Self>, TensorInferenceError> {
        self.simplify(&SimplifySettings {
            epsilon: true,
            ..SimplifySettings::default()
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::{g, mink, p};
    use symbolica::{
        atom::{Atom, AtomCore},
        symbol,
    };

    #[test]
    fn compact_inner_product_gamma_reduces_both_raw_and_chain_operands() {
        use crate::gamma;
        use spenso::{chain, network::tags::SPENSO_TAG};
        use symbolica::function;
        let reps = crate::test_support::test_initialize();
        let left = reps.bis4.to_symbolic([Atom::num(91301)]);
        let inner = reps.bis4.to_symbolic([Atom::num(91302)]);
        let right = reps.bis4.to_symbolic([Atom::num(91303)]);
        let scalar_dot = function!(
            SPENSO_TAG.dot,
            p!(reps.mink4.to_symbolic([])),
            spenso::q!(reps.mink4.to_symbolic([]))
        );
        for dimension in [
            Atom::num(4),
            Atom::var(symbol!("compact_inner_product_dimension")),
        ] {
            let compact = function!(
                LibraryRep::from(spenso::structure::representation::Minkowski {}).symbol(),
                &dimension
            );
            let direct = function!(
                SPENSO_TAG.dot,
                gamma!(&left, &inner, &compact),
                gamma!(&inner, &right, &compact)
            );
            let chained = function!(
                SPENSO_TAG.dot,
                chain!(&left, &inner, gamma!(&compact)),
                chain!(&inner, &right, gamma!(&compact))
            );
            for expression in [direct, chained] {
                let source =
                    SymbolicTensor::<PartialStructure>::infer(&scalar_dot * expression).unwrap();
                let collected = source
                    .simplify_gamma(GammaSimplifySettings {
                        output: GammaOutput::Chains,
                        ..GammaSimplifySettings::default()
                    })
                    .unwrap();
                assert_eq!(collected.root().structure, source.structure);
                let result = source
                    .simplify_gamma(GammaSimplifySettings::default())
                    .unwrap();
                assert_eq!(result.root().structure, source.structure);
                assert_eq!(
                    result.resolved().unwrap().expression,
                    &dimension * &scalar_dot * g!(&left, &right)
                );
                let handle = source.alias_handle().unwrap();
                let aliased = Arc::new(
                    handle
                        .clone()
                        .with_aliases([(handle, source.clone())])
                        .unwrap(),
                );
                let alias_result = aliased
                    .simplify_gamma(GammaSimplifySettings::default())
                    .unwrap();
                assert_eq!(alias_result.root().structure, source.structure);
                assert_eq!(
                    alias_result.resolved().unwrap().expression,
                    result.resolved().unwrap().expression
                );
                let repeated = collected
                    .simplify_gamma(GammaSimplifySettings::default())
                    .unwrap();
                assert_eq!(
                    repeated.resolved().unwrap().expression,
                    result.resolved().unwrap().expression
                );
            }
        }
    }

    #[test]
    fn shared_pipeline_preserves_factored_spectators_and_typed_zero() {
        crate::test_support::test_initialize();
        let x = Atom::var(symbol!("shared_pipeline_x"));
        let y = Atom::var(symbol!("shared_pipeline_y"));
        let source = SymbolicTensor::<PartialStructure>::infer(
            (x + Atom::one()) * (y + Atom::one()) * p!(mink!(4, mu)),
        )
        .unwrap();
        let settings = SimplifySettings::hep();
        assert_eq!(source.simplify(&settings).unwrap().root(), source);
        let zero = source.with_rewritten_expression(Atom::Zero).unwrap();
        assert_eq!(zero.simplify(&settings).unwrap().root(), zero);
    }

    #[test]
    fn shared_pipeline_enforces_its_fixed_point_budget() {
        crate::test_support::test_initialize();
        let source = SymbolicTensor::<PartialStructure>::infer(
            g!(mink!(4, mu), mink!(4, nu)) * p!(mink!(4, nu)),
        )
        .unwrap();
        let mut settings = SimplifySettings::default();
        let result = source.simplify(&settings).unwrap();
        assert_eq!(result.root().expression, p!(mink!(4, mu)));
        assert_eq!(result.root().structure, source.structure);
        let rerun = result.simplify(&settings).unwrap();
        assert!(Arc::ptr_eq(&result, &rerun));
        assert_eq!(rerun.root(), result.root());
        assert_eq!(
            rerun.expression.get_aliases(),
            result.expression.get_aliases()
        );
        settings.max_passes = 1;
        assert!(source.simplify(&settings).is_err());
        settings.max_passes = 0;
        assert!(source.simplify(&settings).is_err());
    }

    #[test]
    fn shared_pipeline_registers_trace_alias_port_changes() {
        use crate::{gamma, test_support::test_initialize};
        use spenso::trace;

        let reps = test_initialize();
        let slots = (0..6)
            .map(|i| reps.mink_d.to_symbolic([Atom::num(93700 + i)]))
            .collect::<Vec<_>>();
        let trace = trace!(reps.bis4.to_symbolic([]);
            slots[..4].iter().map(|slot| gamma!(slot)));
        let traced = SymbolicTensor::<PartialStructure>::infer(trace.clone())
            .unwrap()
            .simplify_gamma(GammaSimplifySettings::default())
            .unwrap();
        let metric = g!(&slots[0], &slots[4]);
        let input =
            SymbolicTensor::<PartialStructure>::infer(&metric * &traced.root().expression).unwrap();
        let source = Arc::new(
            input
                .clone()
                .with_aliases(traced.aliases().unwrap())
                .unwrap(),
        );
        let result = source
            .simplify_gamma(GammaSimplifySettings::default())
            .unwrap();
        let reference_input = SymbolicTensor::<PartialStructure>::infer(metric * trace).unwrap();
        let reference = reference_input
            .simplify_gamma(GammaSimplifySettings::default())
            .unwrap();
        assert_eq!(result.root().structure, input.structure);
        assert_eq!(reference.root().structure, reference_input.structure);
        // The two syntactic inputs establish different logical port orders.
        // Each result retains its own order; compare their algebra by labels.
        assert_eq!(
            result.expanded().unwrap().expression,
            reference.expanded().unwrap().expression
        );
        assert!(result.aliases().unwrap().len() > traced.aliases().unwrap().len());
    }

    #[test]
    fn shared_pipeline_checks_callback_rank_loss_inside_a_relabelled_alias() {
        use symbolica::atom::{AtomView, FunctionBuilder};

        crate::test_support::test_initialize();
        let a = mink!(4, 93801);
        let b = mink!(4, 93803);
        let target = b.clone();
        let head = spenso::tensor_symbol!(
            "shared_pipeline_alias_callback",
            norm = move |node, out| {
                if let AtomView::Fun(function) = node
                    && function.iter().any(|argument| argument == target.as_view())
                {
                    **out = Atom::one();
                }
            }
        );
        let body = SymbolicTensor::<PartialStructure>::infer(
            FunctionBuilder::new(head).add_arg(&a).finish(),
        )
        .unwrap();
        let handle = body.alias_handle().unwrap();
        let root =
            SymbolicTensor::<PartialStructure>::infer(g!(&a, &b) * &handle.expression).unwrap();
        let source = Arc::new(root.with_aliases([(handle, body)]).unwrap());
        assert!(source.simplify(&SimplifySettings::hep()).is_err());
    }
    #[test]
    fn shared_gamma_collection_combines_an_owned_alias_with_an_adjacent_word() {
        use crate::gamma;
        use spenso::{chain, q};

        let reps = crate::test_support::test_initialize();
        let left = reps.bis4.to_symbolic([Atom::num(93901)]);
        let right = reps.bis4.to_symbolic([Atom::num(93902)]);
        let momentum = reps.mink4.to_symbolic([]);
        let p = p!(&momentum);
        let q = q!(&momentum);
        let r = p!(1, &momentum);
        let body = SymbolicTensor::<PartialStructure>::infer(
            chain!(&left, &right, gamma!(&p)) + chain!(&left, &right, gamma!(&q)),
        )
        .unwrap();
        let handle = body.alias_handle().unwrap();
        let spectator = Atom::var(symbol!("shared_gamma_alias_spectator")) + Atom::one();
        let root = SymbolicTensor::<PartialStructure>::infer(
            &spectator * &handle.expression * chain!(&right, &left, gamma!(&r)),
        )
        .unwrap();
        let source = Arc::new(root.clone().with_aliases([(handle, body)]).unwrap());
        let collected = source
            .simplify_gamma(GammaSimplifySettings {
                output: GammaOutput::Chains,
                ..GammaSimplifySettings::default()
            })
            .unwrap();
        assert_eq!(collected.root().structure, root.structure);
        let result = collected
            .simplify_gamma(GammaSimplifySettings::default())
            .unwrap();
        let expected = Atom::num(4) * spectator * (g!(&p, &r) + g!(&q, &r));
        assert_eq!(result.expanded().unwrap().expression, expected.expand());
    }

    #[test]
    fn shared_color_collection_combines_an_owned_alias_with_an_adjacent_generator() {
        use crate::color_t;

        crate::test_support::test_initialize();
        let a_ik = color_t!([8, 94001], [3, 94003], [3, 94005]);
        let b_kj = color_t!([8, 94002], [3, 94005], [3, 94004]);
        let b_ik = color_t!([8, 94002], [3, 94003], [3, 94005]);
        let a_kj = color_t!([8, 94001], [3, 94005], [3, 94004]);
        let body = SymbolicTensor::<PartialStructure>::infer(a_ik * b_kj + b_ik * a_kj).unwrap();
        let handle = body.alias_handle().unwrap();
        let root = SymbolicTensor::<PartialStructure>::infer(
            &handle.expression * color_t!([8, 94002], [3, 94004], [3, 94006]),
        )
        .unwrap();
        let source = Arc::new(root.clone().with_aliases([(handle, body)]).unwrap());
        let result = source
            .simplify_color(ColorSimplifySettings {
                substitute_cof_dimension_invariants: true,
                ..ColorSimplifySettings::default()
            })
            .unwrap();
        assert_eq!(result.root().structure, root.structure);
        // Sum_b {T^a,T^b} T^b = (2 C_F - C_A/2) T^a = 7 T^a / 6.
        let expected = Atom::num(7) * color_t!([8, 94001], [3, 94003], [3, 94006]) / Atom::num(6);
        // The color output may retain its one-generator chain. Compare in the
        // existing explicit tensor representation, preserving the exact coefficient.
        assert_eq!(
            result.expanded().unwrap().expression.undo_single_length(),
            expected
        );
    }
    #[test]
    fn massive_chiral_current_product_registers_all_literal_aliases() {
        use spenso::structure::partial::PartialStructureExt;

        crate::test_support::test_initialize();
        let _ = (
            spenso::tensor_symbol!("python::p1"),
            spenso::tensor_symbol!("python::p2"),
            spenso::tensor_symbol!("python::k1"),
            spenso::tensor_symbol!("python::k2"),
        );
        // Exact factorized spin-sum inputs from the installed FeynCalc tree test.
        let incoming = Atom::parse(
            r#"1/4*(UFO::Me*spenso::g(spenso::bis(4,python::i2),spenso::bis(4,python::i3))
    +spenso::gamma(spenso::bis(4,python::i2),spenso::bis(4,python::i3),
        python::p1(spenso::mink(4))
    )
)*(-UFO::Me*spenso::g(spenso::bis(4,python::i0),spenso::bis(4,python::i1))
        +spenso::gamma(spenso::bis(4,python::i0),spenso::bis(4,python::i1),
            python::p2(spenso::mink(4))
        )
    )*spenso::chain(spenso::bis(4,python::i1),spenso::bis(4,python::i2),
        spenso::gamma(spenso::in,spenso::out,spenso::mink(4,python::mu)),
        spenso::projp(spenso::in,spenso::out)
    )*spenso::chain(spenso::bis(4,python::i3),spenso::bis(4,python::i0),
        spenso::gamma(spenso::in,spenso::out,spenso::mink(4,python::nu)),
        spenso::projp(spenso::in,spenso::out)
    )"#,
            "python",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap();
        let outgoing = Atom::parse(
            r#"(UFO::MM*spenso::g(spenso::bis(4,python::j0),spenso::bis(4,python::j1))
    +spenso::gamma(spenso::bis(4,python::j0),spenso::bis(4,python::j1),
        python::k1(spenso::mink(4))
    )
)*(-UFO::MM*spenso::g(spenso::bis(4,python::j2),spenso::bis(4,python::j3))
        +spenso::gamma(spenso::bis(4,python::j2),spenso::bis(4,python::j3),
            python::k2(spenso::mink(4))
        )
    )*spenso::chain(spenso::bis(4,python::j1),spenso::bis(4,python::j2),
        spenso::gamma(spenso::in,spenso::out,spenso::mink(4,python::mu)),
        spenso::projp(spenso::in,spenso::out)
    )*spenso::chain(spenso::bis(4,python::j3),spenso::bis(4,python::j0),
        spenso::gamma(spenso::in,spenso::out,spenso::mink(4,python::nu)),
        spenso::projp(spenso::in,spenso::out)
    )"#,
            "python",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap();
        let product = &incoming * &outgoing;
        let masses = [symbol!("UFO::Me"), symbol!("UFO::MM")];
        let massless = product.replace_map(|value, _, output| {
            if matches!(value, symbolica::atom::AtomView::Var(var) if masses.contains(&var.get_symbol())) {
                **output = Atom::Zero;
            }
        });
        for (name, rank, input) in [
            ("incoming", 2, incoming),
            ("outgoing", 2, outgoing),
            ("product", 0, product),
            ("massless", 0, massless),
        ] {
            let source = SymbolicTensor::<PartialStructure>::infer(input).unwrap();
            assert_eq!(source.structure.logical_slots().len(), rank, "{name}");
            let result = source
                .simplify_gamma(GammaSimplifySettings::default())
                .unwrap_or_else(|error| panic!("{name}: {error}"));
            assert_eq!(result.root().structure, source.structure, "{name}");
            assert!(
                result.resolved().is_ok(),
                "{name}: all literal aliases must resolve"
            );
        }
    }
}
