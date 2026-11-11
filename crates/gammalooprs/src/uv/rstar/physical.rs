//! Complete original-forest quotients grouped before physical integration.

use std::{
    sync::{
        Arc, LazyLock,
        atomic::{AtomicU64, Ordering},
    },
    time::Instant,
};

use eyre::{Result, ensure, eyre};
use idenso::{CookMode, CookSettings, tensor::SymbolicTensor};
use spenso::structure::{
    partial::PartialStructureExt,
    slot::{SlotMatch, SlotMatcher},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::CoefficientView,
};
use vakint::Vakint;

use crate::{
    model::Model,
    utils::GS,
    uv::{
        UVgenerationSettings, approx::integrated::Integrated, counterterm_identity::CountertermKey,
    },
};

use super::expression::{RStarIntegrand, RootProjector};

static MOMENTUM_HEADS: LazyLock<[Symbol; 7]> = LazyLock::new(|| {
    [
        GS.emr_mom,
        GS.loop_mom,
        GS.external_mom,
        feynkit_graph::symbols::loop_momentum(),
        feynkit_graph::symbols::external_momentum(),
        vakint::symbols::S.k,
        vakint::symbols::S.p,
    ]
});

// This discriminator is only an in-memory grouping key. It never enters a
// persisted graph/cohort identity, and prevents unsupported tensors from
// sharing a representative even when their literal spellings happen to match.
static NEXT_SINGLETON: AtomicU64 = AtomicU64::new(0);

/// Exact complete quotient integrand, up to a momentum-independent scalar.
/// Settings and remaining loop order are part of the equivalence certificate.
/// External incoming/outgoing roles stay fixed; external-leg exchanges and
/// unresolved internal automorphisms may conservatively form separate groups.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct PhysicalUvKey {
    integral: Option<CountertermKey>,
    algebra: Option<Atom>,
    singleton: Option<u64>,
    settings: Arc<str>,
    loops: usize,
}

/// One original proper-UV spinney with its complete child-counterterm product.
/// No loop tensor reduction or physical integration has taken place yet.
pub struct PreparedPhysicalUvTerm {
    remainder: Arc<RStarIntegrand>,
    core: Atom,
    weight: Atom,
    key: PhysicalUvKey,
    epsilon_terms: i64,
}

impl PreparedPhysicalUvTerm {
    pub(super) fn new(
        remainder: RStarIntegrand,
        model: &Model,
        settings: Arc<str>,
        projector: Option<&RootProjector<'_>>,
    ) -> Result<Self> {
        ensure!(
            remainder.region.loops <= 4,
            "physical quotient grouping supports at most four remaining loops"
        );
        let numerator = remainder.physical_numerator();
        let numerator = if let Some(projector) = projector {
            projector(&numerator)?
        } else {
            numerator
        };
        // Projectors may spell Lorentz-slot dimensions as 4-2*epsilon.
        // The physical reducer uses symbolic D slots; restore only metadata,
        // retaining every scalar evanescent projector coefficient inside K.
        let numerator = Integrated::restore_lorentz_dimensions(&numerator);
        let numerator = RStarIntegrand::simplify(&numerator)?;
        let (core, weight) = Self::factor_scalar(&numerator)?;
        let loops = remainder.region.loops;
        let (integral, algebra, singleton) = if loops == 0 {
            // There is no denominator graph to canonize. Exact expression
            // equality is conservative and preserves its original open frame.
            (None, Some(core.clone()), None)
        } else if let Some((integral, canonical_core)) =
            remainder.physical_identity(model, &core)?
        {
            (Some(integral), Some(canonical_core), None)
        } else {
            let singleton = NEXT_SINGLETON
                .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |value| {
                    value.checked_add(1)
                })
                .map_err(|_| eyre!("physical quotient singleton identifiers exhausted"))?;
            (None, None, Some(singleton))
        };
        Ok(Self {
            remainder: Arc::new(remainder),
            core,
            weight,
            key: PhysicalUvKey {
                integral,
                algebra,
                singleton,
                settings,
                loops,
            },
            epsilon_terms: i64::try_from(loops)?,
        })
    }

    pub fn key(&self) -> &PhysicalUvKey {
        &self.key
    }

    /// Attach an original branch factor exactly once. Laurent poles are
    /// allowed: this factor stays inside the subsequent loop integration.
    pub fn with_weight(mut self, weight: &Atom) -> Result<Self> {
        ensure!(
            Self::is_scalar_weight(weight)?,
            "physical quotient weights must be momentum-independent scalars"
        );
        self.weight *= weight;
        Ok(self)
    }

    fn factor_scalar(expression: &Atom) -> Result<(Atom, Atom)> {
        let factors = match expression.as_view() {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            value => vec![value],
        };
        let mut core = Atom::one();
        let mut weight = Atom::one();
        for factor in factors {
            let factor = factor.to_owned();
            if Self::is_scalar_weight(&factor)? {
                weight *= factor;
            } else {
                core *= factor;
            }
        }
        Ok((core, weight))
    }

    fn is_scalar_weight(expression: &Atom) -> Result<bool> {
        if MOMENTUM_HEADS
            .iter()
            .copied()
            .any(|momentum| expression.contains_symbol(momentum))
        {
            return Ok(false);
        }
        // A closed color/Dirac tensor can still own private summed indices.
        // Keep every such factor in the certified core; moving it between
        // isomorphic occurrences without transporting it could alias a dummy.
        let mut slots = SlotMatcher::default();
        let mut tensor_syntax = false;
        expression.visitor(&mut |part| {
            if !matches!(slots.classify(part), SlotMatch::Other) {
                tensor_syntax = true;
                false
            } else {
                !tensor_syntax
            }
        });
        if tensor_syntax {
            return Ok(false);
        }
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        Ok(
            SymbolicTensor::infer(cooking.try_cook(expression.as_view())?)?
                .structure()
                .logical_slots()
                .is_empty(),
        )
    }
}

/// Matching complete quotient cores evaluated once with their summed weights.
///
/// The certificate proves equality of both numerator and denominators under
/// one graph isomorphism. Consequently the representative's tensor frame and
/// momentum basis can be retained without transporting or factorizing sites.
/// This value is never passed to a new UV forest.
pub struct PhysicalUvGroup {
    representative: PreparedPhysicalUvTerm,
    weight: Atom,
    members: usize,
    epsilon_terms: i64,
}

impl PhysicalUvGroup {
    pub fn new(members: Vec<PreparedPhysicalUvTerm>) -> Result<Self> {
        let mut members = members.into_iter();
        let representative = members
            .next()
            .ok_or_else(|| eyre!("a physical quotient group needs at least one member"))?;
        let mut weight = representative.weight.clone();
        let mut count = 1;
        let mut epsilon_terms = representative.epsilon_terms;
        for member in members {
            ensure!(
                member.key == representative.key,
                "physical quotient integrands or settings do not match"
            );
            weight += member.weight;
            epsilon_terms = epsilon_terms.max(member.epsilon_terms);
            count += 1;
        }
        Ok(Self {
            representative,
            weight: Self::collect_weight(&weight),
            members: count,
            epsilon_terms,
        })
    }

    fn collect_weight(weight: &Atom) -> Atom {
        // together() already expands rational numerators, preserving factored
        // denominators. Noninteger powers and inexact numbers have different
        // expansion semantics, so retain the original behavior for those.
        let mut rational = true;
        weight.visitor(&mut |part| {
            if !rational {
                return false;
            }
            match part {
                AtomView::Fun(_) => return false,
                AtomView::Num(number) => {
                    rational = matches!(
                        number.get_coeff_view(),
                        CoefficientView::Natural(..) | CoefficientView::Large(..)
                    );
                }
                AtomView::Pow(power) => rational = i64::try_from(power.get_exp()).is_ok(),
                _ => {}
            }
            rational
        });
        if rational {
            weight.together()
        } else {
            weight.expand().together()
        }
    }

    pub fn members(&self) -> usize {
        self.members
    }

    pub fn loops(&self) -> usize {
        self.representative.key.loops
    }

    pub fn is_zero(&self) -> bool {
        self.weight.is_zero() || self.representative.core.is_zero()
    }

    /// Return the positive physical UV-preparation contribution in the first
    /// member's external frame. Any supplied projector and branch weights are
    /// already included. The caller still owns the final minus sign and K.
    pub fn integrate(&self, vakint: &Vakint, uv_settings: &UVgenerationSettings) -> Result<Atom> {
        ensure!(
            self.representative.key.settings.as_ref()
                == serde_json::to_string(&uv_settings.vakint)?,
            "physical quotient settings changed after preparation"
        );
        let started = Instant::now();
        let zero = self.is_zero();
        let result = if zero {
            Atom::Zero
        } else {
            let mut settings = uv_settings.vakint.true_settings();
            settings.number_of_terms_in_epsilon_expansion = self.epsilon_terms;
            self.representative.remainder.integrate_physical_numerator(
                &self.weight * &self.representative.core,
                vakint,
                uv_settings,
                &settings,
                1,
            )?
        };
        tracing::info!(
            uv = true,
            rstar = true,
            timing = true,
            stage = "rstar_grouped_physical_quotient",
            prepared_quotients = self.members,
            canonical_groups = 1,
            nonzero_integrations = usize::from(!zero),
            physical_backend_calls = usize::from(!zero && self.loops() > 0),
            algebra_only = !zero && self.loops() == 0,
            loops = self.loops(),
            integration_seconds = started.elapsed().as_secs_f64(),
            "Evaluated a group of complete original-forest physical quotients"
        );
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        finalized_runtime_dot,
        graph::{Autogen, Graph, parse::IntoFinalizedRuntimeGraph},
        initialisation::test_initialise,
        numerator::aind::Aind,
        utils::load_generic_model,
        uv::rstar::RStarPlan,
    };
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::{
        network::{library::symbolic::ETS, tags::SPENSO_TAG},
        structure::{
            representation::{Minkowski, RepName},
            slot::ParseableAind,
        },
    };
    use symbolica::function;

    #[test]
    fn physical_uv_weight_collection_preserves_rational_and_opaque_coefficients() {
        for weight in [
            symbolica::parse!("(d+1)^3/(ep^3*(d-4)^5)+(d-1)/(ep*(d-4)^2)"),
            symbolica::parse!("(1+1i)*(d+1)/(d^2+1)+(1-1i)*(d-1)/(d^2+1)"),
            symbolica::parse!("f(x*(y+z))*(d+1)/(d-4)^3-f(x*y+x*z)/(d-4)"),
        ] {
            let original = weight.expand().together();
            let collected = PhysicalUvGroup::collect_weight(&weight);
            assert!((collected - original).expand().together().is_zero());
        }
    }

    #[test]
    fn physical_uv_weight_collection_preserves_noninteger_power_expansion() {
        let weight = symbolica::parse!("(x*(y+z))^(1/2)-(x*y+x*z)^(1/2)");
        let original = weight.expand().together();
        assert!(original.is_zero());
        assert_eq!(PhysicalUvGroup::collect_weight(&weight), original);
        let symbolic = symbolica::parse!("(x*(y+z))^n-(x*y+x*z)^n");
        assert_eq!(
            PhysicalUvGroup::collect_weight(&symbolic),
            symbolic.expand().together()
        );
    }

    fn massless(graph: &mut Graph) {
        let edges = graph
            .iter_edges()
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        for edge in edges {
            graph.underlying[edge].particle =
                graph[edge].particle.clone().override_mass(Some(Atom::Zero));
        }
    }

    fn nested_bubble(model: &Model, outer_weight: Atom) -> Result<Graph> {
        let mut graph: Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1]; node [num=1]; edge [particle="scalar_1" num=1];
                e [style=invis];
                e -> A:0 [id=4 sink="{ufo_order:0}"];
                C:1 -> e [id=5 source="{ufo_order:0}"];
                A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:0}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:1}"];
                A -> C [id=2 source="{ufo_order:3}" sink="{ufo_order:1}"];
                B -> C [id=3 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, model
        )?;
        massless(&mut graph);
        let outer = graph
            .underlying
            .iter_nodes()
            .find_map(|(node, _, data)| (data.name.value == "C").then_some(node))
            .unwrap();
        graph.underlying[outer].num = Autogen::explicit(outer_weight);
        Ok(graph)
    }

    fn quotient(graph: &Graph, model: &Model, marked: usize) -> Result<RStarIntegrand> {
        let settings = super::super::tests::settings();
        let plan = RStarPlan::new(graph, model, EdgeIndex::from(marked), &settings)?;
        let term = plan
            .ordinary_uv
            .iter()
            .find(|term| !term.counterterms.is_empty())
            .unwrap();
        RStarIntegrand::new(
            graph,
            model,
            &term.remainder,
            &vec![Atom::one(); term.counterterms.len()],
            plan.marked_edge,
            Atom::Zero,
        )
    }

    #[test]
    fn physical_uv_quotient_identity_tracks_reversed_momenta_and_consumed_edges() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut first = nested_bubble(&model, Atom::one())?;
        first.underlying[EdgeIndex::from(3)]
            .extra_data
            .vakint_edge_power = Some(2);
        let mut second: Graph = finalized_runtime_dot!(
            digraph H {
                graph [projector=1]; node [num=1]; edge [particle="scalar_1" num=1];
                e [style=invis]; Y; Z; X;
                e -> X:0 [id=0 sink="{ufo_order:0}"];
                Y:1 -> e [id=1 source="{ufo_order:0}"];
                X -> Z [id=2 source="{ufo_order:1}" sink="{ufo_order:0}"];
                Z -> X [id=3 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:2}"];
                Z -> Y [id=4 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:2}"];
                Y -> X [id=5 source="{ufo_order:1}" sink="{ufo_order:3}"];
            }, &model
        )?;
        massless(&mut second);
        second.underlying[EdgeIndex::from(4)]
            .extra_data
            .vakint_edge_power = Some(2);
        let mink = Minkowski {}.new_rep(GS.dim);
        let mu = mink.to_symbolic([Atom::from(Aind::Dummy(550))]);
        let nu = mink.to_symbolic([Atom::from(Aind::Dummy(551))]);
        let first_core = function!(GS.emr_mom, 2, &mu) * function!(GS.emr_mom, 4, &mu);
        let second_core = -function!(GS.emr_mom, 5, &nu) * function!(GS.emr_mom, 0, &nu);
        let first_key = quotient(&first, &model, 0)?
            .physical_identity(&model, &first_core)?
            .expect("closed quotient core has a compact identity");
        let second_key = quotient(&second, &model, 2)?
            .physical_identity(&model, &second_core)?
            .expect("renamed closed quotient core has a compact identity");
        assert_eq!(first_key, second_key);
        // Already consumed local rules are not part of this physical integral.
        first.underlying[EdgeIndex::from(0)].num = Autogen::explicit(Atom::num(37));
        assert_eq!(
            first_key,
            quotient(&first, &model, 0)?
                .physical_identity(&model, &first_core)?
                .expect("consumed edges do not change the supported core")
        );
        // A second scalar contraction needs its own dummy index; squaring the
        // explicit-index syntax would otherwise reuse mu on four tensor ports.
        let squared_core =
            &first_core * function!(GS.emr_mom, 2, &nu) * function!(GS.emr_mom, 4, &nu);
        assert!(
            quotient(&first, &model, 0)?
                .physical_identity(&model, &first_core.pow(2))
                .is_err(),
            "an index on four compatible ports must remain an error"
        );
        assert_ne!(
            first_key,
            quotient(&first, &model, 0)?
                .physical_identity(&model, &squared_core)?
                .expect("independent closed contractions have a compact identity")
        );
        first.underlying[EdgeIndex::from(3)]
            .extra_data
            .vakint_edge_power = Some(3);
        assert_ne!(
            first_key,
            quotient(&first, &model, 0)?
                .physical_identity(&model, &first_core)?
                .expect("a changed denominator power still has a compact identity")
        );
        Ok(())
    }

    #[test]
    fn physical_uv_expanded_numerator_keeps_identity_graph_compact() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let graph = nested_bubble(&model, Atom::one())?;
        let remainder = quotient(&graph, &model, 0)?;
        let mink = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
        let scalar = function!(
            SPENSO_TAG.dot,
            function!(GS.emr_mom, 2, &mink),
            function!(GS.emr_mom, 4, &mink)
        );
        let expanded = (1..=128).fold(Atom::Zero, |sum, power| {
            sum + Atom::num(power) * scalar.pow(power)
        });
        let small = remainder
            .physical_identity(&model, &scalar)?
            .expect("closed scalar core is supported");
        let large = remainder
            .physical_identity(&model, &expanded)?
            .expect("large closed scalar core is supported");
        assert_eq!(small.0, large.0);
        assert_eq!(small.0.incidence_nodes(), large.0.incidence_nodes());
        assert!(large.0.incidence_nodes() < 32);
        assert_ne!(small.1, large.1);
        assert!(matches!(large.1.as_view(), AtomView::Add(sum) if sum.iter().count() == 128));
        let unsupported_momentum = function!(vakint::symbols::S.k, 1, &mink);
        assert!(
            remainder
                .physical_identity(
                    &model,
                    &function!(SPENSO_TAG.dot, &unsupported_momentum, &unsupported_momentum),
                )?
                .is_none(),
            "a momentum convention without a certified edge map cannot share a key"
        );
        Ok(())
    }

    #[test]
    fn physical_uv_open_cores_use_distinct_singletons_and_still_integrate() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let graph = nested_bubble(&model, Atom::one())?;
        let settings = super::super::tests::settings();
        let vakint = Vakint::new()?;
        let mink = Minkowski {}.new_rep(GS.dim);
        let external = function!(
            GS.emr_mom,
            4,
            mink.to_symbolic([Atom::from(Aind::Dummy(721))])
        );
        let project = |numerator: &Atom| Ok(numerator * &external);
        let prepare = || {
            PreparedPhysicalUvTerm::new(
                quotient(&graph, &model, 0)?,
                &model,
                serde_json::to_string(&settings.vakint)?.into(),
                Some(&project),
            )
        };
        let first = prepare()?;
        let second = prepare()?;
        assert!(first.key.singleton.is_some());
        assert_ne!(first.key(), second.key());
        assert!(PhysicalUvGroup::new(vec![first, second]).is_err());
        let mut backend = settings.vakint.true_settings();
        backend.number_of_terms_in_epsilon_expansion = 1;
        let scalar =
            quotient(&graph, &model, 0)?.integrate_physical(&vakint, &settings, &backend, 1)?;
        let grouped = PhysicalUvGroup::new(vec![prepare()?])?.integrate(&vakint, &settings)?;
        assert!((grouped - scalar * external).expand().is_zero());
        Ok(())
    }

    #[test]
    fn physical_uv_powered_propagator_matches_direct_scalar_oracle() -> Result<()> {
        use vakint::symbols::S;

        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut graph: Graph = finalized_runtime_dot!(
            digraph bubble {
                graph [projector=1]; node [num=1]; edge [particle="scalar_1" num=1];
                e [style=invis];
                e -> A:0 [id=2 sink="{ufo_order:0}"];
                B:1 -> e [id=3 source="{ufo_order:0}"];
                A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, &model
        )?;
        massless(&mut graph);
        graph.underlying[EdgeIndex::from(0)]
            .extra_data
            .vakint_edge_power = Some(2);
        let settings = super::super::tests::settings();
        let vakint = Vakint::new()?;
        let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
        let term = plan
            .ordinary_uv
            .iter()
            .find(|term| term.counterterms.is_empty())
            .unwrap();
        let remainder = RStarIntegrand::new(
            &graph,
            &model,
            &term.remainder,
            &[],
            plan.marked_edge,
            Atom::Zero,
        )?;
        let mut backend = settings.vakint.true_settings();
        backend.number_of_terms_in_epsilon_expansion = 1;
        backend.use_dot_product_notation = true;
        let actual = remainder.integrate_physical(&vakint, &settings, &backend, 7001)?;
        // Construct the scalar integral independently of Graph edge metadata.
        let input = function!(
            S.topo,
            function!(S.prop, 1, function!(S.edge, 0, 1), function!(S.k, 1), 0, 2)
                * function!(
                    S.prop,
                    2,
                    function!(S.edge, 1, 0),
                    function!(S.k, 1) + function!(S.p, 2),
                    0,
                    1
                )
        );
        let expected = vakint.evaluate_integral(&backend, input.as_view())?;
        let expected = RStarIntegrand::finish_tensor(expected, 7002)?;
        super::super::tests::assert_same_counterterm(&actual, &expected)?;
        Ok(())
    }

    #[test]
    fn physical_uv_groups_preserve_original_forests_and_laurent_weights() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut settings = super::super::tests::settings();
        settings.vakint.normalization = "GScheme".into();
        let vakint = Vakint::new()?;
        let epsilon = Atom::var(GS.dim_epsilon);
        let mut terms = Vec::new();
        let mut expected = Atom::Zero;
        for weight in [Atom::one() / &epsilon, Atom::num(2) + &epsilon] {
            let graph = nested_bubble(&model, weight)?;
            let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
            expected +=
                plan.counterterm_from_pole_part(&graph, &model, &vakint, &settings, &Atom::Zero)?;
            let prepared =
                plan.prepare_physical_uv_terms(&graph, &model, &vakint, &settings, None)?;
            assert_eq!(prepared.len(), 1);
            terms.extend(prepared);
        }
        assert_eq!(terms[0].key(), terms[1].key());
        let group = PhysicalUvGroup::new(terms)?;
        assert_eq!(group.members(), 2);
        assert_eq!(group.loops(), 1);
        let result = -group
            .integrate(&vakint, &settings)?
            .series(GS.dim_epsilon, 0, -1)?
            .to_atom();
        super::super::tests::assert_same_counterterm(&result, &expected)?;
        // The pole in the outer scalar weight needs a deeper scalar bubble
        // expansion; projecting its unweighted integral first is insufficient.
        assert!(!result.is_zero());
        Ok(())
    }

    #[test]
    fn physical_uv_group_cancellation_keeps_preparation_and_skips_integration() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let settings = super::super::tests::settings();
        let vakint = Vakint::new()?;
        let graph = nested_bubble(&model, Atom::one())?;
        let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
        let mut terms = Vec::new();
        for weight in [Atom::one(), Atom::num(-1)] {
            for term in plan.prepare_physical_uv_terms(&graph, &model, &vakint, &settings, None)? {
                terms.push(term.with_weight(&weight)?);
            }
        }
        let group = PhysicalUvGroup::new(terms)?;
        assert_eq!(group.members(), 2);
        assert!(group.is_zero());
        assert!(group.integrate(&vakint, &settings)?.is_zero());
        let mut changed = settings.clone();
        changed.vakint.additional_normalization = "2".into();
        assert!(group.integrate(&vakint, &changed).is_err());
        Ok(())
    }

    #[test]
    fn physical_uv_scalar_weights_reject_hidden_momenta_and_private_tensor_indices() -> Result<()> {
        use idenso::IndexTooling;

        test_initialise()?;
        let mink = Minkowski {}.new_rep(GS.dim);
        let index = mink.to_symbolic([Aind::Hedge(4, 1).to_atom()]);
        let other = mink.to_symbolic([Aind::Hedge(5, 1).to_atom()]);
        let vector = function!(GS.emr_mom, 0, mink.to_symbolic([]));
        let square = function!(SPENSO_TAG.dot, &vector, &vector);
        let closed_tensor = function!(ETS.metric, &index, &other).pow(2);
        assert!(closed_tensor.has_repeated_lorentz_indices());
        for factor in [
            square.pow(-1),
            function!(Symbol::LOG, square),
            closed_tensor,
        ] {
            assert!(
                !PreparedPhysicalUvTerm::is_scalar_weight(&factor)?,
                "a factor retaining momentum or tensor ports must stay in the core: {factor}"
            );
        }
        // g(mu,mu) normalizes to scalar D during construction, so it no
        // longer owns any private tensor index and may safely be extracted.
        assert!(PreparedPhysicalUvTerm::is_scalar_weight(&function!(
            ETS.metric, &index, &index
        ))?);
        let weight = Atom::num(2) / Atom::var(GS.dim_epsilon).pow(3);
        assert!(PreparedPhysicalUvTerm::is_scalar_weight(&weight)?);
        Ok(())
    }
}
