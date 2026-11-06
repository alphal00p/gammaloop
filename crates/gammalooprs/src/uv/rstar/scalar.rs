#[cfg(test)]
mod tests {
    use crate::{
        finalized_runtime_dot,
        graph::parse::IntoFinalizedRuntimeGraph,
        initialisation::test_initialise,
        utils::{GS, load_generic_model},
        uv::{
            ApproximationType, RenormalizationPrescriptionSettings,
            settings::VakintSettings as GammaVakintSettings,
        },
    };
    use crate::{
        graph::Graph,
        uv::{UVgenerationSettings, rstar::RStarPlan},
    };
    use eyre::Result;
    use linnet::half_edge::involution::EdgeIndex;
    use symbolica::atom::{Atom, AtomCore};
    use vakint::Vakint;

    #[test]
    fn two_loop_counterterm_agrees_for_hard_line_inside_and_outside_uv_bubble() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                edge [particle="scalar_1" num=1];
                e [style=invis];
                e -> A:0 [id=4 sink="{ufo_order:0}"];
                B:1 -> e [id=5 source="{ufo_order:0}"];
                A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"];
                A -> C [id=2 source="{ufo_order:3}" sink="{ufo_order:0}"];
                C -> B [id=3 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:3}"];
            }, "scalars"
        )?;
        let edges: Vec<_> = graph.iter_edges().map(|(_, edge, _)| edge).collect();
        for edge in edges {
            let data = &mut graph.underlying[edge];
            data.particle = data.particle.clone().override_mass(Some(Atom::Zero));
        }
        let model = load_generic_model("scalars");
        let settings = UVgenerationSettings {
            renormalization_prescription: RenormalizationPrescriptionSettings {
                log_divergent: ApproximationType::PolePart,
                ..Default::default()
            },
            vakint: GammaVakintSettings {
                evaluation_methods: vec!["alphaloop".into(), "forcer".into()],
                normalization: "FMFTandMATAD".into(),
                additional_normalization: "1".into(),
                ..Default::default()
            },
            ..Default::default()
        };
        let vakint = Vakint::new()?;
        let inside = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
        let outside = RStarPlan::new(&graph, &model, EdgeIndex::from(2), &settings)?;
        assert!(!inside.ir.is_empty());
        assert!(outside.ir.is_empty());
        let inside_result = inside.counterterm(&graph, &model, &vakint, &settings)?;
        let outside_result = outside.counterterm(&graph, &model, &vakint, &settings)?;
        assert!(
            (&inside_result - &outside_result).expand().is_zero(),
            "inside: {inside_result}; outside: {outside_result}"
        );
        let epsilon = Atom::var(GS.dim_epsilon);
        let expected = Atom::num(1) / (Atom::num(2) * epsilon.pow(2))
            - Atom::num(1) / (Atom::num(2) * epsilon);
        assert!(
            (&inside_result - &expected).expand().is_zero(),
            "R*: {inside_result}; expected: {expected}"
        );
        // This test may evaluate the old two-loop vacuum as an independent
        // oracle; the production R* route only requests one-loop UV children.
        let mut ordinary = inside.uv_forest(&graph, &settings);
        ordinary.integrate(&graph, &model, &vakint, &settings)?;
        let old_result = ordinary.counterterm_for_support(&inside.support)?;
        assert!((&inside_result - old_result).expand().is_zero());
        Ok(())
    }

    #[test]
    fn three_loop_nested_counterterms_agree_with_massive_r() -> Result<()> {
        test_initialise()?;
        // The AB bubble lies inside the logarithmic ABC component. Marking
        // AB exercises a nonempty strict ancestor chain; marking AD does not.
        let mut graph: Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                edge [particle="scalar_1" num=1];
                e [style=invis];
                e -> A:0 [id=6 sink="{ufo_order:0}"];
                C:1 -> e [id=7 source="{ufo_order:0}"];
                A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:0}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:1}"];
                A -> C [id=2 source="{ufo_order:3}" sink="{ufo_order:1}"];
                C -> B [id=3 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:2}"];
                A -> D [id=4 source="{ufo_order:4}" sink="{ufo_order:0}"];
                D -> C [id=5 lmb_id=2 source="{ufo_order:1}" sink="{ufo_order:3}"];
            }, "scalars"
        )?;
        let edges: Vec<_> = graph.iter_edges().map(|(_, edge, _)| edge).collect();
        for edge in edges {
            let data = &mut graph.underlying[edge];
            data.particle = data.particle.clone().override_mass(Some(Atom::Zero));
        }
        let model = load_generic_model("scalars");
        let settings = UVgenerationSettings {
            renormalization_prescription: RenormalizationPrescriptionSettings {
                log_divergent: ApproximationType::PolePart,
                ..Default::default()
            },
            vakint: GammaVakintSettings {
                evaluation_methods: vec!["alphaloop".into(), "matad".into(), "forcer".into()],
                normalization: "FMFTandMATAD".into(),
                additional_normalization: "1".into(),
                ..Default::default()
            },
            ..Default::default()
        };
        let vakint = Vakint::new()?;
        let inside = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
        let outside = RStarPlan::new(&graph, &model, EdgeIndex::from(4), &settings)?;
        assert!(
            inside
                .ir
                .iter()
                .flat_map(|term| &term.forest)
                .any(|forest| { forest.steps.len() > 1 })
        );
        assert!(outside.ir.is_empty());
        assert!(
            inside
                .counterterms
                .iter()
                .any(|call| call.argument.loops == 2)
        );
        let inside_result = inside.counterterm(&graph, &model, &vakint, &settings)?;
        let outside_result = outside.counterterm(&graph, &model, &vakint, &settings)?;
        assert!(
            (&inside_result - &outside_result).expand().is_zero(),
            "inside: {inside_result}; outside: {outside_result}"
        );
        // Only this oracle evaluates the full three-loop massive vacuum.
        // Its nested poles need more one-loop master orders than alphaLoop's
        // table supplies, so use MATAD for the independent massive oracle.
        let mut oracle_settings = settings.clone();
        oracle_settings.vakint.evaluation_methods = vec!["matad".into()];
        oracle_settings.vakint.matad.direct_numerical_substition = false;
        oracle_settings.vakint.matad.substitute_hpls = false;
        let mut ordinary = inside.uv_forest(&graph, &oracle_settings);
        ordinary.integrate(&graph, &model, &vakint, &oracle_settings)?;
        let old_result = ordinary.counterterm_for_support(&inside.support)?;
        assert!(
            (&inside_result - &old_result).expand().is_zero(),
            "R*: {inside_result}; R: {old_result}"
        );
        Ok(())
    }
}
