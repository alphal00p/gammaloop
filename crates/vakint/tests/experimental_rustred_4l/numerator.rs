//! Cheap scalar-transport regressions: no rule search, FORM, or master table.

use rustred::family::IntegralKey;
use rustred::scalar_numerator::{FamilyScalarNumeratorService, ScalarNumeratorLimits};
use symbolica::atom::{Atom, AtomCore};
use vakint::rustred_evaluation::experimental::{
    CandidateLoweredTerm, CandidateReduction, CandidateScalarReduction, prepare_candidate_integrals,
};
use vakint::symbols::S;
use vakint::{Vakint, VakintSettings, vakint_parse as vk_parse};

use super::input::ParentInput;

struct LoweringOnly {
    parent: ParentInput,
    dimension: Atom,
}

impl std::fmt::Debug for LoweringOnly {
    fn fmt(&self, out: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        out.write_str("LoweringOnly(no rule search or application)")
    }
}

impl CandidateScalarReduction for LoweringOnly {
    fn index_count(&self) -> usize {
        self.parent.family.denominator_count()
    }

    fn parent_momenta(&self) -> &[Atom] {
        &self.parent.physical_momenta
    }

    fn dimension(&self) -> &Atom {
        &self.dimension
    }

    fn reduce_unit_mass(&self, _: &IntegralKey) -> Result<CandidateReduction, String> {
        Err("lower-only regression must never apply an IBP".into())
    }

    fn lower_scalar_numerator(
        &self,
        numerator: &Atom,
        base: &IntegralKey,
    ) -> Result<Vec<CandidateLoweredTerm>, String> {
        assert!(
            !numerator.contains_symbol(S.k),
            "Vakint loop-vector components must not reach the family lowerer: {numerator}"
        );
        let loops = (1..=4)
            .map(|axis| vk_parse!(format!("k{axis}")).unwrap())
            .collect::<Vec<_>>();
        let service = FamilyScalarNumeratorService::try_new(
            &self.parent.family,
            S.dot,
            loops.clone(),
            ScalarNumeratorLimits::default(),
        )
        .map_err(|error| error.to_string())?;
        let lowered = service
            .lower(numerator, base)
            .map_err(|error| error.to_string())?;
        Ok(lowered
            .terms()
            .iter()
            .map(|term| {
                for momentum in &loops {
                    assert!(
                        !term
                            .scalar_spectator()
                            .to_canonical_string()
                            .contains(momentum.to_canonical_string().as_str()),
                        "lowered scalar spectator still contains a loop momentum"
                    );
                }
                CandidateLoweredTerm {
                    target: term.integral().clone(),
                    coefficient: term.coefficient().to_expression(),
                    scalar_spectator: term.scalar_spectator().clone(),
                    common_mass_squared_power: term.common_mass_squared_power(),
                }
            })
            .collect())
    }
}

fn setup() -> (Vakint, VakintSettings, LoweringOnly, Atom) {
    Vakint::initialize_vakint_symbols();
    let parent = ParentInput::from_csv(include_str!("../inputs/experimental_four_loop_fg.csv"));
    let dimension = parent.family.dimension().to_expression();
    (
        Vakint::new().unwrap(),
        VakintSettings {
            form_exe_path: "/lowering-regression-must-not-invoke-form".into(),
            use_dot_product_notation: true,
            ..VakintSettings::default()
        },
        LoweringOnly { parent, dimension },
        vk_parse!(
            "topo(prop(1,edge(1,1),k(1),muvsq,2)*prop(2,edge(1,1),k(2),muvsq,1)*prop(3,edge(1,1),k(3),muvsq,1)*prop(4,edge(1,1),k(4),muvsq,1))"
        )
        .unwrap(),
    )
}

#[test]
fn clover_parent_routing_lowers_loop_dots_in_both_scalar_notations() {
    let (vakint, settings, reducer, clover) = setup();
    let mut outputs = Vec::new();
    for scalar in ["dot(k(3),k(3))", "k(3,11)*k(3,11)"] {
        let input = vk_parse!(scalar).unwrap() * &clover;
        let prepared = prepare_candidate_integrals(&vakint, &settings, input.as_view(), &reducer)
            .expect("routed loop dot must lower into scalar integral keys");
        assert!(
            prepared.len() > 1,
            "loop dot must not survive as one spectator"
        );
        outputs.push(
            prepared
                .into_iter()
                .map(|term| term.target)
                .collect::<Vec<_>>(),
        );
    }
    assert_eq!(outputs[0], outputs[1]);
}

#[test]
fn clover_parent_routing_rejects_uncontracted_loop_components() {
    let (vakint, settings, reducer, clover) = setup();
    let input = vk_parse!("k(3,11)").unwrap() * clover;
    let error = prepare_candidate_integrals(&vakint, &settings, input.as_view(), &reducer)
        .expect_err("an uncontracted loop component is not a scalar spectator");
    assert!(error.to_string().contains("retains loop-vector components"));
}
