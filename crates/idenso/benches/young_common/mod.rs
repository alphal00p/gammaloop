use idenso::{
    IndexTooling,
    reference_cases::young::{
        FactoredRiemannProjector, YoungProjectorFixture, distinct_head_contracted_riemann_fixtures,
        fully_projected_distinct_head_riemann_triangle_fixture, young_projector_fixtures,
    },
};
use spenso::structure::abstract_index::AbstractIndex;
use symbolica::atom::{Atom, AtomCore, FunctionBuilder};

pub struct ValidatedCorpus {
    pub projectors: Vec<YoungProjectorFixture>,
    pub product_expansion: Vec<CanonicalizationCase>,
    pub canonicalization: Vec<CanonicalizationCase>,
}

pub struct CanonicalizationCase {
    pub name: &'static str,
    pub expression: Atom,
}

fn tensor(head: symbolica::atom::Symbol, arguments: &[Atom]) -> Atom {
    arguments
        .iter()
        .fold(FunctionBuilder::new(head), |builder, argument| {
            builder.add_arg(argument)
        })
        .finish()
}

fn validate_canonicalization(expression: &Atom) -> Atom {
    let canonical = canonicalize(expression.clone());
    assert_eq!(
        canonicalize(canonical.clone()),
        canonical,
        "Young-projector benchmark output must be a fixed point"
    );
    canonical
}

pub fn validated_corpus() -> ValidatedCorpus {
    let projectors = young_projector_fixtures();
    let riemann = projectors
        .iter()
        .find(|fixture| fixture.projector.tableau().shape() == [2, 2])
        .expect("the Young-projector corpus must contain the Riemann fixture");
    let factored_riemann = FactoredRiemannProjector.project(
        FactoredRiemannProjector::tensor_symbol("young_benchmark_factored_riemann"),
        &riemann.arguments,
    );
    let controls = [
        ("ordered_rank_4", riemann.ordered.clone()),
        (
            "symmetric_rank_4",
            tensor(
                spenso::tensor_symbol!(young_benchmark_symmetric_control; Symmetric),
                &riemann.arguments,
            ),
        ),
        (
            "antisymmetric_rank_4",
            tensor(
                spenso::tensor_symbol!(young_benchmark_antisymmetric_control; Antisymmetric),
                &riemann.arguments,
            ),
        ),
    ];
    let mut canonicalization = Vec::new();

    for fixture in &projectors {
        assert_eq!(
            fixture.projector.tableau().rank(),
            fixture.arguments.len(),
            "fixture {} has the wrong argument count",
            fixture.name
        );
        assert!(fixture.projector.hook_product() > 0);
        assert_eq!(
            fixture.projector.project(fixture.head, &fixture.arguments),
            fixture.projected,
            "fixture {} does not match its projector",
            fixture.name
        );

        // The full ordered `[2, 2]` oracle is retained for compilation and
        // expansion timings, but crosses the current whole-graph edge limit.
        // Its exact six-coset structural reduction supplies canonicalization
        // coverage without changing the global graph budget.
        if fixture.projector.tableau().shape() != [2, 2] {
            assert!(!validate_canonicalization(&fixture.projected).is_zero());
            canonicalization.push(CanonicalizationCase {
                name: fixture.name,
                expression: fixture.projected.clone(),
            });
        }
    }
    assert!(!validate_canonicalization(&factored_riemann).is_zero());
    canonicalization.push(CanonicalizationCase {
        name: "riemann_2_2_factored",
        expression: factored_riemann,
    });

    for fixture in distinct_head_contracted_riemann_fixtures() {
        assert_ne!(fixture.expression, fixture.expression.expand());
        let canonical = validate_canonicalization(&fixture.expression);
        let renamed = canonicalize(fixture.renamed);
        assert!(
            !canonical.is_zero(),
            "fixture {} must be nonzero",
            fixture.name
        );
        assert_eq!(
            canonical, renamed,
            "fixture {} must ignore dummy-index names",
            fixture.name
        );
        canonicalization.push(CanonicalizationCase {
            name: fixture.name,
            expression: fixture.expression,
        });
    }

    let fully_projected_triangle = fully_projected_distinct_head_riemann_triangle_fixture();
    assert_ne!(
        fully_projected_triangle.expression,
        fully_projected_triangle.expression.expand()
    );
    let product_expansion = vec![CanonicalizationCase {
        name: fully_projected_triangle.name,
        expression: fully_projected_triangle.expression,
    }];

    for (name, control) in controls {
        assert!(!validate_canonicalization(&control).is_zero());
        canonicalization.push(CanonicalizationCase {
            name,
            expression: control,
        });
    }

    ValidatedCorpus {
        projectors,
        product_expansion,
        canonicalization,
    }
}

pub fn canonicalize(expression: Atom) -> Atom {
    expression
        .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
        .expect("validated Young-projector fixture must canonicalize")
}
