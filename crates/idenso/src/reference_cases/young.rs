//! Independent normalized Young-projector references for benchmarks and tests.

use itertools::Itertools;
use linnet::permutation::Permutation;
use spenso::{
    network::tags::SPENSO_TAG,
    structure::{
        YoungTableau,
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName},
        slot::IsAbstractSlot,
    },
};
use symbolica::atom::{Atom, FunctionBuilder, Symbol};

/// The normalized, manifestly column-antisymmetric Young projector
/// `C_T R_T / h_T` for a fixed tableau.
///
/// `R_T` symmetrizes slots within rows, `C_T` antisymmetrizes slots within
/// columns, and `h_T` is the product of the diagram's hook lengths. The
/// projector stores only permutations, so the tensor head supplied to
/// [`Self::project`] remains an ordinary ordered Spenso tensor.
#[derive(Clone, Debug)]
pub struct YoungProjector {
    tableau: YoungTableau,
    actions: Vec<SignedAction>,
    hook_product: usize,
}

#[derive(Clone, Debug)]
struct SignedAction {
    permutation: Permutation,
    sign: i8,
}

impl YoungProjector {
    /// Construct the full row/column permutation oracle for `tableau`.
    ///
    /// # Panics
    ///
    /// Panics if the hook product cannot be represented by `usize`.
    pub fn new(tableau: YoungTableau) -> Self {
        let rank = tableau.rank();
        let rows = tableau.rows().map(<[usize]>::to_vec).collect::<Vec<_>>();
        let columns = Self::columns(&tableau);
        let row_group = Self::subgroup(rank, &rows);
        let column_group = Self::subgroup(rank, &columns);
        let actions = row_group
            .iter()
            .flat_map(|row| {
                column_group.iter().map(move |column| SignedAction {
                    // Pullback selectors compose row first, then column:
                    // combined[j] = column[row[j]].
                    permutation: column.compose(row),
                    sign: column.sign(),
                })
            })
            .collect();

        Self {
            hook_product: Self::hook_product_for(tableau.shape()),
            tableau,
            actions,
        }
    }

    pub fn tableau(&self) -> &YoungTableau {
        &self.tableau
    }

    /// Number of distinct signed terms in the unnormalized `C_T R_T` sum.
    pub fn term_count(&self) -> usize {
        self.actions.len()
    }

    pub fn hook_product(&self) -> usize {
        self.hook_product
    }

    /// Apply the exact normalized projector to an ordinary ordered tensor.
    ///
    /// # Panics
    ///
    /// Panics if `arguments` does not contain one argument per tableau box.
    pub fn project(&self, head: Symbol, arguments: &[Atom]) -> Atom {
        assert_eq!(
            arguments.len(),
            self.tableau.rank(),
            "Young projector argument count must equal the tableau rank"
        );

        let numerator = self.actions.iter().fold(Atom::Zero, |sum, action| {
            let term = action
                .permutation
                .iter_slice(arguments)
                .fold(FunctionBuilder::new(head), |builder, argument| {
                    builder.add_arg(argument)
                })
                .finish();
            if action.sign < 0 {
                sum - term
            } else {
                sum + term
            }
        });
        numerator / Atom::num(self.hook_product)
    }

    fn columns(tableau: &YoungTableau) -> Vec<Vec<usize>> {
        let mut columns = vec![Vec::new(); tableau.shape()[0]];
        for row in tableau.rows() {
            for (column, &slot) in row.iter().enumerate() {
                columns[column].push(slot);
            }
        }
        columns
    }

    fn subgroup(rank: usize, groups: &[Vec<usize>]) -> Vec<Permutation> {
        groups
            .iter()
            .map(|group| (0..group.len()).permutations(group.len()))
            .multi_cartesian_product()
            .map(|permutations| {
                let mut map = (0..rank).collect::<Vec<_>>();
                for (group, permutation) in groups.iter().zip(permutations) {
                    for (&target, source) in group.iter().zip(permutation) {
                        map[target] = group[source];
                    }
                }
                Permutation::from_map(map)
            })
            .collect()
    }

    fn hook_product_for(shape: &[usize]) -> usize {
        shape
            .iter()
            .enumerate()
            .flat_map(|(row, &row_length)| {
                (0..row_length).map(move |column| {
                    let below = shape[row + 1..]
                        .iter()
                        .filter(|&&length| length > column)
                        .count();
                    row_length - column + below
                })
            })
            .try_fold(1usize, usize::checked_mul)
            .expect("Young-tableau hook product overflows usize")
    }
}

/// A named, ready-to-benchmark ordered tensor and its explicit projector sum.
#[derive(Clone, Debug)]
pub struct YoungProjectorFixture {
    pub name: &'static str,
    pub projector: YoungProjector,
    pub head: Symbol,
    pub arguments: Vec<Atom>,
    pub ordered: Atom,
    pub projected: Atom,
}

impl YoungProjectorFixture {
    fn new(
        name: &'static str,
        head_name: &str,
        shape: &[usize],
        slot_order: &[usize],
        arguments: &[Atom],
    ) -> Self {
        let projector = YoungProjector::new(
            YoungTableau::new(shape.to_vec(), slot_order.to_vec())
                .expect("static Young-projector fixture must be valid"),
        );
        let head = SPENSO_TAG.tensor_symbol(head_name);
        let arguments = arguments[..projector.tableau().rank()].to_vec();
        let ordered = arguments
            .iter()
            .fold(FunctionBuilder::new(head), |builder, argument| {
                builder.add_arg(argument)
            })
            .finish();
        let projected = projector.project(head, &arguments);
        Self {
            name,
            projector,
            head,
            arguments,
            ordered,
            projected,
        }
    }
}

/// Build the mixed-shape corpus used to establish pre-straightening baselines.
///
/// The fillings make the row and column actions visibly different and cover a
/// rank-three hook, row-heavy and column-heavy rank-four hooks, and the Riemann
/// `[2, 2]` shape.
pub fn young_projector_fixtures() -> Vec<YoungProjectorFixture> {
    crate::representations::initialize();
    let representation = Minkowski {}.new_rep(4);
    let arguments = [
        spenso::s!(young_reference_a),
        spenso::s!(young_reference_b),
        spenso::s!(young_reference_c),
        spenso::s!(young_reference_d),
    ]
    .map(|index| representation.slot::<AbstractIndex, _>(index).to_atom());

    [
        (
            "hook_2_1",
            "young_reference_hook_2_1",
            &[2, 1][..],
            &[0, 2, 1][..],
        ),
        (
            "row_heavy_3_1",
            "young_reference_row_heavy_3_1",
            &[3, 1],
            &[0, 2, 3, 1],
        ),
        (
            "column_heavy_2_1_1",
            "young_reference_column_heavy_2_1_1",
            &[2, 1, 1],
            &[0, 3, 1, 2],
        ),
        (
            "riemann_2_2",
            "young_reference_riemann_2_2",
            &[2, 2],
            &[0, 2, 1, 3],
        ),
    ]
    .into_iter()
    .map(|(name, head, shape, slot_order)| {
        YoungProjectorFixture::new(name, head, shape, slot_order, &arguments)
    })
    .collect()
}

/// The graph-native `[2, 2]` `C_T R_T / h_T` projector for a carrier that is
/// already antisymmetric in the slot pairs `(0, 1)` and `(2, 3)`.
///
/// Reducing the full sixteen actions by the right action of the two structural
/// column groups leaves six signed cosets with weights
/// `(4, 2, -2, -2, 2, 4) / 12`. This is not the four-term `R_T C_T / 12`
/// reduction: under the fixed pullback convention, `C_T` and `R_T` cannot be
/// exchanged. Keeping the columns inert and each projected factor as one
/// unexpanded six-term sum still bounds the graph size for Riemann products.
#[derive(Clone, Copy, Debug, Default)]
pub struct FactoredRiemannProjector;

impl FactoredRiemannProjector {
    const STRUCTURAL_ACTIONS: [(i8, [usize; 4]); 6] = [
        (4, [0, 1, 2, 3]),
        (2, [0, 2, 1, 3]),
        (-2, [0, 3, 1, 2]),
        (-2, [1, 2, 0, 3]),
        (2, [1, 3, 0, 2]),
        (4, [2, 3, 0, 1]),
    ];

    /// Construct the ordinary tagged carrier required by [`Self::project`].
    ///
    /// Linearity is essential: it lifts orientation signs normalized inside a
    /// structural `antisym` group to the sign of the complete tensor factor.
    pub fn tensor_symbol(name: &str) -> Symbol {
        spenso::tensor_symbol!((name); Linear)
    }

    pub const fn term_count(self) -> usize {
        Self::STRUCTURAL_ACTIONS.len()
    }

    /// Build the exact six-term structural reduction without expanding either
    /// column group.
    ///
    /// # Panics
    ///
    /// Panics unless `head` is linear or `arguments` has exactly four slots.
    pub fn project(self, head: Symbol, arguments: &[Atom]) -> Atom {
        assert!(
            head.is_linear(),
            "the factored Riemann carrier must be linear"
        );
        assert_eq!(
            arguments.len(),
            4,
            "the factored Riemann projector requires four slots"
        );

        let numerator =
            Self::STRUCTURAL_ACTIONS
                .iter()
                .fold(Atom::Zero, |sum, &(weight, selector)| {
                    let permuted = selector.map(|slot| arguments[slot].clone());
                    sum + Atom::num(weight) * Self::structural_column_carrier(head, &permuted)
                });
        numerator / Atom::num(12)
    }

    /// Build one unprojected linear carrier with structural antisymmetry sites
    /// on slot pairs `(0, 1)` and `(2, 3)`.
    ///
    /// # Panics
    ///
    /// Panics unless `head` is linear or `arguments` has exactly four slots.
    pub fn structural_column_carrier(head: Symbol, arguments: &[Atom]) -> Atom {
        assert!(
            head.is_linear(),
            "the structural Riemann carrier must be linear"
        );
        assert_eq!(
            arguments.len(),
            4,
            "the structural Riemann carrier requires four slots"
        );
        FunctionBuilder::new(head)
            .add_arg(spenso::antisym!(arguments[0].clone(), arguments[1].clone()))
            .add_arg(spenso::antisym!(arguments[2].clone(), arguments[3].clone()))
            .finish()
    }
}

/// A contracted product of factored Riemann projectors and an isomorphic copy
/// whose dummy indices use disjoint names.
#[derive(Clone, Debug)]
pub struct ContractedRiemannFixture {
    pub name: &'static str,
    pub factor_count: usize,
    pub expression: Atom,
    pub renamed: Atom,
}

impl ContractedRiemannFixture {
    const TWO_FACTOR_TOPOLOGY: [[usize; 4]; 2] = [[0, 1, 2, 3], [0, 2, 1, 3]];
    const TRIANGLE_TOPOLOGY: [[usize; 4]; 3] = [[0, 1, 2, 3], [0, 1, 4, 5], [2, 3, 4, 5]];
    const RICCI_CYCLE_TOPOLOGY: [[usize; 4]; 3] = [[0, 3, 1, 3], [1, 4, 2, 4], [2, 5, 0, 5]];

    fn new<const N: usize>(
        name: &'static str,
        heads: [Symbol; N],
        arguments: &[Atom],
        renamed_arguments: &[Atom],
        topology: [[usize; 4]; N],
    ) -> Self {
        let build = |source: &[Atom]| {
            topology
                .iter()
                .zip(heads)
                .map(|(selector, head)| {
                    let factor_arguments = selector.map(|slot| source[slot].clone());
                    FactoredRiemannProjector.project(head, &factor_arguments)
                })
                .reduce(|product, factor| product * factor)
                .expect("a contracted Riemann fixture must contain a factor")
        };

        Self {
            name,
            factor_count: N,
            expression: build(arguments),
            renamed: build(renamed_arguments),
        }
    }

    fn slots(indices: Vec<Symbol>) -> Vec<Atom> {
        crate::representations::initialize();
        let representation = Minkowski {}.new_rep(4);
        indices
            .into_iter()
            .map(|index| representation.slot::<AbstractIndex, _>(index).to_atom())
            .collect()
    }

    fn two_factor(name: &'static str, heads: [Symbol; 2]) -> Self {
        let arguments = Self::slots(vec![
            spenso::s!(young_two_a),
            spenso::s!(young_two_b),
            spenso::s!(young_two_c),
            spenso::s!(young_two_d),
        ]);
        let renamed = Self::slots(vec![
            spenso::s!(young_two_i),
            spenso::s!(young_two_j),
            spenso::s!(young_two_k),
            spenso::s!(young_two_l),
        ]);

        Self::new(name, heads, &arguments, &renamed, Self::TWO_FACTOR_TOPOLOGY)
    }

    fn triangle(name: &'static str, heads: [Symbol; 3]) -> Self {
        let arguments = Self::slots(vec![
            spenso::s!(young_triangle_a),
            spenso::s!(young_triangle_b),
            spenso::s!(young_triangle_c),
            spenso::s!(young_triangle_d),
            spenso::s!(young_triangle_e),
            spenso::s!(young_triangle_f),
        ]);
        let renamed = Self::slots(vec![
            spenso::s!(young_triangle_i),
            spenso::s!(young_triangle_j),
            spenso::s!(young_triangle_k),
            spenso::s!(young_triangle_l),
            spenso::s!(young_triangle_m),
            spenso::s!(young_triangle_n),
        ]);

        Self::new(name, heads, &arguments, &renamed, Self::TRIANGLE_TOPOLOGY)
    }

    fn ricci_cycle(name: &'static str, heads: [Symbol; 3]) -> Self {
        let arguments = Self::slots(vec![
            spenso::s!(young_ricci_a),
            spenso::s!(young_ricci_b),
            spenso::s!(young_ricci_c),
            spenso::s!(young_ricci_x),
            spenso::s!(young_ricci_y),
            spenso::s!(young_ricci_z),
        ]);
        let renamed = Self::slots(vec![
            spenso::s!(young_ricci_i),
            spenso::s!(young_ricci_j),
            spenso::s!(young_ricci_k),
            spenso::s!(young_ricci_u),
            spenso::s!(young_ricci_v),
            spenso::s!(young_ricci_w),
        ]);

        Self::new(
            name,
            heads,
            &arguments,
            &renamed,
            Self::RICCI_CYCLE_TOPOLOGY,
        )
    }
}

/// Build same-head quadratic and fully projected cubic regression fixtures.
pub fn contracted_riemann_fixtures() -> Vec<ContractedRiemannFixture> {
    let two = FactoredRiemannProjector::tensor_symbol("young_reference_riemann_two_factor");
    let triangle = FactoredRiemannProjector::tensor_symbol("young_reference_riemann_triangle");
    vec![
        ContractedRiemannFixture::two_factor("riemann_two_factor_crossed", [two; 2]),
        ContractedRiemannFixture::triangle("riemann_three_factor_triangle", [triangle; 3]),
    ]
}

/// Build successful timing controls with a distinct carrier head per factor.
///
/// Both factors in the quadratic case are projected. The cubic Ricci cycle
/// projects all three factors but repeats one contraction index inside each
/// factor, reducing each six-coset expression to two nonzero terms.
pub fn distinct_head_contracted_riemann_fixtures() -> Vec<ContractedRiemannFixture> {
    vec![
        ContractedRiemannFixture::two_factor(
            "riemann_two_factor_crossed_distinct_heads",
            [
                FactoredRiemannProjector::tensor_symbol("young_reference_riemann_two_left"),
                FactoredRiemannProjector::tensor_symbol("young_reference_riemann_two_right"),
            ],
        ),
        ContractedRiemannFixture::ricci_cycle(
            "riemann_three_factor_ricci_cycle_distinct_heads",
            [
                FactoredRiemannProjector::tensor_symbol("young_reference_riemann_ricci_first"),
                FactoredRiemannProjector::tensor_symbol("young_reference_riemann_ricci_second"),
                FactoredRiemannProjector::tensor_symbol("young_reference_riemann_ricci_third"),
            ],
        ),
    ]
}

/// Build the fully projected distinct-head triangle that currently exceeds
/// the whole-graph budget at 138 vertices and 203 edges.
pub fn fully_projected_distinct_head_riemann_triangle_fixture() -> ContractedRiemannFixture {
    ContractedRiemannFixture::triangle(
        "riemann_three_factor_triangle_fully_projected_distinct_heads",
        [
            FactoredRiemannProjector::tensor_symbol("young_reference_riemann_full_triangle_first"),
            FactoredRiemannProjector::tensor_symbol("young_reference_riemann_full_triangle_second"),
            FactoredRiemannProjector::tensor_symbol("young_reference_riemann_full_triangle_third"),
        ],
    )
}

#[cfg(test)]
mod tests {
    use spenso::{
        network::tags::SPENSO_TAG,
        structure::{YoungTableau, abstract_index::AbstractIndex},
    };
    use symbolica::{
        atom::{Atom, AtomCore, FunctionBuilder},
        symbol,
    };

    use crate::{IndexTooling, tensor::CanonicalizationError};

    use super::{
        FactoredRiemannProjector, YoungProjector, contracted_riemann_fixtures,
        distinct_head_contracted_riemann_fixtures,
        fully_projected_distinct_head_riemann_triangle_fixture, young_projector_fixtures,
    };

    fn ordered(head: symbolica::atom::Symbol, arguments: &[Atom]) -> Atom {
        arguments
            .iter()
            .fold(FunctionBuilder::new(head), |builder, argument| {
                builder.add_arg(argument)
            })
            .finish()
    }

    fn permute<const N: usize>(arguments: &[Atom; N], selector: [usize; N]) -> [Atom; N] {
        selector.map(|slot| arguments[slot].clone())
    }

    fn assert_successful_contraction(fixture: &super::ContractedRiemannFixture) {
        assert_ne!(fixture.expression, fixture.expression.expand());
        let canonical = fixture
            .expression
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();
        let renamed = fixture
            .renamed
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();
        assert_eq!(canonical, renamed, "fixture {}", fixture.name);
        assert!(!canonical.is_zero(), "fixture {}", fixture.name);
        assert_eq!(
            canonical
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap(),
            canonical,
            "fixture {}",
            fixture.name
        );
    }

    fn assert_current_graph_limit(
        expression: &Atom,
        name: &str,
        requested_vertices: usize,
        requested_edges: usize,
    ) {
        assert!(
            matches!(
                expression.try_canonize::<AbstractIndex>(AbstractIndex::Dummy),
                Err(CanonicalizationError::WholeGraphSizeLimit {
                    requested_vertices: actual_vertices,
                    requested_edges: actual_edges,
                    vertex_limit: 128,
                    edge_limit: 160,
                    ..
                }) if actual_vertices == requested_vertices && actual_edges == requested_edges
            ),
            "fixture {name}"
        );
    }

    #[test]
    fn hook_projector_uses_column_after_row_pullback_convention() {
        let projector = YoungProjector::new(YoungTableau::new(vec![2, 1], vec![0, 2, 1]).unwrap());
        let head = SPENSO_TAG.tensor_symbol("young_reference_convention");
        let arguments = [
            symbol!("young_oracle_a"),
            symbol!("young_oracle_b"),
            symbol!("young_oracle_c"),
        ]
        .map(Atom::var);

        let expected = (ordered(
            head,
            &[
                arguments[0].clone(),
                arguments[1].clone(),
                arguments[2].clone(),
            ],
        ) + ordered(
            head,
            &[
                arguments[2].clone(),
                arguments[1].clone(),
                arguments[0].clone(),
            ],
        ) - ordered(
            head,
            &[
                arguments[1].clone(),
                arguments[0].clone(),
                arguments[2].clone(),
            ],
        ) - ordered(
            head,
            &[
                arguments[2].clone(),
                arguments[0].clone(),
                arguments[1].clone(),
            ],
        )) / Atom::num(3);

        assert_eq!(projector.term_count(), 4);
        assert_eq!(projector.hook_product(), 3);
        assert_eq!(projector.project(head, &arguments), expected);
    }

    #[test]
    fn hook_projector_obeys_column_antisymmetry_and_cyclic_bianchi() {
        let projector = YoungProjector::new(YoungTableau::new(vec![2, 1], vec![0, 2, 1]).unwrap());
        let head = SPENSO_TAG.tensor_symbol("young_reference_hook_identities");
        let arguments = [
            symbol!("young_hook_a"),
            symbol!("young_hook_b"),
            symbol!("young_hook_c"),
        ]
        .map(Atom::var);
        let project = |selector| projector.project(head, &permute(&arguments, selector));
        let canonical = project([0, 1, 2]);

        assert!((canonical.clone() + project([1, 0, 2])).expand().is_zero());
        assert!(
            (canonical + project([1, 2, 0]) + project([2, 0, 1]))
                .expand()
                .is_zero()
        );
    }

    #[test]
    fn riemann_projector_obeys_pair_exchange_and_first_bianchi() {
        let projector =
            YoungProjector::new(YoungTableau::new(vec![2, 2], vec![0, 2, 1, 3]).unwrap());
        let head = SPENSO_TAG.tensor_symbol("young_reference_riemann_identities");
        let arguments = [
            symbol!("young_riemann_a"),
            symbol!("young_riemann_b"),
            symbol!("young_riemann_c"),
            symbol!("young_riemann_d"),
        ]
        .map(Atom::var);
        let project = |selector| projector.project(head, &permute(&arguments, selector));
        let canonical = project([0, 1, 2, 3]).expand();

        assert_eq!(canonical, project([2, 3, 0, 1]).expand());
        assert!(
            (canonical + project([0, 2, 3, 1]) + project([0, 3, 1, 2]))
                .expand()
                .is_zero()
        );
    }

    #[test]
    fn factored_riemann_projector_equals_full_column_row_action() {
        let arguments = young_projector_fixtures()
            .into_iter()
            .find(|fixture| fixture.name == "riemann_2_2")
            .unwrap()
            .arguments;
        let head = FactoredRiemannProjector::tensor_symbol("young_reference_factored_equality");
        let full = YoungProjector::new(YoungTableau::new(vec![2, 2], vec![0, 2, 1, 3]).unwrap());
        let full_structural = full.actions.iter().fold(Atom::Zero, |sum, action| {
            let permuted = action
                .permutation
                .iter_slice(&arguments)
                .cloned()
                .collect::<Vec<_>>();
            let term = FactoredRiemannProjector::structural_column_carrier(head, &permuted);
            if action.sign < 0 {
                sum - term
            } else {
                sum + term
            }
        }) / Atom::num(full.hook_product());
        let reduced = FactoredRiemannProjector.project(head, &arguments);

        assert_eq!(full_structural.expand(), reduced.expand());
        assert_eq!(
            full_structural
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap(),
            reduced
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap()
        );
    }

    #[test]
    fn full_ordered_riemann_projector_preserves_current_graph_limit() {
        let fixture = young_projector_fixtures()
            .into_iter()
            .find(|fixture| fixture.name == "riemann_2_2")
            .unwrap();
        assert_current_graph_limit(&fixture.projected, fixture.name, 104, 163);
    }

    #[test]
    fn factored_riemann_ricci_projection_collapses_to_two_terms() {
        let arguments = young_projector_fixtures()
            .into_iter()
            .find(|fixture| fixture.name == "riemann_2_2")
            .unwrap()
            .arguments;
        let ricci_arguments = [
            arguments[0].clone(),
            arguments[3].clone(),
            arguments[1].clone(),
            arguments[3].clone(),
        ];
        let projected = FactoredRiemannProjector.project(
            FactoredRiemannProjector::tensor_symbol("young_reference_ricci_term_count"),
            &ricci_arguments,
        );

        assert_eq!(projected.expand().nterms(), 2);
    }

    #[test]
    fn same_head_contracted_riemann_networks_preserve_current_outcomes() {
        let fixtures = contracted_riemann_fixtures();
        assert_eq!(
            fixtures
                .iter()
                .map(|fixture| (fixture.name, fixture.factor_count))
                .collect::<Vec<_>>(),
            [
                ("riemann_two_factor_crossed", 2),
                ("riemann_three_factor_triangle", 3),
            ]
        );

        assert_successful_contraction(&fixtures[0]);
        assert_ne!(fixtures[1].expression, fixtures[1].expression.expand());
        assert_current_graph_limit(&fixtures[1].expression, fixtures[1].name, 138, 203);
        assert_current_graph_limit(&fixtures[1].renamed, fixtures[1].name, 138, 203);
    }

    #[test]
    fn distinct_head_contracted_riemann_networks_are_stable_and_nonzero() {
        let fixtures = distinct_head_contracted_riemann_fixtures();
        assert_eq!(
            fixtures
                .iter()
                .map(|fixture| (fixture.name, fixture.factor_count))
                .collect::<Vec<_>>(),
            [
                ("riemann_two_factor_crossed_distinct_heads", 2),
                ("riemann_three_factor_ricci_cycle_distinct_heads", 3),
            ]
        );

        for fixture in &fixtures {
            assert_successful_contraction(fixture);
        }
    }

    #[test]
    fn fully_projected_distinct_head_triangle_preserves_current_graph_limit() {
        let fixture = fully_projected_distinct_head_riemann_triangle_fixture();
        assert_eq!(fixture.factor_count, 3);
        assert_ne!(fixture.expression, fixture.expression.expand());
        assert_current_graph_limit(&fixture.expression, fixture.name, 138, 203);
        assert_current_graph_limit(&fixture.renamed, fixture.name, 138, 203);
    }

    #[test]
    fn mixed_shape_corpus_has_expected_term_and_hook_counts() {
        let fixtures = young_projector_fixtures();
        let counts = fixtures
            .iter()
            .map(|fixture| {
                (
                    fixture.name,
                    fixture.projector.term_count(),
                    fixture.projector.hook_product(),
                )
            })
            .collect::<Vec<_>>();

        assert_eq!(
            counts,
            [
                ("hook_2_1", 4, 3),
                ("row_heavy_3_1", 12, 8),
                ("column_heavy_2_1_1", 12, 8),
                ("riemann_2_2", 16, 12),
            ]
        );
    }

    #[test]
    fn normalized_hook_projector_is_idempotent() {
        let projector = YoungProjector::new(YoungTableau::new(vec![2, 1], vec![0, 2, 1]).unwrap());
        let head = SPENSO_TAG.tensor_symbol("young_reference_idempotence");
        let arguments = [
            symbol!("young_idempotence_a"),
            symbol!("young_idempotence_b"),
            symbol!("young_idempotence_c"),
        ]
        .map(Atom::var);
        let once = projector.project(head, &arguments);
        let twice_numerator = projector.actions.iter().fold(Atom::Zero, |sum, action| {
            let permuted = action
                .permutation
                .iter_slice(&arguments)
                .cloned()
                .collect::<Vec<_>>();
            let term = projector.project(head, &permuted);
            if action.sign < 0 {
                sum - term
            } else {
                sum + term
            }
        });
        let twice = twice_numerator / Atom::num(projector.hook_product());

        assert_eq!(twice.expand(), once.expand());
    }
}
