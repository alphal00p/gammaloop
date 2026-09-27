use spenso::shadowing::symbolica_utils::SpensoPrintSettings;
use spenso::structure::{OrderedStructure, representation::LibraryRep};
use spenso::{
    network::{
        ExecutionResult, Sequential,
        library::{DummyLibrary, function_lib::Wrap},
        parsing::{ParseSettings, ShorthandParsing, StructureInferenceMode},
    },
    structure::slot::{AbsInd, DummyAind, ParseableAind},
};

use symbolica::atom::{Atom, AtomCore, AtomView};

use crate::{
    NetworkToolingError,
    shorthands::{metric::MetricSimplifier, schoonschip::with_settings::SchoonschipWithSettings},
    tensor::{SymbolicNetParse, SymbolicTensor},
};

use super::{
    contraction::{
        ORDER_MIN_LARGEST_OPERAND_BYTES, ORDER_MIN_PRODUCT_BYTES, ORDER_MIN_PRODUCT_TERMS,
        ORDER_SMALLEST_DEGREE_MIN_LARGEST_OPERAND_BYTES, ORDER_SMALLEST_DEGREE_MIN_PRODUCT_BYTES,
        ORDER_SMALLEST_DEGREE_MIN_PRODUCT_TERMS, SchoonschipExpressionOrder,
        SchoonschipLargestDegree, SchoonschipSmallestDegree,
    },
    settings::{
        SchoonschipContractionOrder, SchoonschipMode, SchoonschipSettings, SchoonschipTraversal,
    },
    utils::TRACE_SCHOONSCHIP,
};

pub trait Schoonschip {
    fn schoonschip(&self) -> Atom;

    fn schoonschip_with_settings(&self, settings: &SchoonschipSettings) -> Atom;

    fn normalize_dots(&self) -> Atom;

    fn to_dots(&self) -> Atom;

    fn schoonschip_net<Aind: AbsInd + DummyAind + ParseableAind + 'static>(
        &self,
    ) -> Result<Atom, NetworkToolingError>;

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Result<Atom, NetworkToolingError>;

    fn schoonschip_with_net_full<Aind: AbsInd + DummyAind + ParseableAind + 'static>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        let settings = SchoonschipSettings::default().with_expanded_contracted_sums();
        self.schoonschip_with_net::<true, Aind>(&settings)
    }
}

impl Schoonschip for Atom {
    fn schoonschip(&self) -> Atom {
        self.as_view().schoonschip()
    }

    fn schoonschip_with_settings(&self, settings: &SchoonschipSettings) -> Atom {
        self.as_view().schoonschip_with_settings(settings)
    }

    fn normalize_dots(&self) -> Atom {
        self.as_view().normalize_dots()
    }

    fn to_dots(&self) -> Atom {
        self.as_view().to_dots()
    }

    fn schoonschip_net<Aind: AbsInd + DummyAind + ParseableAind + 'static>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        self.as_view().schoonschip_net::<Aind>()
    }

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Result<Atom, NetworkToolingError> {
        self.as_view()
            .schoonschip_with_net::<EXPANDSUMS, Aind>(settings)
    }
}

pub(super) struct NetworkSchoonschip<'a> {
    settings: &'a SchoonschipSettings,
}

impl NetworkSchoonschip<'_> {
    // Under opaque depth-one parsing, these index-free forms need no network
    // work. Parser-owned syntax and malformed slots must still reach the parser;
    // callers also require dot normalization to leave the expression unchanged.
    pub(super) fn scalar_requires_network(expression: AtomView<'_>) -> bool {
        use spenso::{
            network::{library::symbolic::ETS, parsing::AtomStructureExt, tags::SPENSO_TAG},
            structure::{
                abstract_index::AIND_SYMBOLS,
                representation::{LibraryRep, Representation},
                slot::{SlotMatch, SlotMatcher},
            },
        };
        let parser_heads = [
            SPENSO_TAG.bracket,
            SPENSO_TAG.pure_scalar,
            SPENSO_TAG.dot,
            SPENSO_TAG.chain,
            SPENSO_TAG.trace,
            AIND_SYMBOLS.aind,
        ];
        let mut slots = SlotMatcher::default();
        let mut required = false;
        expression.has_repeated_explicit_indices_with_observer(
            |node, slot| {
                if required {
                    return;
                }
                match slot {
                    SlotMatch::Explicit(_) => required = true,
                    SlotMatch::Opaque => {
                        required = slots.compact_representation(node).is_none()
                            || Representation::<LibraryRep>::try_from(node).is_err();
                    }
                    SlotMatch::Other => {
                        if let AtomView::Fun(function) = node {
                            let head = function.get_symbol();
                            required =
                                parser_heads.contains(&head) || head.has_tag(&SPENSO_TAG.broadcast);
                            if head == ETS.metric {
                                let mut arguments = function.iter().map(|argument| {
                                    let AtomView::Fun(vector) = argument else {
                                        return None;
                                    };
                                    if !vector.get_symbol().has_tag(&SPENSO_TAG.rank1) {
                                        return None;
                                    }
                                    let argument = slots.vector_argument(vector)?;
                                    slots.compact_representation(argument)
                                });
                                required |=
                                    match (arguments.next(), arguments.next(), arguments.next()) {
                                        (Some(Some(left)), Some(Some(right)), None) => {
                                            !left.matches(right)
                                        }
                                        _ => true,
                                    };
                            }
                        }
                    }
                }
            },
            || true,
        );
        required
    }

    fn run<const EXPANDSUMS: bool, Aind>(
        &self,
        view: AtomView<'_>,
    ) -> Result<Atom, NetworkToolingError>
    where
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    {
        let mut current = view.to_owned();
        loop {
            let next = self.apply::<EXPANDSUMS, Aind>(current.as_view())?;
            if self.settings.mode == SchoonschipMode::SinglePass || next == current {
                return Ok(next);
            }
            current = next;
        }
    }

    fn apply<const EXPANDSUMS: bool, Aind>(
        &self,
        view: AtomView<'_>,
    ) -> Result<Atom, NetworkToolingError>
    where
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    {
        if let AtomView::Add(add) = view {
            // Keep term evaluation (and callbacks) in source order, then merge
            // once instead of repeatedly normalizing the growing prefix.
            let terms = add
                .iter()
                .map(|term| self.apply::<EXPANDSUMS, Aind>(term))
                .collect::<Result<Vec<_>, _>>()?;
            return Ok(Atom::add_many(terms));
        }

        let scalar_shortcut =
            self.settings.depth_limit == Some(1) && !Self::scalar_requires_network(view);
        let normalized = view.normalize_dots();
        let scalar_shortcut = scalar_shortcut && normalized.as_view() == view;
        let new = if scalar_shortcut {
            normalized
        } else {
            let view = normalized.as_view();
            match (self.settings.expand_contracted_sums, self.settings.mode) {
                (true, SchoonschipMode::SinglePass) => {
                    self.run_once::<true, false, true, Aind>(view)
                }
                (true, SchoonschipMode::Recursive(SchoonschipTraversal::DepthFirst)) => {
                    self.run_once::<true, true, true, Aind>(view)
                }
                (true, SchoonschipMode::Recursive(SchoonschipTraversal::BreadthFirst)) => {
                    self.run_once::<true, true, false, Aind>(view)
                }
                (false, SchoonschipMode::SinglePass) => {
                    self.run_once::<EXPANDSUMS, false, true, Aind>(view)
                }
                (false, SchoonschipMode::Recursive(SchoonschipTraversal::DepthFirst)) => {
                    self.run_once::<EXPANDSUMS, true, true, Aind>(view)
                }
                (false, SchoonschipMode::Recursive(SchoonschipTraversal::BreadthFirst)) => {
                    self.run_once::<EXPANDSUMS, true, false, Aind>(view)
                }
            }?
        };

        if TRACE_SCHOONSCHIP {
            println!(
                "New: {}",
                new.printer(SpensoPrintSettings::compact().nice_symbolica())
            );
        }

        let normalized = if scalar_shortcut {
            new
        } else {
            new.normalize_dots()
        };
        // Distribute signs and numerical coefficients after each local
        // contraction so equal terms cancel before the next network pass.
        // Symbolic coefficients and products of sums remain factorized.
        Ok(if EXPANDSUMS || self.settings.expand_contracted_sums {
            normalized.expand_num()
        } else {
            normalized
        })
    }

    fn run_once<const EXPANDSUMS: bool, const RECURSE: bool, const DEPTH_FIRST: bool, Aind>(
        &self,
        view: AtomView<'_>,
    ) -> Result<Atom, NetworkToolingError>
    where
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    {
        let mut net = view
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                depth_limit: self.settings.depth_limit,
                take_first_term_from_sum: false,
                shorthand_parsing: ShorthandParsing::Opaque {
                    inference: StructureInferenceMode::Fast,
                },
                parse_composite_scalars_as_tensors: RECURSE,
                ..Default::default()
            })
            .map_err(|error| NetworkToolingError::Parse {
                reason: error.to_string(),
            })?;
        let lib = DummyLibrary::<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>::new();

        let execution = match self.settings.contraction_order {
            SchoonschipContractionOrder::SmallestDegree => net
                .execute::<
                    Sequential,
                    SchoonschipSmallestDegree<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {}),
            SchoonschipContractionOrder::LargestDegree => net
                .execute::<
                    Sequential,
                    SchoonschipLargestDegree<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {}),
            SchoonschipContractionOrder::MinLargestOperandBytes => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_MIN_LARGEST_OPERAND_BYTES,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {}),
            SchoonschipContractionOrder::MinProductTerms => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_MIN_PRODUCT_TERMS,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {}),
            SchoonschipContractionOrder::MinProductBytes => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_MIN_PRODUCT_BYTES,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {}),
            SchoonschipContractionOrder::SmallestDegreeMinLargestOperandBytes => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_SMALLEST_DEGREE_MIN_LARGEST_OPERAND_BYTES,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {}),
            SchoonschipContractionOrder::SmallestDegreeMinProductTerms => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_SMALLEST_DEGREE_MIN_PRODUCT_TERMS,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {}),
            SchoonschipContractionOrder::SmallestDegreeMinProductBytes => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_SMALLEST_DEGREE_MIN_PRODUCT_BYTES,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {}),
        };
        execution.map_err(|error| NetworkToolingError::Execute {
            reason: error.to_string(),
        })?;

        Ok(
            match net
                .result_tensor(&lib)
                .map_err(|error| NetworkToolingError::Result {
                    reason: error.to_string(),
                })? {
                ExecutionResult::One => Atom::num(1),
                ExecutionResult::Zero => Atom::Zero,
                ExecutionResult::Val(tensor) => tensor.expression.clone(),
            },
        )
    }
}

impl Schoonschip for AtomView<'_> {
    fn normalize_dots(&self) -> Atom {
        super::normalize_dots::DotNormalizer::run(*self)
    }

    fn schoonschip(&self) -> Atom {
        self.schoonschip_with_settings(&SchoonschipSettings::default())
    }

    fn schoonschip_with_settings(&self, settings: &SchoonschipSettings) -> Atom {
        SchoonschipWithSettings { settings }.run::<false>(*self)
    }

    fn to_dots(&self) -> Atom {
        // Resolve explicit vector pairs before tensor collection can absorb a
        // registered momentum into an ordinary, untagged vector head.
        let explicit = crate::shorthands::metric::to_dots_impl(*self);
        let simplified = explicit
            .schoonschip_with_settings(&SchoonschipSettings::default().with_rank1_tensors())
            .metric_shorthand_to_dot();
        // Explicit representation slots also identify vectors whose heads
        // were created as ordinary Symbolica symbols without rank-one tags.
        crate::shorthands::metric::to_dots_impl(simplified.as_view())
    }

    fn schoonschip_net<Aind: AbsInd + DummyAind + ParseableAind + 'static>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        self.schoonschip_with_net::<false, Aind>(&SchoonschipSettings::default_network())
    }

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Result<Atom, NetworkToolingError> {
        NetworkSchoonschip { settings }.run::<EXPANDSUMS, Aind>(*self)
    }
}
