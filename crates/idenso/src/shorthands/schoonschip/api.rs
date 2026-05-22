use spenso::{
    network::{
        ExecutionResult, Sequential,
        library::{DummyLibrary, function_lib::Wrap},
        parsing::{ParseSettings, ShorthandParsing, StructureInferenceMode},
    },
    shadowing::symbolica_utils::SpensoPrintSettings,
    structure::slot::{AbsInd, DummyAind, ParseableAind},
};

use symbolica::atom::{Atom, AtomCore, AtomView};

use crate::tensor::{SymbolicNetParse, SymbolicTensor};

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

    fn schoonschip_net<Aind: AbsInd + DummyAind + ParseableAind + 'static>(&self) -> Atom;

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom;
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

    fn schoonschip_net<Aind: AbsInd + DummyAind + ParseableAind + 'static>(&self) -> Atom {
        self.as_view().schoonschip_net::<Aind>()
    }

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom {
        self.as_view()
            .schoonschip_with_net::<EXPANDSUMS, Aind>(settings)
    }
}

struct NetworkSchoonschip<'a> {
    settings: &'a SchoonschipSettings,
}

impl NetworkSchoonschip<'_> {
    fn run<const EXPANDSUMS: bool, Aind>(&self, view: AtomView<'_>) -> Atom
    where
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    {
        let mut current = view.to_owned();
        loop {
            let next = self.apply::<EXPANDSUMS, Aind>(current.as_view());
            if self.settings.mode == SchoonschipMode::SinglePass || next == current {
                return next;
            }
            current = next;
        }
    }

    fn apply<const EXPANDSUMS: bool, Aind>(&self, view: AtomView<'_>) -> Atom
    where
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    {
        if let AtomView::Add(add) = view {
            return add
                .iter()
                .map(|term| self.apply::<EXPANDSUMS, Aind>(term))
                .fold(Atom::Zero, |sum, term| sum + term);
        }

        let new = match (self.settings.expand_contracted_sums, self.settings.mode) {
            (true, SchoonschipMode::SinglePass) => self.run_once::<true, false, true, Aind>(view),
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
        };

        if TRACE_SCHOONSCHIP {
            println!(
                "New: {}",
                new.printer(SpensoPrintSettings::compact().nice_symbolica())
            );
        }

        new.normalize_dots()
    }

    fn run_once<const EXPANDSUMS: bool, const RECURSE: bool, const DEPTH_FIRST: bool, Aind>(
        &self,
        view: AtomView<'_>,
    ) -> Atom
    where
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    {
        let normalized = view.normalize_dots();
        let mut net = normalized
            .as_view()
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                depth_limit: self.settings.depth_limit,
                take_first_term_from_sum: false,
                shorthand_parsing: ShorthandParsing::Opaque {
                    inference: StructureInferenceMode::Fast,
                },
                parse_composite_scalars_as_tensors: RECURSE,
                ..Default::default()
            })
            .unwrap();
        let lib = DummyLibrary::<SymbolicTensor<Aind>>::new();

        match self.settings.contraction_order {
            SchoonschipContractionOrder::SmallestDegree => net
                .execute::<
                    Sequential,
                    SchoonschipSmallestDegree<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
            SchoonschipContractionOrder::LargestDegree => net
                .execute::<
                    Sequential,
                    SchoonschipLargestDegree<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
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
                >(&lib, &Wrap {})
                .unwrap(),
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
                >(&lib, &Wrap {})
                .unwrap(),
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
                >(&lib, &Wrap {})
                .unwrap(),
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
                >(&lib, &Wrap {})
                .unwrap(),
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
                >(&lib, &Wrap {})
                .unwrap(),
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
                >(&lib, &Wrap {})
                .unwrap(),
        };

        match net.result_tensor(&lib).unwrap() {
            ExecutionResult::One => Atom::num(1),
            ExecutionResult::Zero => Atom::Zero,
            ExecutionResult::Val(tensor) => tensor.expression.clone(),
        }
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
        super::with_settings::schoonschip_with_settings(*self, settings)
    }

    fn to_dots(&self) -> Atom {
        self.schoonschip_with_settings(&SchoonschipSettings::default().with_rank1_tensors())
    }

    fn schoonschip_net<Aind: AbsInd + DummyAind + ParseableAind + 'static>(&self) -> Atom {
        self.schoonschip_with_net::<false, Aind>(&SchoonschipSettings::default_network())
    }

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom {
        NetworkSchoonschip { settings }.run::<EXPANDSUMS, Aind>(*self)
    }
}
