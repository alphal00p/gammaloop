pub struct SchoonschipSettings {
    pub(super) depth_limit: Option<usize>,
    pub(super) mode: SchoonschipMode,
    pub(super) parse_inner_products: bool,
    pub(super) expand_contracted_sums: bool,
    pub(super) simplify_chain_like_functions: bool,
    pub(super) contraction_order: SchoonschipContractionOrder,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum SchoonschipMode {
    SinglePass,
    Recursive(SchoonschipTraversal),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SchoonschipTraversal {
    DepthFirst,
    BreadthFirst,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum SchoonschipContractionOrder {
    #[default]
    SmallestDegree,
    LargestDegree,
    MinLargestOperandBytes,
    MinProductTerms,
    MinProductBytes,
    SmallestDegreeMinLargestOperandBytes,
    SmallestDegreeMinProductTerms,
    SmallestDegreeMinProductBytes,
}

impl Default for SchoonschipSettings {
    fn default() -> Self {
        Self::partial()
    }
}

impl SchoonschipSettings {
    pub fn new(depth_limit: Option<usize>) -> Self {
        Self::breadth_first(depth_limit)
    }

    pub fn depth_first(depth_limit: Option<usize>) -> Self {
        Self {
            depth_limit,
            mode: SchoonschipMode::Recursive(SchoonschipTraversal::DepthFirst),
            parse_inner_products: true,
            expand_contracted_sums: false,
            simplify_chain_like_functions: false,
            contraction_order: SchoonschipContractionOrder::default(),
        }
    }

    pub fn breadth_first(depth_limit: Option<usize>) -> Self {
        Self {
            depth_limit,
            mode: SchoonschipMode::Recursive(SchoonschipTraversal::BreadthFirst),
            parse_inner_products: true,
            expand_contracted_sums: false,
            simplify_chain_like_functions: false,
            contraction_order: SchoonschipContractionOrder::default(),
        }
    }

    pub fn single_pass(depth_limit: Option<usize>) -> Self {
        Self {
            depth_limit,
            mode: SchoonschipMode::SinglePass,
            parse_inner_products: true,
            expand_contracted_sums: false,
            simplify_chain_like_functions: false,
            contraction_order: SchoonschipContractionOrder::default(),
        }
    }

    pub fn partial() -> Self {
        Self::new(Some(1))
    }

    pub fn with_depth(depth_limit: usize) -> Self {
        Self::new(Some(depth_limit))
    }

    pub fn breadth_first_with_depth(depth_limit: usize) -> Self {
        Self::breadth_first(Some(depth_limit))
    }

    pub fn full() -> Self {
        Self::single_pass(None)
    }

    pub fn with_expanded_contracted_sums(mut self) -> Self {
        self.expand_contracted_sums = true;
        self
    }

    pub fn with_chain_like_functions(mut self) -> Self {
        self.simplify_chain_like_functions = true;
        self
    }

    pub fn without_chain_like_functions(mut self) -> Self {
        self.simplify_chain_like_functions = false;
        self
    }

    pub fn with_contraction_order(mut self, order: SchoonschipContractionOrder) -> Self {
        self.contraction_order = order;
        self
    }

    pub(super) fn without_parse_inner_products(mut self) -> Self {
        self.parse_inner_products = false;
        self
    }
}
