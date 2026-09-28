/// Options for the local ordered metric/vector rewrite fallback.
/// Factored network contraction and materialization are owned by the typed engine.
pub(crate) struct SchoonschipSettings {
    pub simplify_chain_like_functions: bool,
    pub schoonschip_rank1_tensors: bool,
}

impl Default for SchoonschipSettings {
    fn default() -> Self {
        Self {
            simplify_chain_like_functions: false,
            schoonschip_rank1_tensors: true,
        }
    }
}

impl SchoonschipSettings {
    pub(crate) fn with_chain_like_functions(mut self) -> Self {
        self.simplify_chain_like_functions = true;
        self
    }

    #[cfg(test)]
    pub(crate) fn with_rank1_tensors(mut self) -> Self {
        self.schoonschip_rank1_tensors = true;
        self
    }

    pub(crate) fn without_rank1_tensors(mut self) -> Self {
        self.schoonschip_rank1_tensors = false;
        self
    }
}
