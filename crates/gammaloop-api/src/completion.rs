use clap::builder::ArgExt;
use clap::Arg;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SelectorKind {
    Any,
    Amplitude,
    CrossSection,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SelectRaisedSignatureScope {
    All,
    Massive,
    Massless,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ArgValueCompletion {
    ProcessSelector(SelectorKind),
    IntegrandSelector(SelectorKind),
    SelectedIntegrandTarget,
    SelectedMasterGraph,
    SelectedIntegrandCategory,
    SelectRaisedSignature(SelectRaisedSignatureScope),
    SelectRaisedCutSignature(SelectRaisedSignatureScope),
    SelectCycleSignature,
    SelectVertexSignature,
    SelectParticleSignature,
    Disabled,
}

impl ArgExt for ArgValueCompletion {}

pub(crate) trait CompletionArgExt {
    fn completion_process_selector(self, kind: SelectorKind) -> Self;
    fn completion_integrand_selector(self, kind: SelectorKind) -> Self;
    fn completion_selected_integrand_target(self) -> Self;
    fn completion_selected_master_graph(self) -> Self;
    fn completion_selected_integrand_category(self) -> Self;
    fn completion_select_raised_signature(self, scope: SelectRaisedSignatureScope) -> Self;
    fn completion_select_raised_cut_signature(self, scope: SelectRaisedSignatureScope) -> Self;
    fn completion_select_cycle_signature(self) -> Self;
    fn completion_select_vertex_signature(self) -> Self;
    fn completion_select_particle_signature(self) -> Self;
    fn completion_disable_special_value(self) -> Self;
}

impl CompletionArgExt for Arg {
    fn completion_process_selector(self, kind: SelectorKind) -> Self {
        self.add(ArgValueCompletion::ProcessSelector(kind))
    }

    fn completion_integrand_selector(self, kind: SelectorKind) -> Self {
        self.add(ArgValueCompletion::IntegrandSelector(kind))
    }

    fn completion_selected_integrand_target(self) -> Self {
        self.add(ArgValueCompletion::SelectedIntegrandTarget)
    }

    fn completion_selected_master_graph(self) -> Self {
        self.add(ArgValueCompletion::SelectedMasterGraph)
    }

    fn completion_selected_integrand_category(self) -> Self {
        self.add(ArgValueCompletion::SelectedIntegrandCategory)
    }

    fn completion_select_raised_signature(self, scope: SelectRaisedSignatureScope) -> Self {
        self.add(ArgValueCompletion::SelectRaisedSignature(scope))
    }

    fn completion_select_raised_cut_signature(self, scope: SelectRaisedSignatureScope) -> Self {
        self.add(ArgValueCompletion::SelectRaisedCutSignature(scope))
    }

    fn completion_select_cycle_signature(self) -> Self {
        self.add(ArgValueCompletion::SelectCycleSignature)
    }

    fn completion_select_vertex_signature(self) -> Self {
        self.add(ArgValueCompletion::SelectVertexSignature)
    }

    fn completion_select_particle_signature(self) -> Self {
        self.add(ArgValueCompletion::SelectParticleSignature)
    }

    fn completion_disable_special_value(self) -> Self {
        self.add(ArgValueCompletion::Disabled)
    }
}

pub(crate) fn arg_value_completion(arg: &Arg) -> Option<ArgValueCompletion> {
    arg.get::<ArgValueCompletion>().copied()
}
