use gammalooprs::model::ParameterType;

use crate::integrand_info::IntegrandKind;

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

#[cfg(feature = "cli")]
pub(crate) fn arg_value_completion(arg: &Arg) -> Option<ArgValueCompletion> {
    arg.get::<ArgValueCompletion>().copied()
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModelParameterCompletionEntry {
    pub name: String,
    pub parameter_type: ParameterType,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModelVertexCompletionEntry {
    pub name: String,
    pub particles: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProcessKind {
    Amplitude,
    CrossSection,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProcessCompletionEntry {
    pub id: usize,
    pub name: String,
    pub kind: ProcessKind,
    pub integrand_names: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IntegrandDetailCompletionEntry {
    pub process_name: String,
    pub integrand_name: String,
    pub kind: IntegrandKind,
    pub master_graph_names: Vec<String>,
    pub raised_all_signatures: Vec<String>,
    pub raised_massive_signatures: Vec<String>,
    pub raised_massless_signatures: Vec<String>,
    pub cycle_signatures: Vec<String>,
    pub amplitude_raised_all_signatures: Vec<String>,
    pub amplitude_raised_massive_signatures: Vec<String>,
    pub amplitude_raised_massless_signatures: Vec<String>,
    pub amplitude_cycle_signatures: Vec<String>,
    pub raised_cut_all_signatures: Vec<String>,
    pub raised_cut_massive_signatures: Vec<String>,
    pub raised_cut_massless_signatures: Vec<String>,
    pub categories: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IrProfileCompletionEntry {
    pub process_name: String,
    pub integrand_name: String,
    pub graph_names: Vec<String>,
    pub graph_limit_entries: Vec<String>,
}
