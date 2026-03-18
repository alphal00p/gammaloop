use clap::Arg;
use clap::builder::ArgExt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SelectorKind {
    Any,
    Amplitude,
    CrossSection,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ArgValueCompletion {
    ProcessSelector(SelectorKind),
    IntegrandSelector(SelectorKind),
    SelectedIntegrandTarget,
    Disabled,
}

impl ArgExt for ArgValueCompletion {}

pub(crate) trait CompletionArgExt {
    fn completion_process_selector(self, kind: SelectorKind) -> Self;
    fn completion_integrand_selector(self, kind: SelectorKind) -> Self;
    fn completion_selected_integrand_target(self) -> Self;
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

    fn completion_disable_special_value(self) -> Self {
        self.add(ArgValueCompletion::Disabled)
    }
}

pub(crate) fn arg_value_completion(arg: &Arg) -> Option<ArgValueCompletion> {
    arg.get::<ArgValueCompletion>().copied()
}
