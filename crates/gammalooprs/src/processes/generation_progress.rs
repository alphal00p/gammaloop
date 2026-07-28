use std::{
    cell::RefCell,
    sync::{
        Arc, Mutex, OnceLock,
        atomic::{AtomicU8, Ordering},
    },
    time::Duration,
};

use super::GraphGenerationStats;

thread_local! {
    static GENERATION_PROGRESS_CONTEXT: RefCell<Vec<String>> = const { RefCell::new(Vec::new()) };
}

pub struct GenerationProgressContextGuard;

impl Drop for GenerationProgressContextGuard {
    fn drop(&mut self) {
        GENERATION_PROGRESS_CONTEXT.with(|context| {
            context.borrow_mut().pop();
        });
    }
}

pub fn enter_progress_context(context: impl Into<String>) -> GenerationProgressContextGuard {
    GENERATION_PROGRESS_CONTEXT.with(|stack| stack.borrow_mut().push(context.into()));
    GenerationProgressContextGuard
}

pub fn detailed_progress_message(message: &str) -> String {
    GENERATION_PROGRESS_CONTEXT.with(|context| {
        let context = context.borrow();
        if context.is_empty() {
            message.to_string()
        } else {
            format!("{} / {message}", context.join(" / "))
        }
    })
}

pub fn enter_detailed_progress_span(message: &str) -> Option<tracing::span::EnteredSpan> {
    if !detailed_progress_enabled() {
        return None;
    }

    let progress_message = detailed_progress_message(message);
    Some(
        tracing::info_span!(
            "Generation progress",
            indicatif.pb_show = true,
            indicatif.pb_msg = %progress_message,
        )
        .entered(),
    )
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenerationProgressMode {
    Detailed,
    Aggregate,
}

impl GenerationProgressMode {
    fn as_u8(self) -> u8 {
        match self {
            Self::Detailed => 0,
            Self::Aggregate => 1,
        }
    }

    fn from_u8(value: u8) -> Self {
        match value {
            1 => Self::Aggregate,
            _ => Self::Detailed,
        }
    }
}

static GENERATION_PROGRESS_MODE: AtomicU8 = AtomicU8::new(0);

pub struct GenerationProgressModeGuard {
    previous: u8,
}

impl GenerationProgressModeGuard {
    pub fn set(mode: GenerationProgressMode) -> Self {
        let previous = GENERATION_PROGRESS_MODE.swap(mode.as_u8(), Ordering::Relaxed);
        Self { previous }
    }
}

impl Drop for GenerationProgressModeGuard {
    fn drop(&mut self) {
        GENERATION_PROGRESS_MODE.store(self.previous, Ordering::Relaxed);
    }
}

pub fn current_generation_progress_mode() -> GenerationProgressMode {
    GenerationProgressMode::from_u8(GENERATION_PROGRESS_MODE.load(Ordering::Relaxed))
}

pub fn detailed_progress_enabled() -> bool {
    current_generation_progress_mode() == GenerationProgressMode::Detailed
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenerationProcessKind {
    Amplitude,
    CrossSection,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenerationProgressPhase {
    CutDiscovery,
    GraphGeneration,
    Backend,
}

pub trait GenerationProgressObserver: Send + Sync {
    fn begin_phase(
        &self,
        _phase: GenerationProgressPhase,
        _kind: GenerationProcessKind,
        _process: &str,
        _integrand: &str,
        _total_graphs: usize,
        _total_cuts: Option<usize>,
    ) {
    }

    fn graph_started(
        &self,
        _kind: GenerationProcessKind,
        _integrand: &str,
        _graph: &str,
        _cut_count: Option<usize>,
    ) {
    }

    fn graph_finished(
        &self,
        _kind: GenerationProcessKind,
        _integrand: &str,
        _graph: &str,
        _stats: &GraphGenerationStats,
        _completed_cuts: Option<usize>,
    ) {
    }

    fn cuts_discovered(
        &self,
        _integrand: &str,
        _graph: &str,
        _st_cut_count: usize,
        _valid_cut_count: usize,
    ) {
    }

    fn cut_finished(&self, _integrand: &str, _graph: &str, _cut_count: usize) {}

    fn backend_started(&self, _kind: GenerationProcessKind, _integrand: &str, _graph_count: usize) {
    }

    fn backend_finished(&self, _kind: GenerationProcessKind, _integrand: &str, _elapsed: Duration) {
    }
}

static GENERATION_PROGRESS_OBSERVER: OnceLock<Mutex<Option<Arc<dyn GenerationProgressObserver>>>> =
    OnceLock::new();

fn observer_slot() -> &'static Mutex<Option<Arc<dyn GenerationProgressObserver>>> {
    GENERATION_PROGRESS_OBSERVER.get_or_init(|| Mutex::new(None))
}

fn observer() -> Option<Arc<dyn GenerationProgressObserver>> {
    observer_slot().lock().ok().and_then(|guard| guard.clone())
}

pub struct GenerationProgressObserverGuard {
    previous: Option<Arc<dyn GenerationProgressObserver>>,
}

impl GenerationProgressObserverGuard {
    pub fn set(observer: Arc<dyn GenerationProgressObserver>) -> Self {
        let previous = observer_slot()
            .lock()
            .expect("generation progress observer mutex is poisoned")
            .replace(observer);
        Self { previous }
    }
}

impl Drop for GenerationProgressObserverGuard {
    fn drop(&mut self) {
        *observer_slot()
            .lock()
            .expect("generation progress observer mutex is poisoned") = self.previous.take();
    }
}

pub fn begin_phase(
    phase: GenerationProgressPhase,
    kind: GenerationProcessKind,
    process: &str,
    integrand: &str,
    total_graphs: usize,
    total_cuts: Option<usize>,
) {
    if let Some(observer) = observer() {
        observer.begin_phase(phase, kind, process, integrand, total_graphs, total_cuts);
    }
}

pub fn graph_started(
    kind: GenerationProcessKind,
    integrand: &str,
    graph: &str,
    cut_count: Option<usize>,
) {
    if let Some(observer) = observer() {
        observer.graph_started(kind, integrand, graph, cut_count);
    }
}

pub fn graph_finished(
    kind: GenerationProcessKind,
    integrand: &str,
    graph: &str,
    stats: &GraphGenerationStats,
    completed_cuts: Option<usize>,
) {
    if let Some(observer) = observer() {
        observer.graph_finished(kind, integrand, graph, stats, completed_cuts);
    }
}

pub fn cuts_discovered(integrand: &str, graph: &str, st_cut_count: usize, valid_cut_count: usize) {
    if let Some(observer) = observer() {
        observer.cuts_discovered(integrand, graph, st_cut_count, valid_cut_count);
    }
}

pub fn cut_finished(integrand: &str, graph: &str, cut_count: usize) {
    if let Some(observer) = observer() {
        observer.cut_finished(integrand, graph, cut_count);
    }
}

pub fn backend_started(kind: GenerationProcessKind, integrand: &str, graph_count: usize) {
    if let Some(observer) = observer() {
        observer.backend_started(kind, integrand, graph_count);
    }
}

pub fn backend_finished(kind: GenerationProcessKind, integrand: &str, elapsed: Duration) {
    if let Some(observer) = observer() {
        observer.backend_finished(kind, integrand, elapsed);
    }
}
