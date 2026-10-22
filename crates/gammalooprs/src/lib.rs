#![cfg_attr(feature = "fail-on-warnings", deny(warnings))]
// #![deny(clippy::all)]
// #![warn(clippy::pedantic)]
#![warn(clippy::all)]
// #![warn(clippy::restriction)]
// #![warn(clippy::nursery)]
// #![warn(clippy::cargo)]
// #![feature(min_specialization)]
//

pub mod cff;
pub mod feyngen;
pub mod graph;
pub mod initialisation;
pub mod integrands;
pub mod integrate;
pub mod model;
pub mod momentum;
pub mod numerator;
pub mod observables;
pub mod processes;
pub mod settings;
pub mod subtraction;
pub mod utils;
pub mod uv;
pub use gammaloop_tracing_filter::LogMessage;
use integrands::*;
use model::Model;
use momentum::signature::ExternalSignature;
#[cfg(test)]
use std::env;
#[cfg(test)]
use std::path::{Path, PathBuf};
#[cfg(test)]
use std::sync::OnceLock;
use std::sync::atomic::AtomicBool;
use symbolica::state::HasStateMap;
use symbolica::state::StateMap;
use utils::F;
use utils::FloatLike;

pub static INTERRUPTED: AtomicBool = AtomicBool::new(false);
pub static ITERATION_ABORT_REQUESTED: AtomicBool = AtomicBool::new(false);

pub const GAMMALOOP_NAMESPACE: &str = "GL";
pub const MAX_CORES: usize = 1000;

#[cfg(test)]
pub(crate) fn test_workspace_root() -> &'static Path {
    static WORKSPACE_ROOT: OnceLock<PathBuf> = OnceLock::new();

    WORKSPACE_ROOT
        .get_or_init(|| {
            let current_dir = env::current_dir().ok();
            let current_exe = env::current_exe().ok();
            // Cargo and Nextest remap this at runtime; direct binaries retain the native fallback.
            let manifest_dir = env::var_os("CARGO_MANIFEST_DIR")
                .map(PathBuf::from)
                .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")));

            let mut candidates = Vec::new();
            if let Some(dir) = current_dir.as_deref() {
                candidates.extend(dir.ancestors().map(Path::to_path_buf));
            }
            if let Some(exe) = current_exe.as_deref() {
                candidates.extend(exe.ancestors().skip(1).map(Path::to_path_buf));
            }
            candidates.extend(manifest_dir.ancestors().map(Path::to_path_buf));

            candidates
                .into_iter()
                .find(|path| {
                    path.join("Cargo.toml").is_file()
                        && path.join("tests/resources").is_dir()
                        && path.join("examples/cli").is_dir()
                })
                .unwrap_or_else(|| {
                    panic!(
                        "Failed to locate workspace root from current dir '{:?}', current exe '{:?}', or manifest dir '{}'",
                        current_dir,
                        current_exe,
                        manifest_dir.display()
                    )
                })
        })
        .as_path()
}

pub trait GammaLoopContext: HasStateMap + HasModel {}

pub trait HasModel {
    /// Returns the physics model associated with this computation context.
    fn get_model(&self) -> &Model;
}

impl HasModel for Model {
    fn get_model(&self) -> &Model {
        self
    }
}

impl HasModel for GammaLoopContextContainer<'_> {
    fn get_model(&self) -> &Model {
        self.model
    }
}

#[derive(Clone, Copy)]
pub struct GammaLoopContextContainer<'a> {
    /// Symbolica state map used to resolve symbols during GammaLoop computations.
    pub state_map: &'a StateMap,
    /// Physics model that supplies particles, parameters, and interactions.
    pub model: &'a Model,
}

impl<'a> HasStateMap for GammaLoopContextContainer<'a> {
    fn get_state_map(&self) -> &StateMap {
        self.state_map
    }
}

impl<'a> GammaLoopContext for GammaLoopContextContainer<'a> {}

#[cfg(not(feature = "higher_loops"))]
pub const MAX_LOOP: usize = 3;
#[cfg(feature = "higher_loops")]
pub const MAX_LOOP: usize = 6;

pub fn set_interrupt_handler() {
    INTERRUPTED.store(false, std::sync::atomic::Ordering::Relaxed);
    ITERATION_ABORT_REQUESTED.store(false, std::sync::atomic::Ordering::Relaxed);
    let _ = ctrlc::set_handler(|| {
        INTERRUPTED.store(true, std::sync::atomic::Ordering::Relaxed);
    });
}

#[inline]
pub fn request_interrupt() {
    INTERRUPTED.store(true, std::sync::atomic::Ordering::Relaxed);
}

#[inline]
pub fn request_iteration_abort() {
    ITERATION_ABORT_REQUESTED.store(true, std::sync::atomic::Ordering::Relaxed);
}

#[inline]
pub fn is_interrupt_requested() -> bool {
    INTERRUPTED.load(std::sync::atomic::Ordering::Relaxed)
}

#[inline]
pub(crate) fn is_interrupted() -> bool {
    is_interrupt_requested()
}

#[inline]
pub(crate) fn is_iteration_abort_requested() -> bool {
    ITERATION_ABORT_REQUESTED.load(std::sync::atomic::Ordering::Relaxed)
}

#[inline]
pub(crate) fn set_interrupted(flag: bool) {
    INTERRUPTED.store(flag, std::sync::atomic::Ordering::Relaxed);
}

#[inline]
pub fn clear_interrupt_request() {
    set_interrupted(false);
}

#[inline]
pub(crate) fn clear_iteration_abort_request() {
    ITERATION_ABORT_REQUESTED.store(false, std::sync::atomic::Ordering::Relaxed);
}

#[derive(Clone, Copy, Debug)]
pub enum DependentMomentaConstructor<'a> {
    Amplitude(&'a ExternalSignature),
    //CrossSection {
    //    external_connections: &'a [ExternalConnection],
    //}, // at the moment I assume the first n/2 externals are incoming and the second n/2 are outgoing, the mapping is (0, n/2), (1, n/2+1), (2, n/2+2), ...
    CrossSection,
}
