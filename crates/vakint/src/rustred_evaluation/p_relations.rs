//! The persistent relation cache of the RustRed p library (vk5l WP4, user
//! requirement 2026-11-14): every key a boxed node solves with a recipe at run
//! time is written to a cache that later runs and other processes reuse.
//!
//! Format and lookup are RustRed's (`rustred.relation-entry.v1`,
//! `rustred::propagator::relation_entries`): at a boxed node a key is taken by
//! the table's rules, then the shipped (read-only) relations, then this
//! persistent cache, and only then solved by a recipe, whose relation is
//! inserted. vakint chooses the locations ([`PRelationCacheSettings`]) and
//! reports the counters ([`crate::Vakint::p_relation_cache_statistics`]).
//!
//! Durability: RustRed writes each inserted entry as one line of an
//! append-only segment file owned by the writing process and syncs it at
//! once, so there is nothing to flush at the end of an evaluation or
//! periodically in a long run, and a crash loses at most the line being
//! written (a torn last line is ignored on load). Several processes may share
//! one root.
//!
//! The roots are process-wide in RustRed and bind when a boxed node's recipe
//! dispatch is built (on the node's first use in the process): vakint
//! installs them before it loads or uses a p library, so the settings of the
//! first p evaluation of a process are the ones in force.

use std::path::{Path, PathBuf};
use std::sync::Mutex;

use rustred::propagator::recipes::{
    RecipeCacheRoots, RecipeStatistics, recipe_cache_roots, recipe_statistics,
    set_recipe_cache_roots,
};
use rustred::propagator::relation_entries::{RelationCacheStatistics, relation_cache_statistics};

use crate::VakintSettings;

/// Directory name of the relation cache inside vakint's cache directory.
pub const P_RELATION_CACHE_DIRECTORY: &str = "p_relations";

/// Environment variable naming vakint's cache directory when the settings
/// name none ([`PRelationCacheLocation::Default`]).
pub const VAKINT_CACHE_DIR_VARIABLE: &str = "VAKINT_CACHE_DIR";

/// Where the persistent relation cache lives.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum PRelationCacheLocation {
    /// `<vakint cache directory>/p_relations`: the Ordered step cache's
    /// directory ([`crate::OrderedCacheSettings::directory`]) when one is set,
    /// else `$VAKINT_CACHE_DIR`, else `$XDG_CACHE_HOME/vakint`, else
    /// `$HOME/.cache/vakint`; no cache when none of these exists.
    #[default]
    Default,
    /// This directory (created on first insert).
    Directory(PathBuf),
    /// No persistent cache: solved relations live in the process only.
    Disabled,
}

/// The relation caches of the RustRed p library: where the three kinds of a
/// node's service (plan ARTIFACT_THREE_KINDS) other than its table rules
/// live: EXPLICIT key expressions (stored relations, read-only: the package's
/// `relations/` directory and `shipped` roots), the persistent cache of keys
/// solved by a RECIPE at run time, and their limits.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PRelationCacheSettings {
    /// The read-write persistent cache.
    pub location: PRelationCacheLocation,
    /// Read-only relation roots (shipped or pre-filled relations, e.g. the
    /// optional download), searched before the persistent cache.
    pub shipped: Vec<PathBuf>,
    /// Use the explicit key expressions the loaded package ships in its
    /// `relations/` directory (default); `false` ignores them.
    pub package_relations: bool,
    /// Size cap of the persistent cache in bytes: beyond it no new entry is
    /// written (counted), nothing is evicted; `None`: no cap.
    pub max_bytes: Option<u64>,
    /// Whether relation entries without per-key evidence are used; `None`:
    /// RustRed's default (allowed and reported).
    pub allow_field_entries: Option<bool>,
}

impl Default for PRelationCacheSettings {
    fn default() -> Self {
        Self {
            location: PRelationCacheLocation::Default,
            shipped: Vec::new(),
            package_relations: true,
            max_bytes: None,
            allow_field_entries: None,
        }
    }
}

/// Directory of a package's explicit key expressions (one directory per node
/// with its `MANIFEST.json` and entries, `rustred.node-service-manifest.v1`).
pub const P_PACKAGE_RELATIONS_DIRECTORY: &str = "relations";

impl PRelationCacheSettings {
    /// The persistent cache directory these settings resolve to.
    pub fn resolved_directory(&self, ordered_cache_directory: Option<&Path>) -> Option<PathBuf> {
        match &self.location {
            PRelationCacheLocation::Disabled => None,
            PRelationCacheLocation::Directory(directory) => Some(directory.clone()),
            PRelationCacheLocation::Default => default_cache_directory(ordered_cache_directory)
                .map(|directory| directory.join(P_RELATION_CACHE_DIRECTORY)),
        }
    }

    fn roots(
        &self,
        ordered_cache_directory: Option<&Path>,
        package_relations: Option<&Path>,
    ) -> RecipeCacheRoots {
        let mut prefilled = self.shipped.clone();
        if self.package_relations
            && let Some(root) = package_relations
        {
            prefilled.insert(0, root.to_path_buf());
        }
        RecipeCacheRoots {
            user: self.resolved_directory(ordered_cache_directory),
            prefilled,
            max_bytes: self.max_bytes,
            allow_field: self.allow_field_entries,
        }
    }
}

/// vakint's cache directory (module documentation of [`PRelationCacheLocation::Default`]).
fn default_cache_directory(ordered_cache_directory: Option<&Path>) -> Option<PathBuf> {
    let variable = |name: &str| std::env::var_os(name).filter(|value| !value.is_empty());
    ordered_cache_directory
        .map(Path::to_path_buf)
        .or_else(|| variable(VAKINT_CACHE_DIR_VARIABLE).map(PathBuf::from))
        .or_else(|| variable("XDG_CACHE_HOME").map(|root| PathBuf::from(root).join("vakint")))
        .or_else(|| variable("HOME").map(|home| PathBuf::from(home).join(".cache").join("vakint")))
}

/// The roots vakint installed last (`None`: never installed in this process).
static INSTALLED: Mutex<Option<RecipeCacheRoots>> = Mutex::new(None);

/// Install the relation cache roots of `settings` in RustRed (before a p
/// library is loaded or used). A later call with other roots changes them for
/// the boxed nodes not used yet in this process only (module documentation);
/// such a change is logged.
pub(crate) fn install(settings: &VakintSettings, package_relations: Option<&Path>) {
    rustred::propagator::recipes::set_outside_settings(settings.rustred_p_outside);
    let roots = settings.rustred_p_relation_cache.roots(
        settings.ordered_cache.directory.as_deref(),
        package_relations,
    );
    let mut installed = INSTALLED
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    if installed.as_ref() == Some(&roots) {
        return;
    }
    if let Some(previous) = installed.as_ref() {
        log::info!(
            "RustRed p relation cache: roots change from {previous:?} to {roots:?}; boxed nodes \
             already used in this process keep their roots"
        );
    } else {
        log::debug!("RustRed p relation cache: roots {roots:?}");
    }
    set_recipe_cache_roots(roots.clone());
    *installed = Some(roots);
}

/// Counters and locations of the relation caches of this process
/// ([`crate::Vakint::p_relation_cache_statistics`]).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PRelationCacheStatistics {
    /// The persistent cache directory in force (`None`: disabled).
    pub directory: Option<PathBuf>,
    /// The read-only roots in force.
    pub shipped: Vec<PathBuf>,
    /// RustRed's cache counters since the process started: shipped and
    /// persistent-cache hits, misses, inserts (and those skipped by a size
    /// cap), foreign and corrupt entries, bytes written.
    pub cache: RelationCacheStatistics,
    /// RustRed's recipe counters: solves (and their seconds), memo hits,
    /// failures, solves outside a certified set.
    pub recipes: RecipeStatistics,
    /// Bytes of the persistent cache directory on disk (0 when absent).
    pub bytes_on_disk: u64,
}

pub(crate) fn statistics() -> PRelationCacheStatistics {
    let roots = recipe_cache_roots();
    PRelationCacheStatistics {
        bytes_on_disk: roots.user.as_deref().map_or(0, directory_bytes),
        directory: roots.user,
        shipped: roots.prefilled,
        cache: relation_cache_statistics(),
        recipes: recipe_statistics(),
    }
}

fn directory_bytes(path: &Path) -> u64 {
    let Ok(entries) = std::fs::read_dir(path) else {
        return 0;
    };
    entries
        .filter_map(Result::ok)
        .map(|entry| match entry.file_type() {
            Ok(kind) if kind.is_dir() => directory_bytes(&entry.path()),
            Ok(_) => entry.metadata().map_or(0, |metadata| metadata.len()),
            Err(_) => 0,
        })
        .sum()
}

/// vk5l WP4: how the keys of a reduction were served, by kind (plan
/// ARTIFACT_THREE_KINDS): `table` rule steps, `explicit` key expressions
/// (shipped relation entries), `cached` relations of the persistent cache,
/// `recipe_solved` recipe solves at run time (inside a certified set one per
/// key; outside it one per attempted recipe family, failed attempts
/// included), `refused` typed refusals. Counted per reduction call
/// ([`crate::Vakint::p_dispatch_counters`] for the process).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PDispatchCounters {
    pub table: u64,
    pub explicit: u64,
    pub cached: u64,
    pub recipe_solved: u64,
    pub refused: u64,
}

impl PDispatchCounters {
    /// `self - before`.
    pub fn since(&self, before: &Self) -> Self {
        Self {
            table: self.table - before.table,
            explicit: self.explicit - before.explicit,
            cached: self.cached - before.cached,
            recipe_solved: self.recipe_solved - before.recipe_solved,
            refused: self.refused - before.refused,
        }
    }
}

static TABLE_STEPS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static REFUSALS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// Table rule steps of a reduction (its rule applications).
pub(crate) fn count_table_steps(steps: u64) {
    TABLE_STEPS.fetch_add(steps, std::sync::atomic::Ordering::Relaxed);
}

/// A typed refusal of the p library.
pub(crate) fn count_refusal() {
    REFUSALS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
}

/// The dispatch counters of this process.
pub(crate) fn dispatch_counters() -> PDispatchCounters {
    let cache = relation_cache_statistics();
    let recipes = recipe_statistics();
    PDispatchCounters {
        table: TABLE_STEPS.load(std::sync::atomic::Ordering::Relaxed),
        explicit: cache.shipped_hits,
        cached: cache.cache_hits,
        recipe_solved: recipes.solves + recipes.outside_solves,
        refused: REFUSALS.load(std::sync::atomic::Ordering::Relaxed),
    }
}

/// One request's relation cache counters, for the per-request diagnostics.
pub(crate) fn describe(counters: &RelationCacheStatistics) -> String {
    format!(
        "relation_cache shipped_hits={} cache_hits={} misses={} inserts={} inserts_skipped_cap={} \
         foreign={} corrupt={} bytes_written={}",
        counters.shipped_hits,
        counters.cache_hits,
        counters.misses,
        counters.inserts,
        counters.inserts_skipped_cap,
        counters.foreign,
        counters.corrupt,
        counters.bytes_written
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_default_location_is_next_to_the_ordered_cache() {
        let settings = PRelationCacheSettings::default();
        assert_eq!(
            settings.resolved_directory(Some(Path::new("/x/ordered"))),
            Some(PathBuf::from("/x/ordered/p_relations"))
        );
        let moved = PRelationCacheSettings {
            location: PRelationCacheLocation::Directory("/y/rel".into()),
            ..Default::default()
        };
        assert_eq!(
            moved.resolved_directory(Some(Path::new("/x/ordered"))),
            Some(PathBuf::from("/y/rel"))
        );
        let off = PRelationCacheSettings {
            location: PRelationCacheLocation::Disabled,
            ..Default::default()
        };
        assert_eq!(off.resolved_directory(Some(Path::new("/x/ordered"))), None);
    }
}
