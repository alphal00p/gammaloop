//! The RustRed p library of one engine (M2 plan, sections 3.1 and 6, M2-V3).
//!
//! An engine reduces p-integrals with ONE RustRed library program whose scope
//! (a range of loop counts) covers the request: one program, one census, one
//! step cache for every loop count of its scope. It comes from, in this order:
//!
//! 1. a library injected into the engine ([`Vakint::inject_p_library`], tests
//!    and callers that build their own program, e.g.
//!    `rustred::propagator::testing::one_loop_reference_library()`);
//! 2. the external data the engine loaded: `<root>/p/p-library-<a>-<b>/` under
//!    the external-data root (the layout `rustred p-package` writes), so the
//!    4L release of M3/M4 needs no new loader;
//! 3. the built-in library embedded from `data/p/` (shipped by M2-V4; none in
//!    this build, so 1L-3L p-integrals keep going to Forcer).
//!
//! Directory and embedded libraries go through `rustred_app`'s fail-closed
//! loader (`load_p_library`). Every library, injected ones included, is
//! checked once before its first use ([`PLibraryHandle::new`]): the
//! coefficient context `{d}`, one top presentation per vakint family of its
//! scope whose line slots are vakint's canonical lines, and terminals that map
//! onto `1` or a tabulated Forcer master only. Nothing is ever recomputed from
//! a digest; the cache identity is textual ([`p_cache_directory_name`]).
//!
//! RustRed accessors used beyond the frozen `PLibrary` interface (plan
//! section 4): `terminals()` and `basis_map()` (the load-time terminal check).

use std::collections::BTreeMap;
use std::fmt;
use std::fs;
use std::ops::RangeInclusive;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, PoisonError};

use rustred::propagator::PresentationId;
use rustred::propagator::library::{PBasis, PLibrary};
use rustred::propagator::presentation::Slot;
use rustred::propagator::tag::GammaTag;
use rustred_app::propagator::{PLibraryOptions, PLibrarySource, load_p_library};
use symbolica::atom::Atom;

use super::RustRedEvaluationError;
use super::cache::OrderedCacheProgram;
use crate::Vakint;
use crate::p_integral::{P_FAMILIES, P_INTEGRAL_LOOP_COUNTS, masters, parts};

/// Directory of the p area under an external-data root.
pub(crate) const P_DATA_DIRECTORY: &str = "p";
/// Prefix of a library directory, `p-library-<a>-<b>` for the scope `a..=b`.
pub(crate) const P_LIBRARY_DIRECTORY_PREFIX: &str = "p-library-";
/// Longest cache directory name: below `NAME_MAX` (255) with a margin.
const MAX_CACHE_DIRECTORY_NAME: usize = 200;

fn data_failure(detail: impl fmt::Display) -> RustRedEvaluationError {
    RustRedEvaluationError::ExternalData {
        detail: format!("p library: {detail}"),
    }
}

/// Where an engine's library comes from (reported as the cache's data version).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PLibraryOrigin {
    Injected,
    ExternalData,
    #[allow(dead_code, reason = "M2-V4 wires the built-in Memory source")]
    BuiltIn,
}

impl PLibraryOrigin {
    fn as_str(self) -> &'static str {
        match self {
            Self::Injected => "injected",
            Self::ExternalData => "external-data",
            Self::BuiltIn => "built-in",
        }
    }

    /// Suffix of the cache directory name: an injected library is any program
    /// a caller built, so its identity (scope and versions) does not tell it
    /// from the shipped library of the same scope; its steps get their own
    /// directory. Shipped (external-data or built-in) libraries of the same
    /// identity are the same program and share one.
    fn cache_suffix(self) -> Option<&'static str> {
        match self {
            Self::Injected => Some("injected"),
            Self::ExternalData | Self::BuiltIn => None,
        }
    }
}

/// The top presentation of one vakint family in the library.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PTop {
    pub presentation: PresentationId,
    /// Line slots first (vakint's canonical lines, in order), then ISPs.
    pub line_count: usize,
    pub arity: usize,
}

/// A checked library and what vakint needs to call it.
pub(crate) struct PLibraryHandle {
    pub library: Arc<PLibrary>,
    pub origin: PLibraryOrigin,
    /// One entry per [`P_FAMILIES`] entry; `None` outside the scope.
    tops: Vec<Option<PTop>>,
    /// The library's dimension parameter `d` as an expression.
    pub dimension: Atom,
    /// The step-cache program of this library (one for every loop count).
    pub cache: OrderedCacheProgram,
}

impl fmt::Debug for PLibraryHandle {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("PLibraryHandle")
            .field("scope", &self.library.scope())
            .field("origin", &self.origin)
            .field("cache", &self.cache.label)
            .finish_non_exhaustive()
    }
}

/// `p-library-<a>-<b>_<v>_<v>...`: the cache directory of a library program,
/// a textual rendering of [`PLibrary::identity`] (scope and the value of
/// every `P_*` version constant) followed by RustRed's Ordered step-cache
/// schema and step semantics, each value without its `rustred` prefix, and
/// `_injected` for an injected library ([`PLibraryOrigin::cache_suffix`]). No
/// digest: a version bump names a new directory. Values are made file-name
/// safe; a name longer than 200 bytes is an error (fail closed rather than an
/// unusable path).
pub(crate) fn p_cache_directory_name(
    identity: &str,
    origin: PLibraryOrigin,
) -> Result<String, RustRedEvaluationError> {
    let malformed = || data_failure(format!("unexpected library identity {identity:?}"));
    let rest = identity
        .strip_prefix("p-library:scope=")
        .ok_or_else(malformed)?;
    let (scope, versions) = rest.split_once(':').ok_or_else(malformed)?;
    let mut values = Vec::new();
    for entry in versions.split(',') {
        let (name, value) = entry.split_once('=').ok_or_else(malformed)?;
        if name.is_empty() || value.is_empty() {
            return Err(malformed());
        }
        values.push(value);
    }
    values.push(rustred::reduction::ORDERED_STEP_CACHE_SCHEMA);
    values.push(rustred::reduction::ORDERED_STEP_SEMANTICS);
    values.extend(origin.cache_suffix());
    let safe = |text: &str| -> String {
        let text = text
            .strip_prefix("rustred.")
            .or_else(|| text.strip_prefix("rustred-"))
            .unwrap_or(text);
        text.chars()
            .map(|c| {
                if c.is_ascii_alphanumeric() || matches!(c, '-' | '.') {
                    c
                } else {
                    '-'
                }
            })
            .collect()
    };
    let mut name = format!("{P_LIBRARY_DIRECTORY_PREFIX}{}", safe(scope));
    for value in values {
        name.push('_');
        name.push_str(&safe(value));
    }
    if name.len() > MAX_CACHE_DIRECTORY_NAME {
        return Err(data_failure(format!(
            "the cache directory name of {identity:?} has {} bytes (at most {MAX_CACHE_DIRECTORY_NAME})",
            name.len()
        )));
    }
    Ok(name)
}

impl PLibraryHandle {
    /// Check `library` once for vakint (module documentation).
    pub(crate) fn new(
        library: Arc<PLibrary>,
        origin: PLibraryOrigin,
    ) -> Result<Self, RustRedEvaluationError> {
        let scope = library.scope();
        let scope_usize = usize::from(*scope.start())..=usize::from(*scope.end());
        if scope.is_empty()
            || !P_INTEGRAL_LOOP_COUNTS.contains(scope_usize.start())
            || !P_INTEGRAL_LOOP_COUNTS.contains(scope_usize.end())
        {
            return Err(data_failure(format!(
                "scope {scope:?} is outside the p-integral loop counts {P_INTEGRAL_LOOP_COUNTS:?}"
            )));
        }
        let context = library.coefficient_context();
        if context.parameter_names() != ["d"] {
            return Err(data_failure(format!(
                "coefficient variables {:?}, expected [\"d\"]",
                context.parameter_names()
            )));
        }
        let dimension = context
            .parameter("d")
            .ok_or_else(|| data_failure("no dimension parameter d"))?
            .to_expression();
        let mut tops = Vec::with_capacity(P_FAMILIES.len());
        for family in &P_FAMILIES {
            if !scope_usize.contains(&family.loops) {
                tops.push(None);
                continue;
            }
            tops.push(Some(Self::check_top(&library, family)?));
        }
        Self::check_terminals(&library)?;
        let identity = library.identity();
        let binding = library
            .step_cache_binding()
            .map_err(RustRedEvaluationError::from_native)?;
        let cache = OrderedCacheProgram::textual(
            *scope.end() as usize,
            &format!(
                "{P_LIBRARY_DIRECTORY_PREFIX}{}-{}",
                scope.start(),
                scope.end()
            ),
            origin.as_str(),
            p_cache_directory_name(&identity, origin)?,
            binding,
        );
        Ok(Self {
            library,
            origin,
            tops,
            dimension,
            cache,
        })
    }

    /// The top of `family`: a presentation of the family's loop count whose
    /// line slots are vakint's canonical lines (momenta up to an overall sign,
    /// no eps units), followed by ISP slots only.
    fn check_top(
        library: &PLibrary,
        family: &crate::p_integral::PFamilyData,
    ) -> Result<PTop, RustRedEvaluationError> {
        let id = library
            .top_presentation(family.name)
            .ok_or_else(|| data_failure(format!("no top presentation {}", family.name)))?;
        // A top id is one of the library's own presentations (`PLibrary` builds
        // its tops from them).
        check_top_presentation(library.presentation(id), id, family)
    }
}

/// [`PLibraryHandle::check_top`] on one presentation `id` of `family`.
pub(crate) fn check_top_presentation(
    presentation: &rustred::propagator::presentation::Presentation,
    id: PresentationId,
    family: &crate::p_integral::PFamilyData,
) -> Result<PTop, RustRedEvaluationError> {
    let mismatch =
        |detail: String| data_failure(format!("top {} of the library: {detail}", family.name));
    if usize::from(presentation.loops) != family.loops {
        return Err(mismatch(format!("{} loops", presentation.loops)));
    }
    let lines = parts::canonical_lines(family)
        .map_err(|error| mismatch(format!("canonical lines: {error}")))?;
    let arity = presentation.slots.len();
    if arity != lines.len() + family.isps.len() {
        return Err(mismatch(format!(
            "{arity} slots for {} lines and {} ISPs",
            lines.len(),
            family.isps.len()
        )));
    }
    for (index, slot) in presentation.slots.iter().enumerate() {
        match (slot, lines.get(index)) {
            (Slot::Line(line), Some(canonical)) => {
                let mut expected = canonical.loops.clone();
                expected.push(canonical.external);
                let negated = expected.iter().map(|c| -c).collect::<Vec<_>>();
                if line.eps_units != 0 || (line.momentum != expected && line.momentum != negated) {
                    return Err(mismatch(format!(
                        "line slot {} is {:?} (eps units {}), vakint's line is {expected:?}",
                        index + 1,
                        line.momentum,
                        line.eps_units
                    )));
                }
            }
            (Slot::Isp(_), None) => {}
            _ => {
                return Err(mismatch(format!(
                    "slot {} is not a {}",
                    index + 1,
                    if index < lines.len() { "line" } else { "ISP" }
                )));
            }
        }
    }
    Ok(PTop {
        presentation: id,
        line_count: lines.len(),
        arity,
    })
}

impl PLibraryHandle {
    /// Every terminal maps onto `1` or a tabulated Forcer master: without a
    /// basis map only the unit terminal may exist, and a residual never
    /// reaches vakint (convention C4).
    fn check_terminals(library: &PLibrary) -> Result<(), RustRedEvaluationError> {
        let names = masters::p_master_names().map_err(|error| data_failure(error.to_string()))?;
        let terminals = library.terminals().len();
        let Some(map) = library.basis_map() else {
            return if terminals == 1 {
                Ok(())
            } else {
                Err(data_failure(format!(
                    "{terminals} terminals but no basis map onto the Forcer masters"
                )))
            };
        };
        for id in 0..terminals {
            let terminal = rustred::propagator::PTerminalId(
                u32::try_from(id).map_err(|_| data_failure("terminal index overflow"))?,
            );
            let row = map
                .map(terminal, &GammaTag::one())
                .map_err(RustRedEvaluationError::from_native)?;
            if row.is_empty() {
                return Err(data_failure(format!("terminal {id} has no basis row")));
            }
            for term in row {
                match term.basis {
                    PBasis::Unit => {}
                    PBasis::Master(name) if names.contains(&name.as_str()) => {}
                    PBasis::Master(name) => {
                        return Err(RustRedEvaluationError::MasterCapability {
                            detail: format!(
                                "terminal {id} maps onto the p-master {name}, which the vendored \
                                 Forcer Masters.prc does not tabulate"
                            ),
                        });
                    }
                    PBasis::Residual(residual) => {
                        return Err(RustRedEvaluationError::MasterCapability {
                            detail: format!(
                                "terminal {id} maps onto the residual terminal {}",
                                residual.0
                            ),
                        });
                    }
                }
            }
        }
        Ok(())
    }

    pub(crate) fn scope(&self) -> RangeInclusive<u8> {
        self.library.scope()
    }

    /// The top of the [`P_FAMILIES`] entry `family`.
    pub(crate) fn top(&self, family: usize) -> Result<PTop, RustRedEvaluationError> {
        self.tops.get(family).copied().flatten().ok_or_else(|| {
            RustRedEvaluationError::UnsupportedMatchedFamily {
                detail: format!(
                    "the RustRed p library of scope {:?} has no family {}",
                    self.library.scope(),
                    P_FAMILIES.get(family).map_or("?", |family| family.name)
                ),
            }
        })
    }
}

/// A library directory under an external-data root.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct PLibraryDirectory {
    /// `<root>/p/p-library-<a>-<b>`.
    pub path: PathBuf,
    pub scope: RangeInclusive<u8>,
}

/// `a..=b` of a directory name `p-library-<a>-<b>`.
fn directory_scope(name: &str) -> Option<RangeInclusive<u8>> {
    let (start, end) = name
        .strip_prefix(P_LIBRARY_DIRECTORY_PREFIX)?
        .split_once('-')?;
    let canonical = |text: &str| !text.is_empty() && !text.starts_with('0');
    if !canonical(start) || !canonical(end) {
        return None;
    }
    let (start, end) = (start.parse::<u8>().ok()?, end.parse::<u8>().ok()?);
    (start <= end).then_some(start..=end)
}

/// The library directory of an external-data root: `None` without a `p/`
/// area; an error when the root is missing, `p/` is not a plain directory, or
/// it holds no or several `p-library-<a>-<b>` directories (one engine reduces
/// with one library program).
pub(crate) fn resolve_directory(
    root: &Path,
) -> Result<Option<PLibraryDirectory>, RustRedEvaluationError> {
    if !root.is_dir() {
        return Err(data_failure(format!(
            "external-data root {} is not a directory",
            root.display()
        )));
    }
    let area = root.join(P_DATA_DIRECTORY);
    let metadata = match fs::symlink_metadata(&area) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(error) => return Err(data_failure(format!("{}: {error}", area.display()))),
    };
    if !metadata.is_dir() {
        return Err(data_failure(format!(
            "{} is not a directory (symlinks are refused)",
            area.display()
        )));
    }
    let mut found = Vec::new();
    for entry in
        fs::read_dir(&area).map_err(|error| data_failure(format!("{}: {error}", area.display())))?
    {
        let entry = entry.map_err(|error| data_failure(format!("{}: {error}", area.display())))?;
        let name = entry.file_name();
        let Some(scope) = name.to_str().and_then(directory_scope) else {
            continue;
        };
        let file_type = entry
            .file_type()
            .map_err(|error| data_failure(format!("{}: {error}", entry.path().display())))?;
        if !file_type.is_dir() {
            return Err(data_failure(format!(
                "{} is not a directory (symlinks are refused)",
                entry.path().display()
            )));
        }
        found.push(PLibraryDirectory {
            path: entry.path(),
            scope,
        });
    }
    found.sort_by(|a, b| a.path.cmp(&b.path));
    match found.len() {
        0 => Err(data_failure(format!(
            "{} holds no {P_LIBRARY_DIRECTORY_PREFIX}<a>-<b> library directory",
            area.display()
        ))),
        1 => Ok(found.pop()),
        _ => Err(data_failure(format!(
            "{} holds several library directories ({}); one engine uses one library program",
            area.display(),
            found
                .iter()
                .map(|directory| directory.path.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        ))),
    }
}

/// Load a library directory through `rustred_app`'s fail-closed loader and
/// check that its scope is the one its name declares.
pub(crate) fn load_directory(
    directory: &PLibraryDirectory,
) -> Result<Arc<PLibraryHandle>, RustRedEvaluationError> {
    if !directory.path.is_dir() {
        return Err(data_failure(format!(
            "missing library directory {}",
            directory.path.display()
        )));
    }
    let library = load_p_library(
        PLibrarySource::Directory(directory.path.clone()),
        &PLibraryOptions {
            verify_descent: true,
        },
    )
    .map_err(|error| data_failure(format!("{}: {error}", directory.path.display())))?;
    if library.scope() != directory.scope {
        return Err(data_failure(format!(
            "{} holds a library of scope {:?}",
            directory.path.display(),
            library.scope()
        )));
    }
    PLibraryHandle::new(library, PLibraryOrigin::ExternalData).map(Arc::new)
}

/// Libraries loaded from external-data directories, by directory: loaded once
/// per process, like the vendored data. Failures are not kept.
fn loaded_directories() -> &'static Mutex<BTreeMap<PathBuf, Arc<PLibraryHandle>>> {
    static LOADED: Mutex<BTreeMap<PathBuf, Arc<PLibraryHandle>>> = Mutex::new(BTreeMap::new());
    &LOADED
}

/// Library directories found under external-data roots, by root. Only found
/// directories are kept: a root without a `p/` area costs two `stat` calls
/// per lookup and picks up an area created later.
fn resolved_roots() -> &'static Mutex<BTreeMap<PathBuf, PLibraryDirectory>> {
    static RESOLVED: Mutex<BTreeMap<PathBuf, PLibraryDirectory>> = Mutex::new(BTreeMap::new());
    &RESOLVED
}

/// The built-in library embedded from `data/p/` (M2-V4 wires the Memory
/// source; this build embeds none).
fn built_in() -> Result<Option<Arc<PLibraryHandle>>, RustRedEvaluationError> {
    Ok(None)
}

/// What [`Vakint::p_library_scope`] found without loading anything.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum PDataPresence {
    None,
    Scope(RangeInclusive<u8>),
    /// A `p/` area that cannot be resolved: p-integrals fail closed with its error.
    Broken,
}

impl Vakint {
    /// Use `library` for the p-integrals of this engine (and its later clones),
    /// ahead of external and built-in data. The library is checked first
    /// (coefficient context, family tops, terminals). Not part of the API:
    /// tests and M2 consumers inject e.g.
    /// `rustred::propagator::testing::one_loop_reference_library()`.
    #[doc(hidden)]
    pub fn inject_p_library(&mut self, library: Arc<PLibrary>) -> Result<(), crate::VakintError> {
        self.p_library = Some(Arc::new(PLibraryHandle::new(
            library,
            PLibraryOrigin::Injected,
        )?));
        Ok(())
    }

    /// The external-data root of this engine's loaded data, if it was read
    /// from a directory.
    fn external_data_root(&self) -> Option<PathBuf> {
        let data = self.rustred_data.as_ref()?;
        data.loops.values().find_map(|pack| {
            let manifest = pack.path("manifest.json").ok().flatten()?;
            Some(manifest.parent()?.parent()?.to_owned())
        })
    }

    fn resolved_directory(&self) -> Result<Option<PLibraryDirectory>, RustRedEvaluationError> {
        let Some(root) = self.external_data_root() else {
            return Ok(None);
        };
        if let Some(known) = resolved_roots()
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .get(&root)
        {
            return Ok(Some(known.clone()));
        }
        let resolved = resolve_directory(&root)?;
        if let Some(directory) = &resolved {
            resolved_roots()
                .lock()
                .unwrap_or_else(PoisonError::into_inner)
                .insert(root, directory.clone());
        }
        Ok(resolved)
    }

    /// The scope of the library this engine would use, without loading it.
    pub(crate) fn p_library_presence(&self) -> PDataPresence {
        if let Some(handle) = &self.p_library {
            return PDataPresence::Scope(handle.scope());
        }
        match self.resolved_directory() {
            Ok(Some(directory)) => return PDataPresence::Scope(directory.scope),
            Ok(None) => {}
            Err(_) => return PDataPresence::Broken,
        }
        match built_in() {
            Ok(Some(handle)) => PDataPresence::Scope(handle.scope()),
            Ok(None) => PDataPresence::None,
            Err(_) => PDataPresence::Broken,
        }
    }

    /// The library of this engine (module documentation), loaded on first use.
    pub(crate) fn p_library(&self) -> Result<Arc<PLibraryHandle>, RustRedEvaluationError> {
        if let Some(handle) = &self.p_library {
            return Ok(Arc::clone(handle));
        }
        if let Some(directory) = self.resolved_directory()? {
            if let Some(handle) = loaded_directories()
                .lock()
                .unwrap_or_else(PoisonError::into_inner)
                .get(&directory.path)
            {
                return Ok(Arc::clone(handle));
            }
            let handle = load_directory(&directory)?;
            loaded_directories()
                .lock()
                .unwrap_or_else(PoisonError::into_inner)
                .insert(directory.path.clone(), Arc::clone(&handle));
            return Ok(handle);
        }
        built_in()?.ok_or_else(|| RustRedEvaluationError::UnsupportedMatchedFamily {
            detail: "no RustRed p library is loaded (the built-in 1L-3L library ships with M2-V4; \
                     external data provide p/p-library-<a>-<b>/)"
                .into(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reference() -> Arc<PLibrary> {
        rustred::propagator::testing::one_loop_reference_library().unwrap()
    }

    #[test]
    fn the_reference_library_is_accepted_with_its_bubble_top() {
        let handle = PLibraryHandle::new(reference(), PLibraryOrigin::Injected).unwrap();
        assert_eq!(handle.scope(), 1..=1);
        let p1l = P_FAMILIES.iter().position(|f| f.name == "P1L").unwrap();
        let top = handle.top(p1l).unwrap();
        assert_eq!((top.line_count, top.arity), (2, 2));
        let p2l = P_FAMILIES.iter().position(|f| f.name == "P2L").unwrap();
        assert!(matches!(
            handle.top(p2l),
            Err(RustRedEvaluationError::UnsupportedMatchedFamily { .. })
        ));
        assert_eq!(handle.cache.label, "p-library-1-1");
        assert_eq!(handle.cache.loop_count, 1);
    }

    /// The cache directory renders the library identity textually, with
    /// RustRed's step-cache versions, well below NAME_MAX; an injected
    /// library gets its own directory.
    #[test]
    fn the_cache_directory_is_a_textual_rendering_of_the_identity() {
        let library = reference();
        let identity = library.identity();
        let directory =
            |identity: &str| p_cache_directory_name(identity, PLibraryOrigin::ExternalData);
        let name = directory(&identity).unwrap();
        assert!(name.starts_with("p-library-1-1_p-owner-library-"), "{name}");
        for (_, value) in rustred::propagator::p_versions() {
            let value = value
                .strip_prefix("rustred.")
                .or_else(|| value.strip_prefix("rustred-"))
                .unwrap_or(&value)
                .to_owned();
            assert!(name.contains(&value), "{value} missing in {name}");
        }
        assert!(name.ends_with(&format!(
            "_{}",
            rustred::reduction::ORDERED_STEP_SEMANTICS.trim_start_matches("rustred-")
        )));
        assert_eq!(
            p_cache_directory_name(&identity, PLibraryOrigin::BuiltIn).unwrap(),
            name
        );
        let injected = p_cache_directory_name(&identity, PLibraryOrigin::Injected).unwrap();
        assert_eq!(injected, format!("{name}_injected"));
        assert!(
            injected.len() <= MAX_CACHE_DIRECTORY_NAME,
            "{} bytes",
            injected.len()
        );
        assert!(
            injected
                .chars()
                .all(|c| c.is_ascii_alphanumeric() || matches!(c, '-' | '.' | '_'))
        );
        // A different scope or version is a different directory.
        let other = identity.replacen("scope=1-1", "scope=1-3", 1);
        assert_ne!(directory(&other).unwrap(), name);
        let bumped = identity.replacen("P_CENSUS_VERSION=1", "P_CENSUS_VERSION=2", 1);
        assert_ne!(bumped, identity);
        assert_ne!(directory(&bumped).unwrap(), name);
        for malformed in ["", "p-library:1-3", "p-library:scope=1-3:P_X"] {
            assert!(directory(malformed).is_err(), "{malformed}");
        }
        let long = format!("{identity},P_LONG={}", "x".repeat(200));
        assert!(matches!(
            directory(&long),
            Err(RustRedEvaluationError::ExternalData { .. })
        ));
    }

    /// An injected library and a loaded one of the same identity never share
    /// a step cache, on disk or in memory (one engine, one registry).
    #[test]
    fn an_injected_library_has_its_own_step_cache() {
        let root = tempfile::tempdir().unwrap();
        let settings = crate::VakintSettings {
            ordered_cache: crate::OrderedCacheSettings::read_write(root.path()),
            ..Default::default()
        };
        let engine = Vakint::new().unwrap();
        let injected = PLibraryHandle::new(reference(), PLibraryOrigin::Injected).unwrap();
        let loaded = PLibraryHandle::new(reference(), PLibraryOrigin::ExternalData).unwrap();
        let caches = [&injected, &loaded].map(|handle| {
            engine
                .ordered_step_cache(&settings, &handle.cache)
                .unwrap()
                .unwrap()
        });
        assert!(!Arc::ptr_eq(&caches[0], &caches[1]));
        let directories = engine
            .ordered_cache_statistics()
            .programs
            .into_iter()
            .map(|program| (program.data_version, program.directory.unwrap()))
            .collect::<Vec<_>>();
        let name = |handle: &PLibraryHandle| {
            root.path()
                .join(handle.cache.textual_tables.as_ref().unwrap())
        };
        assert_eq!(
            directories,
            [
                ("external-data".to_owned(), name(&loaded)),
                ("injected".to_owned(), name(&injected)),
            ]
        );
        assert_eq!(
            name(&injected).file_name().unwrap().to_str().unwrap(),
            format!(
                "{}_injected",
                name(&loaded).file_name().unwrap().to_str().unwrap()
            )
        );
    }

    #[test]
    fn library_directory_names_declare_their_scope() {
        assert_eq!(directory_scope("p-library-1-3"), Some(1..=3));
        assert_eq!(directory_scope("p-library-1-4"), Some(1..=4));
        for name in [
            "p-library-3-1",
            "p-library-01-3",
            "p-library-1-",
            "p-library-1",
            "library-1-3",
            "p-library-1-3x",
        ] {
            assert_eq!(directory_scope(name), None, "{name}");
        }
    }

    /// The Directory source of an external-data root (critique-phys SF8):
    /// `<root>/p/p-library-<a>-<b>/`, absent without a `p/` area, typed errors
    /// for a missing root, an ambiguous or empty area and a missing library
    /// directory; at this rustred pin the loader itself still refuses.
    #[test]
    fn directory_sources_resolve_under_the_external_data_root() {
        let root = tempfile::tempdir().unwrap();
        assert_eq!(resolve_directory(root.path()).unwrap(), None);
        let missing = root.path().join("missing-root");
        assert!(matches!(
            resolve_directory(&missing),
            Err(RustRedEvaluationError::ExternalData { ref detail }) if detail.contains("missing-root")
        ));

        let area = root.path().join(P_DATA_DIRECTORY);
        fs::create_dir(&area).unwrap();
        assert!(matches!(
            resolve_directory(root.path()),
            Err(RustRedEvaluationError::ExternalData { ref detail }) if detail.contains("no p-library-")
        ));
        fs::write(area.join("manifest.json"), b"{}").unwrap();
        fs::create_dir(area.join("p-library-1-3")).unwrap();
        let resolved = resolve_directory(root.path()).unwrap().unwrap();
        assert_eq!(resolved.path, area.join("p-library-1-3"));
        assert_eq!(resolved.scope, 1..=3);

        fs::create_dir(area.join("p-library-1-4")).unwrap();
        assert!(matches!(
            resolve_directory(root.path()),
            Err(RustRedEvaluationError::ExternalData { ref detail }) if detail.contains("several")
        ));
        fs::remove_dir(area.join("p-library-1-4")).unwrap();
        fs::write(area.join("p-library-1-4"), b"").unwrap();
        assert!(matches!(
            resolve_directory(root.path()),
            Err(RustRedEvaluationError::ExternalData { ref detail }) if detail.contains("not a directory")
        ));
        fs::remove_file(area.join("p-library-1-4")).unwrap();

        let gone = PLibraryDirectory {
            path: area.join("p-library-2-2"),
            scope: 2..=2,
        };
        assert!(matches!(
            load_directory(&gone),
            Err(RustRedEvaluationError::ExternalData { ref detail }) if detail.contains("missing library directory")
        ));
        // An empty directory is refused by rustred's loader (fail closed).
        assert!(matches!(
            load_directory(&resolved),
            Err(RustRedEvaluationError::ExternalData { ref detail }) if detail.contains("p-library-1-3")
        ));
    }

    /// A `p/` area of the loaded external data takes part in the dispatch:
    /// its scope counts as present data (so a broken library fails closed
    /// instead of silently going to Forcer), an injected library comes first.
    #[test]
    fn an_external_data_p_area_is_found_from_the_loaded_root() {
        let root = tempfile::tempdir().unwrap();
        for (path, bytes) in super::super::data::vendored_files() {
            let target = root.path().join(path);
            fs::create_dir_all(target.parent().unwrap()).unwrap();
            fs::write(target, bytes).unwrap();
        }
        let mut engine = Vakint::new().unwrap();
        engine.load_external_data(root.path(), 3).unwrap();
        assert_eq!(
            engine.external_data_root().unwrap(),
            root.path().canonicalize().unwrap()
        );
        assert_eq!(engine.p_library_presence(), PDataPresence::None);
        fs::create_dir_all(root.path().join("p/p-library-1-3")).unwrap();
        let mut engine = Vakint::new().unwrap();
        engine.load_external_data(root.path(), 3).unwrap();
        assert_eq!(engine.p_library_presence(), PDataPresence::Scope(1..=3));
        assert!(matches!(
            engine.p_library(),
            Err(RustRedEvaluationError::ExternalData { .. })
        ));
        engine.inject_p_library(reference()).unwrap();
        assert_eq!(engine.p_library_presence(), PDataPresence::Scope(1..=1));
        assert_eq!(engine.p_library().unwrap().origin, PLibraryOrigin::Injected);
    }
}
