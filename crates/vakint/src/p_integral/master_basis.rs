//! The five-loop p master basis (user ruling 2026-11-16: the master EXCHANGE
//! basis is the primary route for the missing coefficients).
//!
//! [`PMasterBasis::Original`] (the default) keeps the published GGPPSS masters
//! (arXiv:2104.08272, `MasterExpansions.m`) of the p master catalogue: nothing
//! in this module runs, no file is read and no symbol of it is created. With
//! [`PMasterBasis::PoleBounded`] the finalizer rewrites a resolved result,
//! exactly in `ε` and with the masters still symbolic, after the terminal map
//! ([`super::master_catalogue::PMasterSource::resolve_terminals`]) and before
//! the `ε` expansion, into the pole-bounded basis: every exchanged master
//! `PMaster(old)` becomes `Σ_n B(ε) PMaster(new_n) + Σ_k A(ε) PMaster(int_k)`
//! (GGPPSS `anc/BasisChange.m`: 31 masters exchanged for `ggppss_new(1..31)`;
//! further exchanges of ours, hep-ph/0601165 procedure, come as further
//! layers, applied in order). No reduction is redone.
//!
//! The data is five-loop data and ships with the five-loop package, in its
//! masters directory `<package>/masters/` (layout agreed with the package
//! assembler, lane vk5l-wp6):
//!
//! ```text
//! MANIFEST.json                              schema vk5l.p-masters-manifest.v1: catalogue_version,
//!                                            rotation_version (equal), basis_version, layers (apply
//!                                            order), optional layer_versions {layer: data_version}
//!                                            (a layer not listed carries basis_version), every
//!                                            file with role, sha256 and bytes
//! original/ggppss_catalogue.<v>.json.gz      the GGPPSS masters and terminal rows (catalogue schema)
//! pole-bounded/masters.<v>.json.gz           the new masters (schema vakint.p-master-basis-masters.v1:
//!                                            catalogue master entries + definitions)
//! basis-change/<NN>-<name>.json.gz           one layer each (schema vakint.p-master-basis-change.v1;
//!                                            later layers may exchange masters an earlier one
//!                                            introduced, and introduce masters of their own,
//!                                            e.g. `our(i)` of our hep-ph/0601165 exchanges; the
//!                                            masters file holds the new masters of every layer):
//!   { "schema", "version": 1, "data_version", "variable": "eps", "dimension": "d = 4 - 2*eps",
//!     "normalisation": "g-scheme ..", "layer": { "name", "source", "direction", "new_masters",
//!       "exchanges": [ { "old": "ggppss_int(116)", "terms": [ { "master": "ggppss_new(1)",
//!                                                    "num": ["3118348800", ..], "den": [..] } ] } ] } }
//! ```
//!
//! `num` / `den` are the integer coefficients of the numerator and the
//! denominator polynomial in `ε`, ascending powers. A term may carry
//! `"tag": [[m, e], ..]` (the factor `Π PGamma(m)^e`) and `"master": null` (a
//! closed part): our exchanges (helper VH 2 `our_exchanges`) have both. Every file vakint reads is
//! checked against the manifest's sha256 and size; a different
//! catalogue/rotation version is refused. The new masters' values (catalogue
//! schema: exact through their known order, unknown from the next one) extend
//! the catalogue only when the basis is selected; an unknown order is the
//! named `Oep(k, PMaster(..))` tail, never 0.
//!
//! The directory is the loaded five-loop package's `masters/`
//! ([`P_PACKAGE_MASTERS_DIRECTORY`]); [`P_MASTERS_DIRECTORY_VARIABLE`] names
//! another one. Selecting the pole-bounded basis without one is an error
//! naming the expected path. Every result finalized in the pole-bounded basis
//! logs its [`PMasterBasisRecord`] (basis version, sha256 of the manifest,
//! layers); [`crate::Vakint::p_master_basis_record`] returns it for result
//! metadata.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, PoisonError};

use serde::Deserialize;
use sha2::{Digest, Sha256};
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::id::Replacement;

use super::master_catalogue::PMasterCatalogue;
use super::{p_gamma_symbol, p_master_symbol};
use crate::VakintError;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};

/// The schema of a master basis change layer file.
pub const P_MASTER_BASIS_SCHEMA: &str = "vakint.p-master-basis-change.v1";
/// The schema of the pole-bounded masters file.
pub const P_MASTER_BASIS_MASTERS_SCHEMA: &str = "vakint.p-master-basis-masters.v1";
/// The schema of the package masters manifest.
pub const P_MASTERS_MANIFEST_SCHEMA: &str = "vk5l.p-masters-manifest.v1";
/// The format version of these files this vakint reads.
pub const P_MASTER_BASIS_VERSION: u32 = 1;
/// The masters directory of a five-loop package.
pub const P_PACKAGE_MASTERS_DIRECTORY: &str = "masters";
/// Its manifest.
pub const P_MASTERS_MANIFEST_FILE: &str = "MANIFEST.json";
/// Environment variable naming a masters directory (instead of the loaded
/// package's).
pub const P_MASTERS_DIRECTORY_VARIABLE: &str = "VAKINT_P_MASTERS_DIRECTORY";

/// The master basis of five-loop p results
/// ([`crate::VakintSettings::p_master_basis`]).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum PMasterBasis {
    /// The published GGPPSS masters of the catalogue (unchanged behaviour).
    #[default]
    Original,
    /// GGPPSS's new basis plus our further exchanges: no reduction
    /// coefficient has a pole deeper than its master's known order.
    PoleBounded,
}

impl std::fmt::Display for PMasterBasis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Original => write!(f, "original"),
            Self::PoleBounded => write!(f, "pole-bounded"),
        }
    }
}

/// What a result finalized in the pole-bounded basis records about it.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
pub struct PMasterBasisRecord {
    /// `pole-bounded`.
    pub basis: String,
    /// The basis version (manifest `basis_version`).
    pub data_version: String,
    /// sha256 of the masters manifest (which binds every file by sha256).
    pub sha256: String,
    /// The layers applied, in order.
    pub layers: Vec<String>,
    /// The masters directory read.
    pub directory: String,
}

#[derive(Deserialize)]
struct Manifest {
    schema: String,
    catalogue_version: String,
    rotation_version: String,
    basis_version: String,
    layers: Vec<String>,
    /// The data version of a layer when it differs from `basis_version` (a
    /// layer kept byte-identical when a later one is added).
    #[serde(default)]
    layer_versions: BTreeMap<String, String>,
    files: Vec<ManifestFile>,
}

#[derive(Deserialize)]
struct ManifestFile {
    path: String,
    role: String,
    sha256: String,
    bytes: u64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct LayerDocument {
    schema: String,
    version: u32,
    data_version: String,
    variable: String,
    #[allow(dead_code)]
    dimension: String,
    normalisation: String,
    layer: LayerFile,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct MastersDocument {
    schema: String,
    version: u32,
    data_version: String,
    #[allow(dead_code)]
    normalisation: String,
    masters: serde_json::Value,
    #[serde(default)]
    #[allow(dead_code)]
    legend: serde_json::Value,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct LayerFile {
    name: String,
    #[allow(dead_code)]
    source: String,
    #[allow(dead_code)]
    direction: String,
    new_masters: Vec<String>,
    exchanges: Vec<ExchangeFile>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ExchangeFile {
    old: String,
    terms: Vec<TermFile>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct TermFile {
    /// `null`: a closed part (no master).
    master: Option<String>,
    num: Vec<String>,
    den: Vec<String>,
    /// `[[m, e], ..]`: the factor `Π PGamma(m)^e`, `PGamma(m) = Γ(1 + m ε)`.
    #[serde(default)]
    tag: Vec<[i64; 2]>,
}

/// One term of an exchange: `num(ε)/den(ε) · Π PGamma(m)^e · PMaster(master)`
/// (no master: a closed part).
#[derive(Debug, Clone)]
pub(crate) struct Term {
    pub(crate) master: Option<String>,
    pub(crate) num: Vec<String>,
    pub(crate) den: Vec<String>,
    pub(crate) tag: Vec<[i64; 2]>,
}

/// One exchange: `old = Σ coefficient(ε)·master` (catalogue ids).
#[derive(Debug, Clone)]
pub(crate) struct Exchange {
    pub(crate) old: String,
    /// Numerators and denominators: integer coefficients ascending in ε.
    pub(crate) terms: Vec<Term>,
}

/// One layer: its name, the masters it introduces, its exchanges.
#[derive(Debug, Clone)]
pub(crate) struct Layer {
    pub(crate) name: String,
    pub(crate) new_masters: Vec<String>,
    pub(crate) exchanges: Vec<Exchange>,
}

/// A parsed basis change (no symbol created yet).
#[derive(Debug)]
pub struct PMasterBasisData {
    source_name: String,
    data_version: String,
    sha256: String,
    pub(crate) layers: Vec<Layer>,
    masters: serde_json::Value,
    /// The catalogue extended by the new masters, per base catalogue (its
    /// data version and master count).
    extended: Mutex<BTreeMap<(String, usize), Arc<PMasterCatalogue>>>,
    /// The rotation rules per epsilon symbol name and base catalogue.
    rules: Mutex<RuleCache>,
}

/// Rewrite rules `(PMaster(old), value)` per (epsilon symbol name, catalogue
/// data version, catalogue master count).
type RuleCache = BTreeMap<(String, String, usize), Arc<Vec<(Atom, Atom)>>>;

fn invalid(source: &str, detail: impl std::fmt::Display) -> VakintError {
    VakintError::EvaluationError(format!("p master basis {source}: {detail}"))
}

impl PMasterBasisData {
    /// Read the masters directory `directory` (its manifest, then every
    /// basis-change layer in order and the pole-bounded masters, each checked
    /// against the manifest's sha256 and size).
    pub fn from_directory(directory: &Path) -> Result<Self, VakintError> {
        let source_name = directory.display().to_string();
        let source = source_name.as_str();
        let manifest_path = directory.join(P_MASTERS_MANIFEST_FILE);
        let manifest_bytes = std::fs::read(&manifest_path).map_err(|error| {
            invalid(
                source,
                format!(
                    "the pole-bounded master basis needs {} ({error})",
                    manifest_path.display()
                ),
            )
        })?;
        let manifest: Manifest = serde_json::from_slice(&manifest_bytes)
            .map_err(|error| invalid(source, format!("{P_MASTERS_MANIFEST_FILE}: {error}")))?;
        if manifest.schema != P_MASTERS_MANIFEST_SCHEMA {
            return Err(invalid(
                source,
                format!("manifest schema {}", manifest.schema),
            ));
        }
        if manifest.catalogue_version != manifest.rotation_version {
            return Err(invalid(
                source,
                format!(
                    "rotation_version {} differs from catalogue_version {}",
                    manifest.rotation_version, manifest.catalogue_version
                ),
            ));
        }
        let read = |path: &str, role: &str| -> Result<String, VakintError> {
            let entry = manifest
                .files
                .iter()
                .find(|f| f.path == path && f.role == role)
                .ok_or_else(|| {
                    invalid(source, format!("the manifest lists no {role} file {path}"))
                })?;
            let bytes = std::fs::read(directory.join(path))
                .map_err(|error| invalid(source, format!("{path}: {error}")))?;
            let sha256 = format!("{:x}", Sha256::digest(&bytes));
            if sha256 != entry.sha256 || bytes.len() as u64 != entry.bytes {
                return Err(invalid(
                    source,
                    format!(
                        "{path}: sha256 {sha256} ({} bytes) differs from the manifest's {} ({} bytes)",
                        bytes.len(),
                        entry.sha256,
                        entry.bytes
                    ),
                ));
            }
            let mut text = String::new();
            flate2::read::GzDecoder::new(bytes.as_slice())
                .read_to_string(&mut text)
                .map_err(|error| invalid(source, format!("{path}: not gzip JSON: {error}")))?;
            Ok(text)
        };
        let mut layers = Vec::new();
        for name in &manifest.layers {
            let path = format!("basis-change/{name}.json.gz");
            let document: LayerDocument = serde_json::from_str(&read(&path, "basis-change")?)
                .map_err(|error| invalid(source, format!("{path}: not a basis change: {error}")))?;
            if document.schema != P_MASTER_BASIS_SCHEMA
                || document.version != P_MASTER_BASIS_VERSION
            {
                return Err(invalid(
                    source,
                    format!(
                        "{path}: schema {} version {}, this vakint reads {P_MASTER_BASIS_SCHEMA} version {P_MASTER_BASIS_VERSION}",
                        document.schema, document.version
                    ),
                ));
            }
            let data_version = manifest
                .layer_versions
                .get(name)
                .unwrap_or(&manifest.basis_version);
            if document.data_version != *data_version
                || document.variable != "eps"
                || document.normalisation.split_whitespace().next() != Some("g-scheme")
                || document.layer.name != *name
            {
                return Err(invalid(
                    source,
                    format!(
                        "{path}: data {} / variable {} / normalisation {} / layer {} do not match the manifest ({}, eps, g-scheme, {name})",
                        document.data_version,
                        document.variable,
                        document.normalisation,
                        document.layer.name,
                        data_version
                    ),
                ));
            }
            let mut exchanges = Vec::new();
            for exchange in document.layer.exchanges {
                let terms = exchange
                    .terms
                    .into_iter()
                    .map(|term| {
                        if term.den.iter().all(|c| c.trim_start_matches('-') == "0") {
                            Err(invalid(
                                source,
                                format!(
                                    "{}: zero denominator of {}",
                                    exchange.old,
                                    term.master.as_deref().unwrap_or("its closed part")
                                ),
                            ))
                        } else {
                            Ok(Term {
                                master: term.master,
                                num: term.num,
                                den: term.den,
                                tag: term.tag,
                            })
                        }
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                exchanges.push(Exchange {
                    old: exchange.old,
                    terms,
                });
            }
            layers.push(Layer {
                name: name.clone(),
                new_masters: document.layer.new_masters,
                exchanges,
            });
        }
        let masters_path = manifest
            .files
            .iter()
            .find(|f| f.role == "pole-bounded-masters")
            .map(|f| f.path.clone())
            .ok_or_else(|| invalid(source, "the manifest lists no pole-bounded-masters file"))?;
        let masters: MastersDocument =
            serde_json::from_str(&read(&masters_path, "pole-bounded-masters")?).map_err(
                |error| {
                    invalid(
                        source,
                        format!("{masters_path}: not a masters file: {error}"),
                    )
                },
            )?;
        if masters.schema != P_MASTER_BASIS_MASTERS_SCHEMA
            || masters.version != P_MASTER_BASIS_VERSION
            || masters.data_version != manifest.basis_version
        {
            return Err(invalid(
                source,
                format!(
                    "{masters_path}: schema {} version {} data {} (expected {P_MASTER_BASIS_MASTERS_SCHEMA} {P_MASTER_BASIS_VERSION} {})",
                    masters.schema, masters.version, masters.data_version, manifest.basis_version
                ),
            ));
        }
        let ids: Vec<&str> = masters
            .masters
            .as_array()
            .map(|list| list.iter().filter_map(|m| m["id"].as_str()).collect())
            .unwrap_or_default();
        for layer in &layers {
            if let Some(missing) = layer
                .new_masters
                .iter()
                .find(|id| !ids.contains(&id.as_str()))
            {
                return Err(invalid(
                    source,
                    format!(
                        "layer {}: new master {missing} is not in {masters_path}",
                        layer.name
                    ),
                ));
            }
        }
        Ok(Self {
            source_name,
            data_version: manifest.basis_version,
            sha256: format!("{:x}", Sha256::digest(&manifest_bytes)),
            layers,
            masters: masters.masters,
            extended: Mutex::new(BTreeMap::new()),
            rules: Mutex::new(BTreeMap::new()),
        })
    }

    /// The record of this basis (result metadata).
    pub fn record(&self) -> PMasterBasisRecord {
        PMasterBasisRecord {
            basis: PMasterBasis::PoleBounded.to_string(),
            data_version: self.data_version.clone(),
            sha256: self.sha256.clone(),
            layers: self.layers.iter().map(|layer| layer.name.clone()).collect(),
            directory: self.source_name.clone(),
        }
    }

    /// `catalogue` with this basis's new masters (cached per catalogue).
    pub(crate) fn extended_catalogue(
        &self,
        catalogue: &PMasterCatalogue,
    ) -> Result<Arc<PMasterCatalogue>, VakintError> {
        let key = (catalogue.data_version().to_owned(), catalogue.counts().0);
        let mut cache = self.extended.lock().unwrap_or_else(PoisonError::into_inner);
        if let Some(extended) = cache.get(&key) {
            return Ok(Arc::clone(extended));
        }
        let extended = Arc::new(
            catalogue
                .with_masters_json(&self.source_name, &self.masters)
                .map_err(VakintError::from)?,
        );
        cache.insert(key, Arc::clone(&extended));
        Ok(extended)
    }

    /// The rewrite rules `PMaster(old) -> Σ coefficient(epsilon)·PMaster(..)`
    /// of every layer, in order, for `catalogue` (already extended).
    fn rules(
        &self,
        catalogue: &PMasterCatalogue,
        epsilon: Symbol,
    ) -> Result<Arc<Vec<(Atom, Atom)>>, VakintError> {
        let key = (
            epsilon.get_name().to_owned(),
            catalogue.data_version().to_owned(),
            catalogue.counts().0,
        );
        let mut cache = self.rules.lock().unwrap_or_else(PoisonError::into_inner);
        if let Some(rules) = cache.get(&key) {
            return Ok(Arc::clone(rules));
        }
        let symbol = |id: &str| {
            catalogue
                .symbol_of_id(id)
                .map(str::to_owned)
                .ok_or_else(|| {
                    invalid(
                        &self.source_name,
                        format!(
                            "master {id} is not in the catalogue {}",
                            catalogue.data_version()
                        ),
                    )
                })
        };
        let placeholder = Atom::var(vk_symbol!("p_master_basis_eps"));
        let epsilon_atom = Atom::var(epsilon);
        let polynomial = |coefficients: &[String]| -> Result<Atom, VakintError> {
            let text = coefficients
                .iter()
                .enumerate()
                .filter(|(_, c)| c.as_str() != "0")
                .map(|(power, c)| format!("({c})*p_master_basis_eps^{power}"))
                .collect::<Vec<_>>()
                .join("+");
            let text = if text.is_empty() {
                "0".to_owned()
            } else {
                text
            };
            Ok(vk_parse!(&text)
                .map_err(|error| invalid(&self.source_name, error))?
                .replace(placeholder.to_pattern())
                .with(epsilon_atom.to_pattern()))
        };
        let mut rules = Vec::new();
        for layer in &self.layers {
            for exchange in &layer.exchanges {
                let old =
                    p_master_symbol().call(&[Atom::var(vk_symbol!(&symbol(&exchange.old)?))][..]);
                let mut value = Atom::Zero;
                for term in &exchange.terms {
                    let mut factor = polynomial(&term.num)? / polynomial(&term.den)?;
                    for [m, e] in &term.tag {
                        factor *= p_gamma_symbol()
                            .call(&[Atom::num(*m)][..])
                            .pow(Atom::num(*e));
                    }
                    if let Some(master) = &term.master {
                        factor *=
                            p_master_symbol().call(&[Atom::var(vk_symbol!(&symbol(master)?))][..]);
                    }
                    value += factor;
                }
                rules.push((old, value));
            }
        }
        let rules = Arc::new(rules);
        cache.insert(key, Arc::clone(&rules));
        Ok(rules)
    }

    /// `expression` (exact in `epsilon`, masters `PMaster(..)` symbolic) in
    /// the pole-bounded basis: every layer's exchanges substituted, in order.
    /// `catalogue` must be [`Self::extended_catalogue`].
    pub(crate) fn rotate(
        &self,
        expression: &Atom,
        catalogue: &PMasterCatalogue,
        epsilon: Symbol,
    ) -> Result<Atom, VakintError> {
        let rules = self.rules(catalogue, epsilon)?;
        let mut rotated = expression.clone();
        // One layer's exchanges never occur on its own right-hand sides
        // (checked by the tests), so a layer is one simultaneous rewrite;
        // consecutive layers are applied one after the other.
        let mut start = 0;
        for layer in &self.layers {
            let count = layer.exchanges.len();
            let layer = &rules[start..start + count];
            start += count;
            if !layer.iter().any(|(old, _)| contains(&rotated, old)) {
                continue;
            }
            let replacements = layer
                .iter()
                .map(|(old, new)| Replacement::new(old.to_pattern(), new.clone()))
                .collect::<Vec<_>>();
            rotated = rotated.replace_multiple(&replacements);
        }
        Ok(rotated)
    }
}

fn contains(expression: &Atom, part: &Atom) -> bool {
    let mut found = false;
    expression.visitor(&mut |view: AtomView| {
        if found {
            return false;
        }
        if view == part.as_view() {
            found = true;
            return false;
        }
        true
    });
    found
}

/// The masters directory of a five-loop package directory.
pub(crate) fn package_masters_directory(package: &Path) -> PathBuf {
    package.join(P_PACKAGE_MASTERS_DIRECTORY)
}

/// The basis change of `directory` ([`P_MASTERS_DIRECTORY_VARIABLE`] first,
/// else the loaded five-loop package's masters directory), parsed once per
/// directory. Without either: an error naming what is expected.
pub fn load(directory: Option<&Path>) -> Result<Arc<PMasterBasisData>, VakintError> {
    static LOADED: Mutex<BTreeMap<PathBuf, Arc<PMasterBasisData>>> = Mutex::new(BTreeMap::new());
    let directory = match std::env::var_os(P_MASTERS_DIRECTORY_VARIABLE) {
        Some(value) => PathBuf::from(value),
        None => directory.map(Path::to_owned).ok_or_else(|| {
            invalid(
                "(none)",
                format!(
                    "the pole-bounded master basis needs the five-loop package's masters directory \
                     <package>/{P_PACKAGE_MASTERS_DIRECTORY}/{P_MASTERS_MANIFEST_FILE} (no five-loop \
                     package directory is loaded) or {P_MASTERS_DIRECTORY_VARIABLE}"
                ),
            )
        })?,
    };
    let mut loaded = LOADED.lock().unwrap_or_else(PoisonError::into_inner);
    if let Some(data) = loaded.get(&directory) {
        return Ok(Arc::clone(data));
    }
    let data = Arc::new(PMasterBasisData::from_directory(&directory)?);
    loaded.insert(directory, Arc::clone(&data));
    Ok(data)
}

#[cfg(test)]
#[path = "master_basis_tests.rs"]
mod tests;
