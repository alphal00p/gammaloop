//! The hand-off of one R* remainder to vakint without expanding it.
//!
//! Today the adaptor of gammalooprs expands every product `P_i × D_j` of the
//! contracted numerator `P` and the Taylor coefficients `D_j` of the
//! denominator product into monomials, groups them per topology and hands
//! vakint one fully expanded numerator per topology. At six loops that
//! expansion is the whole cost (630 truncated to two momenta per vertex:
//! 413 kB of input, 1.0·10⁷ summands, 31 GB).
//!
//! A [`RemainderBoundary`] carries the same information factored:
//!
//! * `blocks`: closed scalars in loop-momentum dot products (contracted by
//!   idenso, no open index, no external momentum `p`), shared by every
//!   denominator power pattern of the remainder;
//! * `patterns`: one cut topology (`topo(prop(..))`, the line powers of the
//!   pattern) and one closing power per denominator power pattern, with its
//!   numerator `Σ_t blocks[t] · R_t`; the `R_t` are small (polynomials in the
//!   spectators and in a few loop dots);
//! * the closing line (`momentum`, `mass_squared`) and the cut routing, common
//!   to all patterns.
//!
//! [`RemainderBoundary::evaluate`] reads every pattern with the direct reader
//! of `p_lowering` (Horner over shared loop factors, numeric split, machine
//! integer slot substitution), absorbs all patterns with the same closing line
//! into one key plan and reduces it once. A pattern the reader declines either
//! fails with a typed error ([`BoundaryOptions::strict`]) or falls back to the
//! ordinary path, which expands that pattern's numerator.
//!
//! Block reuse: every block is read (Horner over its loop factors into the
//! exact canonical polynomial) once per distinct cut routing and reused by all
//! patterns; a pattern's polynomial is `Σ_t read(block_t) · read(R_t)` times
//! the cut factor, then classified, slot-substituted and emitted with its own
//! line powers (`VAKINT_P_BOUNDARY_BLOCKS=off` reads each assembled pattern
//! numerator instead). Chunked reading of one large block with budget
//! checkpoints is next; it does not change this type.
use symbolica::atom::{Atom, AtomCore};

use super::massive::{DirectClosureInput, LoweredClosure, MassiveClosure};
use crate::{Vakint, VakintError, VakintSettings};

/// One denominator power pattern of a remainder.
#[derive(Clone, Debug)]
pub struct BoundaryPattern {
    /// `topo(prop(id, edge(..), momentum, mass, power), ...)` of the cut graph
    /// with this pattern's line powers, in the cut loop momenta `k(1..)` and
    /// the closing head (as [`DirectClosureInput::topology`]).
    pub topology: Atom,
    /// Power of the closing massive line for this pattern.
    pub closure_power: i64,
    /// `(block index, coefficient)`: the pattern's numerator is
    /// `Σ blocks[index] · coefficient`. Coefficients may contain spectators
    /// (ε, `dot(p,p)`, masses, insertion coefficients) and loop dots.
    pub terms: Vec<(usize, Atom)>,
}

/// One R* remainder × massive closure, handed to vakint factored.
#[derive(Clone, Debug)]
pub struct RemainderBoundary {
    /// Free-form provenance (graph, colour branch, term); not used in values.
    pub label: String,
    /// The closing momentum `p(i)` (see [`MassiveClosure::momentum`]).
    pub momentum: Atom,
    /// The closing mass squared (see [`MassiveClosure::mass_squared`]).
    pub mass_squared: Atom,
    /// Every vacuum loop momentum and its image in the cut momenta, common to
    /// all patterns (as [`DirectClosureInput::routing`]).
    pub routing: Vec<(Atom, Atom)>,
    /// Shared closed scalar blocks in dot products.
    pub blocks: Vec<Atom>,
    /// The denominator power patterns.
    pub patterns: Vec<BoundaryPattern>,
}

/// How [`RemainderBoundary::evaluate`] treats a pattern the reader declines.
#[derive(Clone, Copy, Debug, Default)]
pub struct BoundaryOptions {
    /// Fail with an error naming the reason instead of expanding the pattern
    /// on the ordinary path. Also enabled by `VAKINT_P_DIRECT_STRICT=1`.
    pub strict: bool,
}

impl BoundaryOptions {
    fn strict(&self) -> bool {
        self.strict || std::env::var("VAKINT_P_DIRECT_STRICT").as_deref() == Ok("1")
    }
}

impl RemainderBoundary {
    /// A checked boundary: at least one pattern, every block index in range,
    /// no empty pattern.
    pub fn new(
        label: impl Into<String>,
        momentum: Atom,
        mass_squared: Atom,
        routing: Vec<(Atom, Atom)>,
        blocks: Vec<Atom>,
        patterns: Vec<BoundaryPattern>,
    ) -> Result<Self, VakintError> {
        if patterns.is_empty() {
            return Err(super::invalid("a remainder boundary without patterns"));
        }
        for (index, pattern) in patterns.iter().enumerate() {
            if pattern.terms.is_empty() {
                return Err(super::invalid(format!(
                    "pattern {index} of a remainder boundary has no terms"
                )));
            }
            if let Some((block, _)) = pattern
                .terms
                .iter()
                .find(|(block, _)| *block >= blocks.len())
            {
                return Err(super::invalid(format!(
                    "pattern {index} refers to block {block} of {}",
                    blocks.len()
                )));
            }
        }
        Ok(Self {
            label: label.into(),
            momentum,
            mass_squared,
            routing,
            blocks,
            patterns,
        })
    }

    /// The boundary of one existing direct closure input: one block (its
    /// numerator), one pattern (its topology) with coefficient one.
    pub fn from_direct(closure: &MassiveClosure, input: &DirectClosureInput) -> Self {
        Self {
            label: String::new(),
            momentum: closure.momentum.clone(),
            mass_squared: closure.mass_squared.clone(),
            routing: input.routing.clone(),
            blocks: vec![input.numerator.clone()],
            patterns: vec![BoundaryPattern {
                topology: input.topology.clone(),
                closure_power: closure.power,
                terms: vec![(0, Atom::num(1))],
            }],
        }
    }

    /// The closing line of `pattern`.
    pub fn closure(&self, pattern: &BoundaryPattern) -> MassiveClosure {
        MassiveClosure {
            momentum: self.momentum.clone(),
            mass_squared: self.mass_squared.clone(),
            power: pattern.closure_power,
        }
    }

    /// The direct closure input of `pattern` (its numerator assembled from
    /// the blocks, not expanded).
    pub fn direct_input(&self, pattern: &BoundaryPattern) -> DirectClosureInput {
        let numerator = Atom::add_many(
            pattern
                .terms
                .iter()
                .map(|(block, coefficient)| &self.blocks[*block] * coefficient),
        );
        DirectClosureInput {
            topology: pattern.topology.clone(),
            routing: self.routing.clone(),
            numerator,
        }
    }

    /// Every pattern as an ordinary closure input (`numerator · topo`,
    /// expanded): the reference the direct reading must reproduce.
    pub fn ordinary(&self) -> Vec<(MassiveClosure, Atom)> {
        self.patterns
            .iter()
            .map(|pattern| (self.closure(pattern), self.direct_input(pattern).ordinary()))
            .collect()
    }

    /// Read every pattern and absorb those with the same closing line into one
    /// key plan each. Returns the key plans and the patterns the reader
    /// declined (an error naming the reason instead when strict).
    pub fn lower(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        options: BoundaryOptions,
    ) -> Result<(Vec<LoweredClosure>, Vec<usize>), VakintError> {
        let (classes, stats) = self.lower_with_stats(vakint, settings, options)?;
        Ok((
            classes,
            stats.declined.iter().map(|(index, _)| *index).collect(),
        ))
    }

    /// [`Self::lower`] plus the block read counts `(reads, reuses)`.
    pub fn lower_counted(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        options: BoundaryOptions,
    ) -> Result<(Vec<LoweredClosure>, Vec<usize>, (usize, usize)), VakintError> {
        let (classes, stats) = self.lower_with_stats(vakint, settings, options)?;
        Ok((
            classes,
            stats.declined.iter().map(|(index, _)| *index).collect(),
            (stats.block_reads, stats.block_reuses),
        ))
    }

    /// [`Self::lower`] with its [`BoundaryStats`] (declines with their
    /// reasons, block reads / reuses, stage seconds). Every block is read once
    /// per distinct cut routing and reused by the other patterns
    /// (`VAKINT_P_BOUNDARY_BLOCKS=off`: each pattern's assembled numerator is
    /// read on its own).
    pub fn lower_with_stats(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        options: BoundaryOptions,
    ) -> Result<(Vec<LoweredClosure>, BoundaryStats), VakintError> {
        let mut cache = crate::rustred_evaluation::p_lowering::BlockReadCache::default();
        let mut direct = crate::rustred_evaluation::p_lowering::DirectStats::default();
        let mut stats = BoundaryStats {
            patterns: self.patterns.len(),
            ..BoundaryStats::default()
        };
        let mut classes = Vec::<LoweredClosure>::new();
        for index in 0..self.patterns.len() {
            // A checkpoint between patterns: the resource budget and the
            // execution limits of the settings apply to every pattern.
            settings.validate_execution()?;
            crate::execution::checkpoint().map_err(VakintError::from)?;
            match self.lower_one(vakint, settings, index, &mut cache, &mut direct)? {
                Ok(lowered) => {
                    stats.lowered += 1;
                    let mut lowered = Some(lowered);
                    for class in classes.iter_mut() {
                        match class.absorb(lowered.take().unwrap()) {
                            Ok(()) => break,
                            Err(back) => lowered = Some(back),
                        }
                    }
                    if let Some(lowered) = lowered {
                        classes.push(lowered);
                    }
                }
                Err(reason) if options.strict() => {
                    return Err(super::invalid(format!(
                        "remainder boundary {:?}: the direct reader declined pattern {index}: {reason} (strict mode)",
                        self.label
                    )));
                }
                Err(reason) => stats.declined.push((index, reason.to_string())),
            }
        }
        stats.key_plans = classes.len();
        stats.absorb_direct(&direct, &cache);
        Ok((classes, stats))
    }

    /// The folded RustRed input of one pattern alone (for check-mode pairs:
    /// compare its [`LoweredClosure::fingerprint`] with another lowering).
    /// `Ok(Err(reason))`: the direct reader declines.
    pub fn lower_pattern(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        index: usize,
    ) -> Result<Result<LoweredClosure, String>, VakintError> {
        if index >= self.patterns.len() {
            return Err(super::invalid(format!(
                "pattern {index} of {}",
                self.patterns.len()
            )));
        }
        let mut cache = crate::rustred_evaluation::p_lowering::BlockReadCache::default();
        let mut direct = crate::rustred_evaluation::p_lowering::DirectStats::default();
        Ok(self
            .lower_one(vakint, settings, index, &mut cache, &mut direct)?
            .map_err(str::to_string))
    }

    fn lower_one(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        index: usize,
        cache: &mut crate::rustred_evaluation::p_lowering::BlockReadCache,
        direct: &mut crate::rustred_evaluation::p_lowering::DirectStats,
    ) -> Result<Result<LoweredClosure, &'static str>, VakintError> {
        let pattern = &self.patterns[index];
        let closure = self.closure(pattern);
        let reuse = !matches!(
            std::env::var("VAKINT_P_BOUNDARY_BLOCKS").as_deref(),
            Ok("off" | "0")
        );
        if reuse {
            closure.lower_direct_blocks(
                vakint,
                settings,
                &pattern.topology,
                &self.routing,
                &self.blocks,
                &pattern.terms,
                cache,
                direct,
            )
        } else {
            closure.lower_direct_explained(vakint, settings, &self.direct_input(pattern), direct)
        }
    }

    /// The value of the remainder: Σ over patterns of the closed p-integral,
    /// one reduction per closing line.
    pub fn evaluate(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        options: BoundaryOptions,
    ) -> Result<Atom, VakintError> {
        Ok(self.evaluate_with_stats(vakint, settings, options)?.0)
    }

    /// [`Self::evaluate`] with its [`BoundaryStats`] (including the reduction
    /// and finalisation seconds of the key plans and the seconds of declined
    /// patterns on the ordinary path).
    pub fn evaluate_with_stats(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        options: BoundaryOptions,
    ) -> Result<(Atom, BoundaryStats), VakintError> {
        // Lane unknown-masters: strict refusal once, on the sum of the patterns.
        if let Some(_boundary) = crate::public_result_boundary() {
            let inner = crate::unknown_master_coefficients_carried(settings);
            let (value, stats) = self.evaluate_with_stats(vakint, &inner, options)?;
            vakint.refuse_unknown_master_coefficients(settings, &value)?;
            return Ok((value, stats));
        }
        settings.with_resource_budget(|| {
            let (classes, mut stats) = self.lower_with_stats(vakint, settings, options)?;
            let mut result = Atom::Zero;
            let reducing = std::time::Instant::now();
            for class in classes {
                result += class.evaluate(vakint, settings)?;
            }
            stats.reduce_finalize_seconds = reducing.elapsed().as_secs_f64();
            let ordinary = std::time::Instant::now();
            for (index, _) in &stats.declined {
                let pattern = &self.patterns[*index];
                result += MassiveClosure::evaluate_batch(
                    vakint,
                    settings,
                    &[(self.closure(pattern), self.direct_input(pattern).ordinary())],
                )?
                .pop()
                .unwrap_or_default();
            }
            stats.ordinary_seconds = ordinary.elapsed().as_secs_f64();
            Ok((result.expand(), stats))
        })
    }
}

/// What one boundary lowering / evaluation did (R9 of the producer design).
#[derive(Clone, Debug, Default)]
pub struct BoundaryStats {
    /// Patterns of the boundary, patterns lowered directly, key plans formed.
    pub patterns: usize,
    pub lowered: usize,
    pub key_plans: usize,
    /// Declined patterns and the reader's reason (ordinary path, non-strict).
    pub declined: Vec<(usize, String)>,
    /// Block reads performed / served from the cache.
    pub block_reads: usize,
    pub block_reuses: usize,
    /// Canonical polynomial terms, distinct loop monomials and lowered key
    /// entries over all patterns.
    pub canonical_terms: usize,
    pub distinct_monomials: usize,
    pub lowered_entries: usize,
    /// Stage seconds summed over patterns: topology match, numerator read,
    /// classification, slot substitution, channel finish.
    pub match_seconds: f64,
    pub read_seconds: f64,
    pub classify_seconds: f64,
    pub slot_seconds: f64,
    pub finish_seconds: f64,
    /// Reduction and finalisation of the key plans; ordinary path of declines.
    pub reduce_finalize_seconds: f64,
    pub ordinary_seconds: f64,
}

impl BoundaryStats {
    fn absorb_direct(
        &mut self,
        direct: &crate::rustred_evaluation::p_lowering::DirectStats,
        cache: &crate::rustred_evaluation::p_lowering::BlockReadCache,
    ) {
        self.block_reads = cache.reads;
        self.block_reuses = cache.hits;
        self.canonical_terms = direct.canonical_terms;
        self.distinct_monomials = direct.distinct_monomials;
        self.lowered_entries = direct.lowered;
        self.match_seconds = direct.match_seconds;
        self.read_seconds = direct.read_seconds;
        self.classify_seconds = direct.classify_seconds;
        self.slot_seconds = direct.scatter_seconds;
        self.finish_seconds = direct.finish_seconds;
    }
}

/// One `pint6.remainder-boundary.v1` record read back ([`read_boundary_record`]).
#[derive(Clone, Debug)]
pub struct BoundaryRecord {
    pub boundary: RemainderBoundary,
    /// `epsilon_terms` and `p_loops` of the header, when present.
    pub epsilon_terms: Option<i64>,
    pub p_loops: Option<i64>,
    /// Per pattern, the v1 closure sequence it replaces (`ref_closure`).
    pub ref_closures: Vec<Option<String>>,
}

/// Schema line of a boundary header.
pub const BOUNDARY_SCHEMA: &str = "pint6.remainder-boundary.v1";
/// Schema line of a block file.
pub const BOUNDARY_BLOCK_SCHEMA: &str = "pint6.boundary-block.v1";
/// Summands added per chunk while a block file is streamed.
const BLOCK_CHUNK: usize = 65_536;

fn parse_atom(text: &str) -> Result<Atom, VakintError> {
    crate::utils::vakint_macros::vk_parse!(text)
        .map_err(|error| super::invalid(format!("cannot parse {:.80}: {error}", text)))
}

/// What a header declares about one block file; every declared item is
/// checked against the file ([`read_boundary_block_checked`]).
#[derive(Clone, Debug, Default)]
pub struct BlockExpectation {
    /// The block id the header references at this position: the file's own
    /// `block` line must name it.
    pub id: Option<usize>,
    /// `summands=` of the header's block line.
    pub summands: Option<usize>,
    /// `text_bytes=` of the header's block line: the bytes of the decompressed
    /// file (its three header lines and the summand lines, with newlines).
    pub text_bytes: Option<u64>,
    /// sha256 (lower-case hex) of the file as stored, from the record's
    /// `SHA256SUMS`.
    pub sha256: Option<String>,
}

/// Stream one block file (`pint6.boundary-block.v1`, zstd unless the name
/// does not end in `.zst`): a summand per line, added in chunks of 65,536;
/// the summand count of the file header must match the lines read.
pub fn read_boundary_block(path: &std::path::Path) -> Result<Atom, VakintError> {
    read_boundary_block_checked(path, &BlockExpectation::default())
}

/// The sha256 of a file, lower-case hex.
fn file_sha256(path: &std::path::Path) -> Result<String, VakintError> {
    use sha2::Digest;
    use std::io::Read;
    let io = |error: std::io::Error| super::invalid(format!("{}: {error}", path.display()));
    let mut file = std::fs::File::open(path).map_err(io)?;
    let mut hasher = sha2::Sha256::new();
    let mut buffer = vec![0u8; 1 << 20];
    loop {
        let n = file.read(&mut buffer).map_err(io)?;
        if n == 0 {
            break;
        }
        hasher.update(&buffer[..n]);
    }
    Ok(hasher
        .finalize()
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect())
}

/// [`read_boundary_block`] with the header's declarations: a block file
/// whose digest, embedded block id, summand count or text size differs from
/// what the header declares is an error naming the file and the item.
pub fn read_boundary_block_checked(
    path: &std::path::Path,
    expect: &BlockExpectation,
) -> Result<Atom, VakintError> {
    use std::io::BufRead;
    let io = |error: std::io::Error| super::invalid(format!("{}: {error}", path.display()));
    if let Some(digest) = &expect.sha256 {
        let found = file_sha256(path)?;
        if &found != digest {
            return Err(super::invalid(format!(
                "{}: sha256 {found} differs from SHA256SUMS {digest}",
                path.display()
            )));
        }
    }
    let mut child = None;
    let reader: Box<dyn BufRead> = if path.to_string_lossy().ends_with(".zst") {
        let mut process = std::process::Command::new("zstd")
            .args(["-qdc"])
            .arg(path)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(io)?;
        let stdout = process.stdout.take().unwrap();
        child = Some(process);
        Box::new(std::io::BufReader::with_capacity(1 << 20, stdout))
    } else {
        Box::new(std::io::BufReader::with_capacity(
            1 << 20,
            std::fs::File::open(path).map_err(io)?,
        ))
    };
    let mut lines = reader.lines();
    let mut header = |key: &str| -> Result<String, VakintError> {
        let line = lines
            .next()
            .ok_or_else(|| super::invalid(format!("{}: no {key} line", path.display())))?
            .map_err(io)?;
        Ok(line)
    };
    let schema = header("schema")?;
    if schema.trim() != format!("schema {BOUNDARY_BLOCK_SCHEMA}") {
        return Err(super::invalid(format!(
            "{}: unknown block schema {schema:?}",
            path.display()
        )));
    }
    let block_line = header("block")?;
    let id: usize = block_line
        .strip_prefix("block ")
        .and_then(|id| id.trim().parse().ok())
        .ok_or_else(|| {
            super::invalid(format!("{}: bad block line {block_line:?}", path.display()))
        })?;
    if let Some(expected) = expect.id
        && id != expected
    {
        return Err(super::invalid(format!(
            "{}: the file holds block {id}, the header references block {expected} here",
            path.display()
        )));
    }
    let summands_line = header("summands")?;
    let summands: usize = summands_line
        .strip_prefix("summands ")
        .and_then(|count| count.trim().parse().ok())
        .ok_or_else(|| super::invalid(format!("{}: bad summands line", path.display())))?;
    if let Some(expected) = expect.summands
        && summands != expected
    {
        return Err(super::invalid(format!(
            "{}: the file declares {summands} summands, the header {expected}",
            path.display()
        )));
    }
    let mut partial = Vec::new();
    let mut chunk = Vec::with_capacity(BLOCK_CHUNK.min(summands));
    let mut read = 0usize;
    // text_bytes = the decompressed file: the three header lines and the
    // summand lines, each with its newline.
    let mut text_bytes = (schema.len() + block_line.len() + summands_line.len() + 3) as u64;
    for line in lines {
        let line = line.map_err(io)?;
        text_bytes += line.len() as u64 + 1;
        if line.is_empty() {
            continue;
        }
        chunk.push(parse_atom(&line)?);
        read += 1;
        if chunk.len() == BLOCK_CHUNK {
            partial.push(Atom::add_many(chunk.iter()));
            chunk.clear();
            crate::execution::checkpoint().map_err(VakintError::from)?;
        }
    }
    if !chunk.is_empty() {
        partial.push(Atom::add_many(chunk.iter()));
    }
    if let Some(mut child) = child {
        let status = child.wait().map_err(io)?;
        if !status.success() {
            return Err(super::invalid(format!("{}: zstd failed", path.display())));
        }
    }
    if read != summands {
        return Err(super::invalid(format!(
            "{}: {read} summands read, the header says {summands}",
            path.display()
        )));
    }
    if let Some(expected) = expect.text_bytes
        && text_bytes != expected
    {
        return Err(super::invalid(format!(
            "{}: {text_bytes} bytes of decompressed text read, the header declares text_bytes={expected}",
            path.display()
        )));
    }
    Ok(Atom::add_many(partial.iter()))
}

/// `SHA256SUMS` beside a header (`<hex>  <file name>` lines), if present.
fn read_sums(
    directory: &std::path::Path,
) -> Result<Option<std::collections::BTreeMap<String, String>>, VakintError> {
    let path = directory.join("SHA256SUMS");
    if !path.is_file() {
        return Ok(None);
    }
    let text = std::fs::read_to_string(&path)
        .map_err(|error| super::invalid(format!("{}: {error}", path.display())))?;
    let mut sums = std::collections::BTreeMap::new();
    for line in text.lines().filter(|line| !line.trim().is_empty()) {
        let (digest, name) = line
            .split_once(char::is_whitespace)
            .ok_or_else(|| super::invalid(format!("{}: bad line {line:?}", path.display())))?;
        let name = name.trim_start().trim_start_matches('*').to_string();
        if sums
            .insert(name.clone(), digest.to_ascii_lowercase())
            .is_some()
        {
            return Err(super::invalid(format!(
                "{}: {name} listed twice",
                path.display()
            )));
        }
    }
    Ok(Some(sums))
}

/// Read a `pint6.remainder-boundary.v1` header and its block files (same
/// directory), streaming every block ([`read_boundary_block_checked`]). Only
/// the keys vakint needs are interpreted: closure_momentum,
/// closure_mass_squared, epsilon_terms, p_loops, label, route, block /
/// block_atom, pattern (closure_power, ref_closure), pattern_topology,
/// pattern_term; the others are informational.
///
/// Integrity (a truncated or mixed-up dump must never give a value): the
/// header must end with its `end` line; every declared cardinality must equal
/// what is read (`blocks`, `patterns`, `routing`, `classes`, a block's
/// `summands=` and `text_bytes=`, a pattern's `terms=` and `lines=` against its
/// `pattern_prop` rows); no block, pattern or pattern_topology is declared
/// twice; every pattern has one topology and every pattern_term names a
/// declared pattern and an existing block; the classes, when declared,
/// partition the patterns; a block file names the block the header references
/// for it; with a `SHA256SUMS` beside the header, the header and every block
/// file read must be listed and match.
pub fn read_boundary_record(path: &std::path::Path) -> Result<BoundaryRecord, VakintError> {
    // Vakint's symbols (with their attributes) before any parse.
    crate::Vakint::initialize_vakint_symbols();
    let bad = |what: String| super::invalid(format!("{}: {what}", path.display()));
    let directory = path.parent().unwrap_or(std::path::Path::new("."));
    let sums = read_sums(directory)?;
    let name_of = |file: &std::path::Path| {
        file.file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string()
    };
    if let Some(sums) = &sums {
        let name = name_of(path);
        let digest = sums
            .get(&name)
            .ok_or_else(|| bad("not listed in SHA256SUMS".to_string()))?;
        let found = file_sha256(path)?;
        if &found != digest {
            return Err(bad(format!(
                "sha256 {found} differs from SHA256SUMS {digest}"
            )));
        }
    }
    let text = std::fs::read_to_string(path)
        .map_err(|error| super::invalid(format!("{}: {error}", path.display())))?;
    let mut lines = text.lines();
    if lines.next().map(str::trim) != Some(format!("schema {BOUNDARY_SCHEMA}").as_str()) {
        return Err(super::invalid(format!(
            "{}: not a {BOUNDARY_SCHEMA} header",
            path.display()
        )));
    }
    let mut label = String::new();
    let (mut momentum, mut mass) = (None, None);
    let (mut epsilon_terms, mut p_loops) = (None, None);
    let mut routing = Vec::new();
    let mut block_files = std::collections::BTreeMap::<usize, (String, BlockExpectation)>::new();
    let mut block_atoms = std::collections::BTreeMap::<usize, Atom>::new();
    let mut powers = std::collections::BTreeMap::<usize, (i64, Option<String>)>::new();
    let mut declared_terms = std::collections::BTreeMap::<usize, usize>::new();
    let mut declared_lines = std::collections::BTreeMap::<usize, usize>::new();
    let mut props = std::collections::BTreeMap::<usize, usize>::new();
    let mut topologies = std::collections::BTreeMap::<usize, Atom>::new();
    let mut terms = std::collections::BTreeMap::<usize, Vec<(usize, Atom)>>::new();
    let mut class_patterns = Vec::<Vec<usize>>::new();
    let (mut n_blocks, mut n_patterns, mut n_routing, mut n_classes) = (None, None, None, None);
    let mut ended = false;
    let field = |fields: &[&str], name: &str| -> Option<String> {
        fields
            .iter()
            .find_map(|field| field.strip_prefix(name)?.strip_prefix('='))
            .map(str::to_string)
    };
    let index = |text: &str| -> Result<usize, VakintError> {
        text.trim()
            .parse()
            .map_err(|_| super::invalid(format!("{}: bad index {text:?}", path.display())))
    };
    let count = |name: &str, text: &str| -> Result<usize, VakintError> {
        text.trim()
            .parse()
            .map_err(|_| super::invalid(format!("{}: bad {name} count {text:?}", path.display())))
    };
    // Key policy (Astra acceptance 2026-11-18): a single-valued field occurs
    // at most once (a repeat is an error even when the values agree); an
    // indexed row occurs at most once per index; a field name occurs at most
    // once inside a row; a row of the wrong arity and an unknown key are
    // errors. INFORMATIONAL single-valued keys are accepted and not
    // interpreted.
    const SINGLE: &[&str] = &[
        "label",
        "closure_momentum",
        "closure_mass_squared",
        "epsilon_terms",
        "p_loops",
        "blocks",
        "patterns",
        "routing",
        "classes",
        "end",
    ];
    const INFORMATIONAL: &[&str] = &[
        "sequence",
        "context",
        "integration_id",
        "vacuum_loops",
        "evaluated",
        "taylor_degree",
        "colour_basis",
        "colour_branch",
        "colour_coefficient",
        "colour_coefficient_symbolic",
        "normalisation",
        "counts",
        "zero_patterns",
    ];
    let mut seen_single = std::collections::BTreeSet::<&str>::new();
    let mut seen_rows = std::collections::BTreeSet::<(String, usize, Option<usize>)>::new();
    for line in lines {
        if ended {
            if line.trim().is_empty() {
                continue;
            }
            return Err(bad(format!("a line after the end line: {:.80}", line)));
        }
        let (key, rest) = line.split_once([' ', '\t']).unwrap_or((line.trim(), ""));
        let fields = rest.split('\t').collect::<Vec<_>>();
        if let Some(single) = SINGLE.iter().chain(INFORMATIONAL).find(|k| **k == key)
            && !seen_single.insert(*single)
        {
            return Err(bad(format!("field {key} occurs more than once")));
        }
        {
            let mut names = std::collections::BTreeSet::new();
            for name in fields
                .iter()
                .skip(1)
                .filter_map(|f| f.split_once('=').map(|(n, _)| n))
            {
                if !name.is_empty()
                    && name.chars().all(|c| c.is_ascii_alphanumeric() || c == '_')
                    && !names.insert(name)
                {
                    return Err(bad(format!(
                        "{key} row: field {name} occurs more than once"
                    )));
                }
            }
        }
        let arity = |n: usize| -> Result<(), VakintError> {
            if fields.len() == n {
                Ok(())
            } else {
                Err(bad(format!(
                    "{key} row with {} fields, {n} expected: {:.80}",
                    fields.len(),
                    line
                )))
            }
        };
        let mut row = |second: Option<usize>| -> Result<(), VakintError> {
            let first = index(fields[0])?;
            if !seen_rows.insert((key.to_string(), first, second)) {
                return Err(bad(format!(
                    "{key} {first}{} occurs more than once",
                    second.map_or(String::new(), |s| format!(" {s}"))
                )));
            }
            Ok(())
        };
        match key {
            "end" => ended = true,
            "label" => label = rest.to_string(),
            "closure_momentum" => momentum = Some(parse_atom(rest)?),
            "closure_mass_squared" => mass = Some(parse_atom(rest)?),
            "epsilon_terms" => epsilon_terms = Some(count("epsilon_terms", rest)? as i64),
            "p_loops" => p_loops = Some(count("p_loops", rest)? as i64),
            "blocks" => n_blocks = Some(count("blocks", rest)?),
            "patterns" => n_patterns = Some(count("patterns", rest)?),
            "routing" => n_routing = Some(count("routing", rest)?),
            "classes" => n_classes = Some(count("classes", rest)?),
            k if INFORMATIONAL.contains(&k) => {}
            "route" => {
                arity(2)?;
                routing.push((parse_atom(fields[0])?, parse_atom(fields[1])?))
            }
            "block" => {
                row(None)?;
                let g = index(fields[0])?;
                let file = field(&fields, "file").unwrap_or_else(|| "inline".into());
                let expect = BlockExpectation {
                    id: Some(g),
                    summands: field(&fields, "summands")
                        .map(|n| count("summands", &n))
                        .transpose()?,
                    text_bytes: field(&fields, "text_bytes")
                        .map(|n| count("text_bytes", &n).map(|n| n as u64))
                        .transpose()?,
                    sha256: None,
                };
                block_files.insert(g, (file, expect));
            }
            "block_atom" => {
                arity(2)?;
                row(None)?;
                block_atoms.insert(index(fields[0])?, parse_atom(fields[1])?);
            }
            "pattern" => {
                row(None)?;
                let pi = index(fields[0])?;
                let power = field(&fields, "closure_power")
                    .and_then(|power| power.parse().ok())
                    .ok_or_else(|| bad(format!("pattern {pi} without closure_power")))?;
                let reference = field(&fields, "ref_closure").filter(|seq| seq != "none");
                powers.insert(pi, (power, reference));
                if let Some(n) = field(&fields, "terms") {
                    declared_terms.insert(pi, count("terms", &n)?);
                }
                if let Some(n) = field(&fields, "lines") {
                    declared_lines.insert(pi, count("lines", &n)?);
                }
            }
            "pattern_topology" => {
                arity(2)?;
                row(None)?;
                topologies.insert(index(fields[0])?, parse_atom(fields[1])?);
            }
            "pattern_integral" | "vacuum" | "class" => row(None)?,
            "pattern_prop" => {
                if fields.len() < 2 {
                    return Err(bad(format!("{key} row too short: {:.80}", line)));
                }
                row(Some(index(fields[1])?))?;
                *props.entry(index(fields[0])?).or_default() += 1;
            }
            "pattern_term" => {
                // No uniqueness: a pattern_term row has no row index; several
                // terms of one pattern on the same block are additive (the
                // writer emits them); the pattern's terms= count guards the
                // row set.
                arity(3)?;
                terms
                    .entry(index(fields[0])?)
                    .or_default()
                    .push((index(fields[1])?, parse_atom(fields[2])?));
            }
            other => return Err(bad(format!("unknown field {other:?}"))),
        }
        if key == "class" {
            let members = field(&fields, "patterns").unwrap_or_default();
            let members = members
                .split(',')
                .filter(|m| !m.trim().is_empty())
                .map(index)
                .collect::<Result<Vec<_>, _>>()?;
            class_patterns.push(members);
        }
    }
    if !ended {
        return Err(bad("no end line: the header is truncated".into()));
    }
    let check = |name: &str, declared: Option<usize>, read: usize| -> Result<(), VakintError> {
        match declared {
            Some(declared) if declared != read => Err(bad(format!(
                "{read} {name} read, the header declares {declared}"
            ))),
            _ => Ok(()),
        }
    };
    check("blocks", n_blocks, block_files.len())?;
    check("patterns", n_patterns, powers.len())?;
    check("routes", n_routing, routing.len())?;
    check("classes", n_classes, class_patterns.len())?;
    if powers.is_empty() {
        return Err(bad("no pattern".into()));
    }
    if let Some(pi) = topologies.keys().find(|pi| !powers.contains_key(pi)) {
        return Err(bad(format!(
            "pattern_topology {pi} of an undeclared pattern"
        )));
    }
    if let Some(pi) = terms.keys().find(|pi| !powers.contains_key(pi)) {
        return Err(bad(format!("pattern_term of an undeclared pattern {pi}")));
    }
    for (pi, declared) in &declared_terms {
        let read = terms.get(pi).map_or(0, Vec::len);
        if read != *declared {
            return Err(bad(format!(
                "pattern {pi}: {read} pattern_term rows read, the pattern declares terms={declared}"
            )));
        }
    }
    if !props.is_empty() {
        for (pi, declared) in &declared_lines {
            let read = props.get(pi).copied().unwrap_or(0);
            if read != *declared {
                return Err(bad(format!(
                    "pattern {pi}: {read} pattern_prop rows read, the pattern declares lines={declared}"
                )));
            }
        }
    }
    if !class_patterns.is_empty() {
        let mut seen = std::collections::BTreeSet::new();
        for pi in class_patterns.iter().flatten() {
            if !seen.insert(*pi) {
                return Err(bad(format!("pattern {pi} is in two classes")));
            }
        }
        if seen.iter().ne(powers.keys()) {
            return Err(bad("the classes do not partition the patterns".into()));
        }
    }
    let mut blocks = Vec::with_capacity(block_files.len());
    for (position, (g, (file, expect))) in block_files.iter().enumerate() {
        if *g != position {
            return Err(bad("blocks are not numbered 0..".into()));
        }
        let block = if file == "inline" {
            block_atoms
                .remove(g)
                .ok_or_else(|| bad(format!("inline block {g} without block_atom")))?
        } else {
            let block_path = directory.join(file);
            let mut expect = expect.clone();
            if let Some(sums) = &sums {
                expect.sha256 = Some(
                    sums.get(&name_of(&block_path))
                        .ok_or_else(|| bad(format!("block file {file} not listed in SHA256SUMS")))?
                        .clone(),
                );
            }
            read_boundary_block_checked(&block_path, &expect)?
        };
        blocks.push(block);
    }
    if let Some(g) = block_atoms.keys().next() {
        return Err(bad(format!("block_atom {g} of no inline block")));
    }
    if let Some((pi, g)) = terms
        .iter()
        .flat_map(|(pi, rows)| rows.iter().map(move |(g, _)| (pi, g)))
        .find(|(_, g)| **g >= blocks.len())
    {
        return Err(bad(format!(
            "pattern {pi}: a pattern_term names block {g} of {}",
            blocks.len()
        )));
    }
    let mut patterns = Vec::with_capacity(powers.len());
    let mut ref_closures = Vec::with_capacity(powers.len());
    for (pi, (closure_power, reference)) in powers {
        let topology = topologies
            .remove(&pi)
            .ok_or_else(|| bad(format!("pattern {pi} without pattern_topology")))?;
        patterns.push(BoundaryPattern {
            topology,
            closure_power,
            terms: terms.remove(&pi).unwrap_or_default(),
        });
        ref_closures.push(reference);
    }
    let missing = |name: &str| super::invalid(format!("{}: no {name}", path.display()));
    let boundary = RemainderBoundary::new(
        label,
        momentum.ok_or_else(|| missing("closure_momentum"))?,
        mass.ok_or_else(|| missing("closure_mass_squared"))?,
        routing,
        blocks,
        patterns,
    )?;
    Ok(BoundaryRecord {
        boundary,
        epsilon_terms,
        p_loops,
        ref_closures,
    })
}

/// Write `boundary` as a `pint6.remainder-boundary.v1` header (`path`) with
/// blocks inline when `inline(block)` holds and otherwise one block file per
/// block beside it (zstd through the system binary when `compress`). The
/// reference writer for [`read_boundary_record`]; the producer writes the full
/// header of its design (section 6b).
pub fn write_boundary_record(
    boundary: &RemainderBoundary,
    path: &std::path::Path,
    epsilon_terms: i64,
    compress: bool,
    inline: impl Fn(&Atom) -> bool,
) -> Result<(), VakintError> {
    use std::io::Write;
    let io = |error: std::io::Error| super::invalid(format!("{}: {error}", path.display()));
    let one_line = |atom: &Atom| -> Result<String, VakintError> {
        let text = atom.to_canonical_string();
        if text.contains('\n') {
            return Err(super::invalid("an atom spans several lines"));
        }
        Ok(text)
    };
    let stem = path
        .file_name()
        .unwrap()
        .to_string_lossy()
        .trim_end_matches(".txt")
        .to_string();
    let mut header = String::new();
    header.push_str(&format!("schema {BOUNDARY_SCHEMA}\n"));
    header.push_str(&format!("label {}\n", boundary.label));
    header.push_str(&format!("epsilon_terms {epsilon_terms}\n"));
    header.push_str(&format!(
        "closure_momentum {}\n",
        one_line(&boundary.momentum)?
    ));
    header.push_str(&format!(
        "closure_mass_squared {}\n",
        one_line(&boundary.mass_squared)?
    ));
    header.push_str(&format!("routing {}\n", boundary.routing.len()));
    for (momentum, image) in &boundary.routing {
        header.push_str(&format!(
            "route {}\t{}\n",
            one_line(momentum)?,
            one_line(image)?
        ));
    }
    header.push_str(&format!("blocks {}\n", boundary.blocks.len()));
    for (g, block) in boundary.blocks.iter().enumerate() {
        let summands = match block.as_view() {
            symbolica::atom::AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
            other => vec![other],
        };
        if inline(block) {
            header.push_str(&format!(
                "block {g}\tsummands={}\tfile=inline\n",
                summands.len()
            ));
            header.push_str(&format!("block_atom {g}\t{}\n", one_line(block)?));
            continue;
        }
        let name = format!(
            "{stem}.block-{g:03}.txt{}",
            if compress { ".zst" } else { "" }
        );
        header.push_str(&format!(
            "block {g}\tsummands={}\tfile={name}\n",
            summands.len()
        ));
        let target = path.with_file_name(&name);
        let mut body = Vec::<u8>::new();
        writeln!(
            body,
            "schema {BOUNDARY_BLOCK_SCHEMA}\nblock {g}\nsummands {}",
            summands.len()
        )
        .map_err(io)?;
        for summand in summands {
            writeln!(body, "{}", one_line(&summand.to_owned())?).map_err(io)?;
        }
        if compress {
            let mut process = std::process::Command::new("zstd")
                .args(["-q", "-3", "-T1", "-f", "-o"])
                .arg(&target)
                .stdin(std::process::Stdio::piped())
                .spawn()
                .map_err(io)?;
            process.stdin.take().unwrap().write_all(&body).map_err(io)?;
            if !process.wait().map_err(io)?.success() {
                return Err(super::invalid("zstd failed"));
            }
        } else {
            std::fs::write(&target, body).map_err(io)?;
        }
    }
    header.push_str(&format!("patterns {}\n", boundary.patterns.len()));
    for (pi, pattern) in boundary.patterns.iter().enumerate() {
        header.push_str(&format!(
            "pattern {pi}\tclosure_power={}\tterms={}\tref_closure=none\n",
            pattern.closure_power,
            pattern.terms.len()
        ));
        header.push_str(&format!(
            "pattern_topology {pi}\t{}\n",
            one_line(&pattern.topology)?
        ));
        for (g, coefficient) in &pattern.terms {
            header.push_str(&format!(
                "pattern_term {pi}\t{g}\t{}\n",
                one_line(coefficient)?
            ));
        }
    }
    header.push_str("end\n");
    std::fs::write(path, header).map_err(io)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::vakint_macros::vk_parse;

    fn pattern(topology: &str, terms: Vec<(usize, &str)>) -> BoundaryPattern {
        BoundaryPattern {
            topology: vk_parse!(topology).unwrap(),
            closure_power: 1,
            terms: terms
                .into_iter()
                .map(|(block, coefficient)| (block, vk_parse!(coefficient).unwrap()))
                .collect(),
        }
    }

    #[test]
    fn remainder_boundary_constructor_checks_patterns_and_blocks() {
        // Vakint's symbols (with their tensor attributes) before any parse:
        // a bare vk_parse!("p(1)") would otherwise define vakint::p first.
        crate::Vakint::initialize_vakint_symbols();
        let momentum = vk_parse!("p(1)").unwrap();
        let mass = vk_parse!("mUV^2").unwrap();
        let blocks = vec![
            vk_parse!("dot(k(1),k(1))").unwrap(),
            vk_parse!("1").unwrap(),
        ];
        let topology = "topo(prop(1,edge(1,2),k(1),0,1))";
        let ok = RemainderBoundary::new(
            "t",
            momentum.clone(),
            mass.clone(),
            Vec::new(),
            blocks.clone(),
            vec![pattern(topology, vec![(0, "1"), (1, "ε")])],
        )
        .unwrap();
        let input = ok.direct_input(&ok.patterns[0]);
        assert_eq!(
            input.numerator,
            vk_parse!("dot(k(1),k(1))+ε").unwrap(),
            "numerator = Σ blocks[t] · R_t"
        );
        assert_eq!(ok.closure(&ok.patterns[0]).power, 1);
        assert!(
            RemainderBoundary::new(
                "t",
                momentum.clone(),
                mass.clone(),
                Vec::new(),
                blocks.clone(),
                vec![]
            )
            .is_err()
        );
        assert!(
            RemainderBoundary::new(
                "t",
                momentum.clone(),
                mass.clone(),
                Vec::new(),
                blocks.clone(),
                vec![pattern(topology, vec![(2, "1")])],
            )
            .is_err()
        );
        assert!(
            RemainderBoundary::new(
                "t",
                momentum,
                mass,
                Vec::new(),
                blocks,
                vec![pattern(topology, vec![])],
            )
            .is_err()
        );
    }

    #[test]
    fn remainder_boundary_record_round_trips() {
        crate::Vakint::initialize_vakint_symbols();
        let directory = std::env::temp_dir().join(format!(
            "vakint-boundary-record-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&directory).unwrap();
        let compress = std::process::Command::new("zstd")
            .arg("--version")
            .output()
            .is_ok_and(|output| output.status.success());
        let boundary = RemainderBoundary::new(
            "graph=t\tterm=1",
            vk_parse!("p(1)").unwrap(),
            vk_parse!("mUV^2").unwrap(),
            vec![
                (vk_parse!("k(1)").unwrap(), vk_parse!("k(1)").unwrap()),
                (vk_parse!("k(2)").unwrap(), vk_parse!("k(1)-k(2)").unwrap()),
                (vk_parse!("p(1)").unwrap(), vk_parse!("p(1)").unwrap()),
            ],
            vec![
                vk_parse!("dot(k(1),k(2))").unwrap(),
                vk_parse!(
                    "(4-2*ep)*dot(k(1),k(1))^2+3/7*dot(k(2),k(2))*dot(k(1),k(2))+dot(k(2),k(2))"
                )
                .unwrap(),
            ],
            vec![
                pattern(
                    "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,1))",
                    vec![(0, "1"), (1, "dot(p(1),p(1))")],
                ),
                pattern(
                    "topo(prop(1,edge(1,2),k(1),0,2)*prop(2,edge(2,1),k(1)-p(1),0,1))",
                    vec![(1, "1/(4-2*ep)")],
                ),
            ],
        )
        .unwrap();
        let header = directory.join("boundary-1-000000.txt");
        write_boundary_record(&boundary, &header, 5, compress, |block| {
            block.as_view().get_byte_size() < 40
        })
        .unwrap();
        let record = read_boundary_record(&header).unwrap();
        let back = &record.boundary;
        assert_eq!(record.epsilon_terms, Some(5));
        assert_eq!(back.label, boundary.label);
        assert_eq!(back.momentum, boundary.momentum);
        assert_eq!(back.mass_squared, boundary.mass_squared);
        assert_eq!(back.routing, boundary.routing);
        assert_eq!(back.blocks, boundary.blocks);
        assert_eq!(back.patterns.len(), 2);
        for (left, right) in back.patterns.iter().zip(&boundary.patterns) {
            assert_eq!(left.topology, right.topology);
            assert_eq!(left.closure_power, right.closure_power);
            assert_eq!(left.terms, right.terms);
        }
        assert_eq!(record.ref_closures, vec![None, None]);
        // A short block file is refused.
        let block_file = std::fs::read_dir(&directory)
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .find(|path| path.to_string_lossy().contains(".block-"))
            .expect("one block file");
        if !block_file.to_string_lossy().ends_with(".zst") {
            let text = std::fs::read_to_string(&block_file).unwrap();
            let short = text
                .lines()
                .take(text.lines().count() - 1)
                .collect::<Vec<_>>()
                .join("\n");
            std::fs::write(&block_file, short).unwrap();
            assert!(read_boundary_record(&header).is_err());
        }
        let _ = std::fs::remove_dir_all(&directory);
    }

    /// vk5l-physrun (Astra review 2026-11-18): a damaged dump never gives a
    /// value. Two-block bubble written by the reference writer; every
    /// damage (a missing pattern_term row, interchanged block files, a
    /// truncated header or block file, a duplicated or extra row, a count
    /// that disagrees, a digest that disagrees) is refused with an error
    /// naming the record; the intact record (with and without SHA256SUMS)
    /// reads back.
    #[test]
    fn damaged_boundary_records_are_refused() {
        crate::Vakint::initialize_vakint_symbols();
        let directory =
            std::env::temp_dir().join(format!("vakint-boundary-integrity-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&directory);
        std::fs::create_dir_all(&directory).unwrap();
        let boundary = RemainderBoundary::new(
            "integrity",
            vk_parse!("p(1)").unwrap(),
            vk_parse!("1").unwrap(),
            vec![
                (vk_parse!("k(1)").unwrap(), vk_parse!("k(1)").unwrap()),
                (vk_parse!("p(1)").unwrap(), vk_parse!("p(1)").unwrap()),
            ],
            vec![vk_parse!("1").unwrap(), vk_parse!("3").unwrap()],
            vec![pattern(
                "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,1))",
                vec![(0, "ε^2"), (1, "2*ε^2")],
            )],
        )
        .unwrap();
        let header = directory.join("boundary-1-000000.txt");
        let block = |g: usize| directory.join(format!("boundary-1-000000.block-{g:03}.txt"));
        let write = || {
            for entry in std::fs::read_dir(&directory).unwrap() {
                std::fs::remove_file(entry.unwrap().path()).unwrap();
            }
            write_boundary_record(&boundary, &header, 3, false, |_| false).unwrap();
        };
        let text = || std::fs::read_to_string(&header).unwrap();
        let error = |case: &str| {
            let error = read_boundary_record(&header)
                .err()
                .unwrap_or_else(|| panic!("{case}: a damaged record was read"));
            let message = format!("{error:?}");
            assert!(message.contains("boundary-1-000000"), "{case}: {message}");
            message
        };
        write();
        let intact = read_boundary_record(&header).unwrap();
        assert_eq!(intact.boundary.patterns[0].terms.len(), 2);
        // (a) the last pattern_term row deleted, terms=2 kept.
        write();
        let full = text();
        let rows: Vec<&str> = full.lines().collect();
        let last = rows
            .iter()
            .rposition(|l| l.starts_with("pattern_term"))
            .unwrap();
        let kept: Vec<&str> = rows
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != last)
            .map(|(_, l)| *l)
            .collect();
        std::fs::write(&header, kept.join("\n") + "\n").unwrap();
        assert!(error("missing term").contains("pattern_term rows"));
        // (b) the two block files interchanged (each keeps its embedded id).
        write();
        let (b0, b1) = (
            std::fs::read(block(0)).unwrap(),
            std::fs::read(block(1)).unwrap(),
        );
        std::fs::write(block(0), &b1).unwrap();
        std::fs::write(block(1), &b0).unwrap();
        assert!(error("swapped blocks").contains("references block"));
        // (c) a truncated header: cut before its end line, and inside the rows.
        write();
        let full = text();
        std::fs::write(&header, full.trim_end().trim_end_matches("end")).unwrap();
        assert!(error("no end line").contains("truncated"));
        std::fs::write(&header, &full[..full.find("pattern_term").unwrap()]).unwrap();
        error("truncated in the rows");
        // (d) a truncated block file.
        write();
        let body = std::fs::read_to_string(block(1)).unwrap();
        std::fs::write(
            block(1),
            body.lines().take(3).collect::<Vec<_>>().join("\n"),
        )
        .unwrap();
        error("truncated block");
        // (e) counts that disagree and duplicated rows.
        for (case, from, to) in [
            ("blocks count", "blocks 2", "blocks 3"),
            ("patterns count", "patterns 1", "patterns 2"),
            ("routing count", "routing 2", "routing 1"),
            (
                "block summands",
                "summands=1\tfile=boundary-1-000000.block-001",
                "summands=2\tfile=boundary-1-000000.block-001",
            ),
        ] {
            write();
            let damaged = text().replacen(from, to, 1);
            assert_ne!(damaged, text(), "{case}: fixture");
            std::fs::write(&header, damaged).unwrap();
            error(case);
        }
        write();
        let full = text();
        let topology = full
            .lines()
            .find(|l| l.starts_with("pattern_topology"))
            .unwrap()
            .to_string();
        std::fs::write(
            &header,
            full.replacen(&topology, &format!("{topology}\n{topology}"), 1),
        )
        .unwrap();
        assert!(error("duplicate topology").contains("more than once"));
        // (g) single-valued fields at most once, unknown fields refused
        // (Astra acceptance: two contradictory closure_mass_squared lines
        // were read, the last one winning).
        for (case, extra, needle) in [
            (
                "duplicate contradictory",
                "closure_mass_squared 7",
                "closure_mass_squared occurs more than once",
            ),
            (
                "duplicate identical",
                "closure_mass_squared 1",
                "closure_mass_squared occurs more than once",
            ),
            (
                "duplicate epsilon_terms",
                "epsilon_terms 3",
                "epsilon_terms occurs more than once",
            ),
            (
                "duplicate route count",
                "routing 2",
                "routing occurs more than once",
            ),
            ("unknown field", "closure_mass 1", "unknown field"),
        ] {
            write();
            let full = text();
            std::fs::write(
                &header,
                full.replacen("blocks 2", &format!("{extra}\nblocks 2"), 1),
            )
            .unwrap();
            let message = error(case);
            assert!(message.contains(needle), "{case}: {message}");
        }
        write();
        let full = text();
        let row = full
            .lines()
            .find(|l| l.starts_with("pattern 0"))
            .unwrap()
            .to_string();
        std::fs::write(
            &header,
            full.replacen(&row, &format!("{row}\tclosure_power=2"), 1),
        )
        .unwrap();
        assert!(
            error("field repeated inside a row").contains("closure_power occurs more than once")
        );
        write();
        let full = text();
        let term = full
            .lines()
            .find(|l| l.starts_with("pattern_term"))
            .unwrap()
            .to_string();
        std::fs::write(&header, full.replacen(&term, &format!("{term}\textra"), 1)).unwrap();
        assert!(error("row arity").contains("fields, 3 expected"));
        // Additive rows (Astra acceptance): a second row of the same pattern on
        // the same block (here the identical row: the terms add), terms=
        // updated: accepted.
        write();
        let full = text();
        let term = full
            .lines()
            .find(|l| l.starts_with("pattern_term 0\t1"))
            .unwrap()
            .to_string();
        let split = format!("{term}\n{term}");
        std::fs::write(
            &header,
            full.replacen(&term, &split, 1)
                .replacen("terms=2", "terms=3", 1),
        )
        .unwrap();
        let record =
            read_boundary_record(&header).expect("additive rows on one block are legitimate");
        assert_eq!(record.boundary.patterns[0].terms.len(), 3);
        // (f) SHA256SUMS: the intact record reads; a changed block is refused.
        write();
        let sums = [header.clone(), block(0), block(1)]
            .iter()
            .map(|path| {
                format!(
                    "{}  {}",
                    file_sha256(path).unwrap(),
                    path.file_name().unwrap().to_string_lossy()
                )
            })
            .collect::<Vec<_>>()
            .join("\n");
        std::fs::write(directory.join("SHA256SUMS"), sums + "\n").unwrap();
        read_boundary_record(&header).unwrap();
        let body = std::fs::read_to_string(block(1)).unwrap();
        std::fs::write(block(1), body.replace("\n3\n", "\n5\n")).unwrap();
        assert!(error("digest").contains("SHA256SUMS"));
        let _ = std::fs::remove_dir_all(&directory);
    }

    /// vk5l-physrun: read (and so check) the boundary records listed in
    /// VAKINT_TEST_BOUNDARY_CHECK (colon-separated header paths) and print
    /// what was read; skipped when unset.
    #[test]
    fn check_boundary_records_from_env() {
        let Some(list) = std::env::var_os("VAKINT_TEST_BOUNDARY_CHECK") else {
            return;
        };
        for header in std::env::split_paths(&list) {
            let started = std::time::Instant::now();
            let record = read_boundary_record(&header)
                .unwrap_or_else(|error| panic!("{}: REFUSED {error:?}", header.display()));
            let terms: usize = record.boundary.patterns.iter().map(|p| p.terms.len()).sum();
            eprintln!(
                "BOUNDARY_CHECK OK {}: blocks {}, patterns {}, pattern_term rows {terms}, routes {}, {:.1} s",
                header.display(),
                record.boundary.blocks.len(),
                record.boundary.patterns.len(),
                record.boundary.routing.len(),
                started.elapsed().as_secs_f64()
            );
        }
    }

    /// vk5l-physrun (Astra acceptance 2026-11-18): everything the writer can
    /// write, the reader accepts and reproduces. Pseudo-random boundaries
    /// (fixed seeds): 1-4 blocks, inline or in files, 1-4 patterns, terms on
    /// repeated blocks with their own coefficients (additive), identical
    /// repeated rows, compressed and plain block files.
    #[test]
    fn everything_written_reads_back() {
        crate::Vakint::initialize_vakint_symbols();
        let directory =
            std::env::temp_dir().join(format!("vakint-boundary-roundtrip-{}", std::process::id()));
        let mut state = 0x9e37_79b9_7f4a_7c15_u64;
        let mut next = |bound: u64| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) % bound
        };
        let block_pool = [
            "1",
            "3",
            "dot(k(1),k(1))",
            "dot(k(1),p(1))+2",
            "5/7*dot(p(1),p(1))",
        ];
        let coefficient_pool = ["1", "ε", "2*ε^2", "-1/3", "(4-2*ε)^(-1)"];
        let topology_pool = [
            "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,1))",
            "topo(prop(1,edge(1,2),k(1),0,2)*prop(2,edge(2,1),k(1)-p(1),0,1))",
            "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,2))",
        ];
        for case in 0..24 {
            let _ = std::fs::remove_dir_all(&directory);
            std::fs::create_dir_all(&directory).unwrap();
            let blocks: Vec<Atom> = (0..1 + next(4))
                .map(|_| vk_parse!(block_pool[next(5) as usize]).unwrap())
                .collect();
            let patterns: Vec<BoundaryPattern> = (0..1 + next(4))
                .map(|_| {
                    let mut terms: Vec<(usize, Atom)> = (0..1 + next(5))
                        .map(|_| {
                            (
                                next(blocks.len() as u64) as usize,
                                vk_parse!(coefficient_pool[next(5) as usize]).unwrap(),
                            )
                        })
                        .collect();
                    if next(3) == 0 {
                        // an identical repeated row (additive)
                        terms.push(terms[0].clone());
                    }
                    BoundaryPattern {
                        topology: vk_parse!(topology_pool[next(3) as usize]).unwrap(),
                        closure_power: 1 + next(2) as i64,
                        terms,
                    }
                })
                .collect();
            let boundary = RemainderBoundary::new(
                format!("roundtrip {case}"),
                vk_parse!("p(1)").unwrap(),
                vk_parse!("mUV^2").unwrap(),
                vec![
                    (vk_parse!("k(1)").unwrap(), vk_parse!("k(1)").unwrap()),
                    (vk_parse!("p(1)").unwrap(), vk_parse!("p(1)").unwrap()),
                ],
                blocks,
                patterns,
            )
            .unwrap();
            let header = directory.join("boundary-7-000000.txt");
            let inline_mode = next(3);
            let compress = case % 4 == 3;
            write_boundary_record(&boundary, &header, 4, compress, |block| match inline_mode {
                0 => true,
                1 => false,
                _ => block.as_view().get_byte_size() < 40,
            })
            .unwrap();
            let record = read_boundary_record(&header).unwrap_or_else(|error| {
                panic!("case {case}: the reader refused the writer's record: {error:?}")
            });
            let back = &record.boundary;
            assert_eq!(record.epsilon_terms, Some(4), "case {case}");
            assert_eq!(back.label, boundary.label, "case {case}");
            assert_eq!(back.mass_squared, boundary.mass_squared, "case {case}");
            assert_eq!(back.routing, boundary.routing, "case {case}");
            assert_eq!(back.blocks, boundary.blocks, "case {case}");
            assert_eq!(back.patterns.len(), boundary.patterns.len(), "case {case}");
            for (left, right) in back.patterns.iter().zip(&boundary.patterns) {
                assert_eq!(left.topology, right.topology, "case {case}");
                assert_eq!(left.closure_power, right.closure_power, "case {case}");
                assert_eq!(left.terms, right.terms, "case {case}");
            }
        }
        let _ = std::fs::remove_dir_all(&directory);
    }

    #[test]
    fn remainder_boundary_from_direct_round_trips() {
        crate::Vakint::initialize_vakint_symbols();
        let closure = MassiveClosure {
            momentum: vk_parse!("p(1)").unwrap(),
            mass_squared: vk_parse!("mUV^2").unwrap(),
            power: 2,
        };
        let input = DirectClosureInput {
            topology: vk_parse!("topo(prop(1,edge(1,2),k(1),0,1))").unwrap(),
            routing: vec![(vk_parse!("k(1)").unwrap(), vk_parse!("k(1)").unwrap())],
            numerator: vk_parse!("dot(k(1),p(1))^2").unwrap(),
        };
        let boundary = RemainderBoundary::from_direct(&closure, &input);
        let pattern = &boundary.patterns[0];
        let back = boundary.direct_input(pattern);
        assert_eq!(back.numerator, input.numerator);
        assert_eq!(back.topology, input.topology);
        assert_eq!(boundary.closure(pattern).power, 2);
        assert_eq!(boundary.ordinary()[0].1, input.ordinary());
    }
}
