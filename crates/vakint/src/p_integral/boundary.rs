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

/// Stream one block file (`pint6.boundary-block.v1`, zstd unless the name
/// does not end in `.zst`): a summand per line, added in chunks of 65,536;
/// the summand count of the file header must match the lines read.
pub fn read_boundary_block(path: &std::path::Path) -> Result<Atom, VakintError> {
    use std::io::BufRead;
    let io = |error: std::io::Error| super::invalid(format!("{}: {error}", path.display()));
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
    header("block")?;
    let summands: usize = header("summands")?
        .strip_prefix("summands ")
        .and_then(|count| count.trim().parse().ok())
        .ok_or_else(|| super::invalid(format!("{}: bad summands line", path.display())))?;
    let mut partial = Vec::new();
    let mut chunk = Vec::with_capacity(BLOCK_CHUNK.min(summands));
    let mut read = 0usize;
    for line in lines {
        let line = line.map_err(io)?;
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
    Ok(Atom::add_many(partial.iter()))
}

/// Read a `pint6.remainder-boundary.v1` header and its block files (same
/// directory), streaming every block ([`read_boundary_block`]). Only the keys
/// vakint needs are interpreted: closure_momentum, closure_mass_squared,
/// epsilon_terms, p_loops, label, route, block / block_atom, pattern
/// (closure_power, ref_closure), pattern_topology, pattern_term; the others
/// are informational.
pub fn read_boundary_record(path: &std::path::Path) -> Result<BoundaryRecord, VakintError> {
    // Vakint's symbols (with their attributes) before any parse.
    crate::Vakint::initialize_vakint_symbols();
    let text = std::fs::read_to_string(path)
        .map_err(|error| super::invalid(format!("{}: {error}", path.display())))?;
    let directory = path.parent().unwrap_or(std::path::Path::new("."));
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
    let mut block_files = std::collections::BTreeMap::<usize, String>::new();
    let mut block_atoms = std::collections::BTreeMap::<usize, Atom>::new();
    let mut powers = std::collections::BTreeMap::<usize, (i64, Option<String>)>::new();
    let mut topologies = std::collections::BTreeMap::<usize, Atom>::new();
    let mut terms = std::collections::BTreeMap::<usize, Vec<(usize, Atom)>>::new();
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
    for line in lines {
        let (key, rest) = line.split_once([' ', '\t']).unwrap_or((line.trim(), ""));
        let fields = rest.split('\t').collect::<Vec<_>>();
        match key {
            "label" => label = rest.to_string(),
            "closure_momentum" => momentum = Some(parse_atom(rest)?),
            "closure_mass_squared" => mass = Some(parse_atom(rest)?),
            "epsilon_terms" => epsilon_terms = rest.trim().parse().ok(),
            "p_loops" => p_loops = rest.trim().parse().ok(),
            "route" if fields.len() == 2 => {
                routing.push((parse_atom(fields[0])?, parse_atom(fields[1])?))
            }
            "block" => {
                let g = index(fields[0])?;
                let file = field(&fields, "file").unwrap_or_else(|| "inline".into());
                block_files.insert(g, file);
            }
            "block_atom" if fields.len() == 2 => {
                block_atoms.insert(index(fields[0])?, parse_atom(fields[1])?);
            }
            "pattern" => {
                let pi = index(fields[0])?;
                let power = field(&fields, "closure_power")
                    .and_then(|power| power.parse().ok())
                    .ok_or_else(|| {
                        super::invalid(format!(
                            "{}: pattern {pi} without closure_power",
                            path.display()
                        ))
                    })?;
                let reference = field(&fields, "ref_closure").filter(|seq| seq != "none");
                powers.insert(pi, (power, reference));
            }
            "pattern_topology" if fields.len() == 2 => {
                topologies.insert(index(fields[0])?, parse_atom(fields[1])?);
            }
            "pattern_term" if fields.len() == 3 => {
                terms
                    .entry(index(fields[0])?)
                    .or_default()
                    .push((index(fields[1])?, parse_atom(fields[2])?));
            }
            _ => {}
        }
    }
    let mut blocks = Vec::with_capacity(block_files.len());
    for (position, (g, file)) in block_files.iter().enumerate() {
        if *g != position {
            return Err(super::invalid(format!(
                "{}: blocks are not numbered 0..",
                path.display()
            )));
        }
        let block = if file == "inline" {
            block_atoms.remove(g).ok_or_else(|| {
                super::invalid(format!(
                    "{}: inline block {g} without block_atom",
                    path.display()
                ))
            })?
        } else {
            read_boundary_block(&directory.join(file))?
        };
        blocks.push(block);
    }
    let mut patterns = Vec::with_capacity(powers.len());
    let mut ref_closures = Vec::with_capacity(powers.len());
    for (pi, (closure_power, reference)) in powers {
        let topology = topologies.remove(&pi).ok_or_else(|| {
            super::invalid(format!(
                "{}: pattern {pi} without pattern_topology",
                path.display()
            ))
        })?;
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
