//! The versioned text form `pint6.factorized-numerator.v1` of a
//! [`FactorizedNumerator`] (N/notes/TEXTFORM-v2.md documents it).
//!
//! Conventions kept from `pint6.remainder-boundary.v1`: a `schema` first
//! line, `key<TAB>fields` lines (the reader also accepts a space after the
//! key, as v1's), `name=value` fields, atoms in Symbolica's canonical string
//! form (`to_canonical_string`), a closing `end` line. Different from v1, by
//! design: the line order is fixed (sections in a fixed order, entries sorted
//! by id or canonical string), every section is counted, a `lines` assertion
//! closes the record, an unknown key is an error (v1 ignores it), and a line
//! carries at most one atom, always as its last field.

use super::{
    ClosureSpec, ColourBranchInfo, DirectRecord, Endpoint, FACTORIZED_SCHEMA, FACTORIZED_SCHEMA_V3,
    FactorizedEdge, FactorizedError, FactorizedKind, FactorizedNumerator, FactorizedVertex,
    Insertion, LineKind, Prefactor, Projector, Result, SoftLeg, TaylorOrders, TaylorSpec,
};
use symbolica::atom::{Atom, AtomCore, AtomView};

fn one_atom(atom: &Atom, what: &str) -> Result<String> {
    let text = atom.to_canonical_string();
    if text.contains(['\n', '\r', '\t']) {
        return Err(FactorizedError::TextSyntax {
            line: 0,
            detail: format!("{what}: the canonical string contains a newline or a TAB"),
        });
    }
    Ok(text)
}

fn token(text: &str, what: &str) -> Result<()> {
    if text.is_empty() || text.contains(['\n', '\r', '\t', '=', ',']) {
        return Err(FactorizedError::TextSyntax {
            line: 0,
            detail: format!("{what} {text:?} is empty or contains TAB, newline, '=' or ','"),
        });
    }
    Ok(())
}

fn endpoint(end: Endpoint) -> String {
    match end {
        Endpoint::Vertex(v) => v.to_string(),
        Endpoint::External => "ext".into(),
    }
}

fn hedge(h: Option<usize>) -> String {
    h.map_or_else(|| "-".into(), |h| h.to_string())
}

fn optional_edge(e: Option<usize>) -> String {
    e.map_or_else(|| "none".into(), |e| e.to_string())
}

fn list<T: ToString>(values: &[T]) -> String {
    values
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

/// The ports a closure record leaves open: those its projector does not close
/// (the `open` lines of v3, derived), sorted canonical strings.
fn open_slots(value: &FactorizedNumerator) -> Result<Vec<String>> {
    let mut open = value
        .open_slots()?
        .iter()
        .map(|p| one_atom(p, "open slot"))
        .collect::<Result<Vec<_>>>()?;
    open.sort();
    Ok(open)
}

/// Write one record: v1 for a value without new features (unchanged), v3
/// for every other (kind closure, insertions, a direct record).
pub(super) fn write_record(value: &FactorizedNumerator) -> Result<String> {
    let closure = match &value.kind {
        FactorizedKind::Closure(spec) => Some(spec),
        FactorizedKind::Direct => None,
    };
    if value.label.contains(['\n', '\r', '\t']) {
        return Err(FactorizedError::TextSyntax {
            line: 0,
            detail: "the label contains a TAB or a newline".into(),
        });
    }
    let mut lines = Vec::<String>::new();
    match closure {
        None if value.has_new_features() => {
            lines.push(format!("schema\t{FACTORIZED_SCHEMA_V3}"));
            lines.push(format!("label\t{}", value.label));
            lines.push("kind\tdirect\texternals=conserved".into());
        }
        None => {
            lines.push(format!("schema\t{FACTORIZED_SCHEMA}"));
            lines.push(format!("label\t{}", value.label));
            lines.push("kind\tdirect".into());
        }
        Some(spec) => {
            lines.push(format!("schema\t{FACTORIZED_SCHEMA_V3}"));
            lines.push(format!("label\t{}", value.label));
            lines.push(if spec.independent_externals {
                "kind\tclosure\texternals=independent".into()
            } else {
                "kind\tclosure\texternals=conserved".into()
            });
        }
    }
    lines.push(format!(
        "dimension\t{}",
        one_atom(&Atom::var(value.dimension), "dimension")?
    ));
    lines.push(format!(
        "momenta\tloops={}\texternals={}",
        value.loops(),
        value.externals()
    ));
    for (i, e) in value.loop_edges.iter().enumerate() {
        lines.push(format!("loop\t{}\tedge={}", i + 1, optional_edge(*e)));
    }
    for (j, e) in value.external_edges.iter().enumerate() {
        lines.push(format!("external\t{}\tedge={}", j + 1, optional_edge(*e)));
    }
    let mut edges = value.edges.iter().collect::<Vec<_>>();
    edges.sort_by_key(|e| e.id);
    lines.push(format!("edges\t{}", edges.len()));
    for edge in edges {
        token(&edge.particle, "particle")?;
        lines.push(format!(
            "edge\t{}\tsource={}\ttarget={}\thedges={},{}\tparticle={}\tline={}\tpower={}\tmomentum={}\t{}",
            edge.id,
            endpoint(edge.source),
            endpoint(edge.target),
            hedge(edge.source_hedge),
            hedge(edge.target_hedge),
            edge.particle,
            edge.line.as_str(),
            edge.power,
            list(&edge.momentum),
            one_atom(&edge.mass_squared, "mass")?
        ));
        lines.push(format!(
            "edge_factor\t{}\t{}",
            edge.id,
            one_atom(&edge.factor, "edge factor")?
        ));
    }
    let mut vertices = value.vertices.iter().collect::<Vec<_>>();
    vertices.sort_by_key(|v| v.id);
    lines.push(format!("vertices\t{}", vertices.len()));
    for vertex in vertices {
        token(&vertex.rule, "rule")?;
        lines.push(format!(
            "vertex\t{}\trule={}\t{}",
            vertex.id,
            vertex.rule,
            one_atom(&vertex.factor, "vertex factor")?
        ));
    }
    let mut ports = value
        .ports
        .iter()
        .map(|p| one_atom(p, "port"))
        .collect::<Result<Vec<_>>>()?;
    ports.sort();
    lines.push(format!("ports\t{}", ports.len()));
    for port in ports {
        lines.push(format!("port\t{port}"));
    }
    match &value.projector {
        None => lines.push("projector\tnone".into()),
        Some(projector) => lines.push(format!(
            "projector\tp2_power={}\t{}",
            projector.p2_power,
            one_atom(&projector.tensor, "projector")?
        )),
    }
    if value.has_new_features() {
        let open = open_slots(value)?;
        lines.push(format!("open\t{}", open.len()));
        for slot in open {
            lines.push(format!("open_slot\t{slot}"));
        }
    }
    let p = &value.prefactor;
    lines.push(format!(
        "prefactor\tloop_sign={}\tnf_power={}",
        p.loop_sign, p.nf_power
    ));
    lines.push(format!(
        "prefactor_overall\t{}",
        one_atom(&p.overall, "prefactor")?
    ));
    lines.push(format!(
        "prefactor_scalars\t{}",
        one_atom(&p.scalars, "prefactor")?
    ));
    lines.push(format!(
        "prefactor_colour\t{}",
        one_atom(&p.colour, "prefactor")?
    ));
    match &value.colour_branch {
        None => lines.push("colour_branch\tnone".into()),
        Some(info) => lines.push(format!(
            "colour_branch\t{}\tof={}\trows={}",
            info.index,
            info.count,
            list(&info.rows)
        )),
    }
    if closure.is_none() && value.has_new_features() {
        // a direct value with new features is a v3 record; one with
        // insertions but no provenance writes an empty one
        let empty = DirectRecord {
            original_cycles: None,
            depth: None,
            weights: Vec::new(),
        };
        write_direct_section(value.record.as_ref().unwrap_or(&empty), &mut lines)?;
    }
    if let Some(spec) = closure {
        write_closure_section(spec, value, &mut lines)?;
    }
    let mut insertions = value.insertions.iter().collect::<Vec<_>>();
    insertions.sort_by_key(|i| i.vertex);
    lines.push(format!("insertions\t{}", insertions.len()));
    for insertion in insertions {
        lines.push(format!(
            "insertion\t{}\thedges={}\tcycles={},{}\t{}",
            insertion.vertex,
            list(&insertion.hedges),
            insertion.cycles.0,
            insertion.cycles.1,
            one_atom(&insertion.coefficient, "insertion coefficient")?
        ));
        lines.push(format!(
            "insertion_factor\t{}\t{}",
            insertion.vertex,
            one_atom(&insertion.factor, "insertion factor")?
        ));
    }
    lines.push(format!("lines\t{}", lines.len() + 1));
    lines.push("end".into());
    let mut text = lines.join("\n");
    text.push('\n');
    Ok(text)
}

/// The closure section of a v3 DIRECT record: `closure none`, `taylor none`,
/// `softs 0`, then its provenance.
fn write_direct_section(record: &DirectRecord, lines: &mut Vec<String>) -> Result<()> {
    lines.push("closure\tnone".into());
    lines.push("taylor\tnone".into());
    lines.push("softs\t0".into());
    match record.original_cycles {
        None => lines.push("original_cycles\tnone".into()),
        Some((f, g)) => lines.push(format!("original_cycles\tfermion={f}\tghost={g}")),
    }
    if !record.weights.is_empty() {
        lines.push(format!("weights\t{}", record.weights.len()));
    }
    for (name, by_vakint, atom) in &record.weights {
        token(name, "weight name")?;
        lines.push(format!(
            "weight\t{name}\tapplied_by={}\t{}",
            if *by_vakint { "vakint" } else { "beta" },
            one_atom(atom, "weight")?
        ));
    }
    lines.push(match record.depth {
        None => "depth\tnone".into(),
        Some(depth) => format!("depth\t{depth}"),
    });
    Ok(())
}

/// The closure section of a v3 DIRECT record (read): `closure none`, `taylor
/// none`, `softs 0`, `original_cycles`, optional `weights`, `depth`.
fn read_direct_section(reader: &mut Reader<'_>) -> Result<DirectRecord> {
    for key in ["closure", "taylor"] {
        let line = reader.expect(key)?;
        if line.fields.as_slice() != ["none"] {
            return Err(reader.syntax(line.number, format!("kind direct with `{key}` not none")));
        }
    }
    let (line, n) = count_section(reader, "softs")?;
    if n != 0 {
        return Err(reader.syntax(line, "kind direct with soft legs"));
    }
    let line = reader.expect("original_cycles")?;
    let original_cycles = if line.fields.as_slice() == ["none"] {
        None
    } else {
        Some((
            number(line.number, field(&line, "fermion")?, "fermion cycles")?,
            number(line.number, field(&line, "ghost")?, "ghost cycles")?,
        ))
    };
    let n = if reader.peek_key() == Some("weights") {
        count_section(reader, "weights")?.1
    } else {
        0
    };
    let mut weights = Vec::with_capacity(reader.bounded(n));
    for _ in 0..n {
        let line = reader.expect("weight")?;
        let by = match field(&line, "applied_by")? {
            "vakint" => true,
            "beta" => false,
            other => return Err(reader.syntax(line.number, format!("applied_by={other}"))),
        };
        weights.push((
            line.fields.first().copied().unwrap_or("").to_string(),
            by,
            last_atom(&line, 2)?,
        ));
    }
    let line = reader.expect("depth")?;
    let depth = match line.fields.as_slice() {
        ["none"] => None,
        [n] => Some(number(line.number, n, "depth")?),
        _ => return Err(reader.syntax(line.number, "depth<TAB><n>")),
    };
    Ok(DirectRecord {
        original_cycles,
        depth,
        weights,
    })
}

/// The closure section of v3 (`closure` .. `depth`).
fn write_closure_section(
    spec: &ClosureSpec,
    value: &FactorizedNumerator,
    lines: &mut Vec<String>,
) -> Result<()> {
    let _ = value;
    lines.push(format!(
        "closure\tedge={}\tpower={}\t{}",
        spec.edge,
        spec.power,
        one_atom(&spec.mass_squared, "closure mass")?
    ));
    match &spec.taylor {
        None => lines.push("taylor\tnone".into()),
        Some(taylor) => lines.push(format!(
            "taylor\tdegree={}\tcompensate_mass={}\tscaled_edges={}\torders={}\t{}",
            taylor.degree,
            u8::from(taylor.compensate_mass),
            list(&taylor.scaled_edges),
            taylor.orders.as_str(),
            one_atom(&spec.auxiliary_mass_squared, "auxiliary mass")?
        )),
    }
    lines.push(format!("softs\t{}", spec.softs.len()));
    let mut softs = spec.softs.iter().collect::<Vec<_>>();
    softs.sort_by_key(|s| s.external);
    for soft in softs {
        lines.push(format!(
            "soft\t{}\tedge={}\trole={}",
            soft.external,
            soft.edge,
            if soft.root { "root" } else { "soft" }
        ));
    }
    match spec.original_cycles {
        None => lines.push("original_cycles\tnone".into()),
        Some((f, g)) => lines.push(format!("original_cycles\tfermion={f}\tghost={g}")),
    }
    if !spec.weights.is_empty() {
        lines.push(format!("weights\t{}", spec.weights.len()));
    }
    for (name, by_vakint, atom) in &spec.weights {
        token(name, "weight name")?;
        lines.push(format!(
            "weight\t{name}\tapplied_by={}\t{}",
            if *by_vakint { "vakint" } else { "beta" },
            one_atom(atom, "weight")?
        ));
    }
    lines.push(match spec.depth {
        None => "depth\tnone".into(),
        Some(depth) => format!("depth\t{depth}"),
    });
    Ok(())
}

/// Keys of the beta's v2 dialect that a v3 closure record may still carry:
/// read and ignored (vakint derives the cut), one warning per record.
const IGNORED_KEYS: [&str; 6] = [
    "closure_momentum",
    "external_square",
    "routing",
    "route",
    "cut_topology",
    "hard",
];

fn skip_ignored(reader: &mut Reader<'_>, ignored: &mut Vec<String>) -> Result<()> {
    while let Some(key) = reader.peek_key() {
        if !IGNORED_KEYS.contains(&key) {
            break;
        }
        reader.raw()?;
        if !ignored.iter().any(|k| k == key) {
            ignored.push(key.to_string());
        }
    }
    Ok(())
}

/// The closure section of a v3 record (after `colour_branch`).
fn read_closure_section(
    reader: &mut Reader<'_>,
    independent_externals: bool,
) -> Result<(ClosureSpec, Vec<String>)> {
    let mut ignored = Vec::new();
    skip_ignored(reader, &mut ignored)?;
    let line = reader.expect("closure")?;
    if line.fields.as_slice() == ["none"] {
        return Err(reader.syntax(line.number, "kind closure with `closure none`"));
    }
    let edge = number(line.number, field(&line, "edge")?, "closure edge")?;
    let power = number(line.number, field(&line, "power")?, "closure power")?;
    let mass_squared = last_atom(&line, 2)?;
    skip_ignored(reader, &mut ignored)?;
    let line = reader.expect("taylor")?;
    let (taylor, auxiliary_mass_squared) = if line.fields.as_slice() == ["none"] {
        (None, mass_squared.clone())
    } else {
        let compensate = match field(&line, "compensate_mass")? {
            "0" => false,
            "1" => true,
            other => return Err(reader.syntax(line.number, format!("compensate_mass={other}"))),
        };
        let orders = field(&line, "orders")?;
        let orders = TaylorOrders::parse(orders)
            .ok_or_else(|| reader.syntax(line.number, format!("orders={orders}")))?;
        (
            Some(TaylorSpec {
                degree: number(line.number, field(&line, "degree")?, "degree")?,
                scaled_edges: integers(line.number, field(&line, "scaled_edges")?, "scaled edge")?,
                compensate_mass: compensate,
                orders,
            }),
            last_atom(&line, 4)?,
        )
    };
    skip_ignored(reader, &mut ignored)?;
    let (_, n) = count_section(reader, "softs")?;
    let mut softs = Vec::with_capacity(reader.bounded(n));
    for _ in 0..n {
        let line = reader.expect("soft")?;
        let role = field(&line, "role")?;
        softs.push(SoftLeg {
            external: number(
                line.number,
                line.fields.first().copied().unwrap_or(""),
                "soft index",
            )?,
            edge: number(line.number, field(&line, "edge")?, "soft edge")?,
            root: match role {
                "root" => true,
                "soft" => false,
                other => return Err(reader.syntax(line.number, format!("role={other}"))),
            },
        });
    }
    let line = reader.expect("original_cycles")?;
    let original_cycles = if line.fields.as_slice() == ["none"] {
        None
    } else {
        Some((
            number(line.number, field(&line, "fermion")?, "fermion cycles")?,
            number(line.number, field(&line, "ghost")?, "ghost cycles")?,
        ))
    };
    // the weights section is optional provenance (absent = none)
    let n = if reader.peek_key() == Some("weights") {
        count_section(reader, "weights")?.1
    } else {
        0
    };
    let mut weights = Vec::with_capacity(reader.bounded(n));
    for _ in 0..n {
        let line = reader.expect("weight")?;
        let by = match field(&line, "applied_by")? {
            "vakint" => true,
            "beta" => false,
            other => return Err(reader.syntax(line.number, format!("applied_by={other}"))),
        };
        weights.push((
            line.fields.first().copied().unwrap_or("").to_string(),
            by,
            last_atom(&line, 2)?,
        ));
    }
    let line = reader.expect("depth")?;
    let depth = match line.fields.as_slice() {
        ["none"] => None,
        [n] => Some(number(line.number, n, "depth")?),
        _ => return Err(reader.syntax(line.number, "depth<TAB><n>")),
    };
    skip_ignored(reader, &mut ignored)?;
    Ok((
        ClosureSpec {
            edge,
            mass_squared,
            auxiliary_mass_squared,
            power,
            independent_externals,
            taylor,
            softs,
            original_cycles,
            depth,
            weights,
        },
        ignored,
    ))
}

/// A line split into its key and TAB-separated fields.
struct Line<'a> {
    number: usize,
    key: &'a str,
    fields: Vec<&'a str>,
}

#[derive(Clone)]
struct Reader<'a> {
    lines: std::iter::Peekable<std::iter::Enumerate<std::str::Lines<'a>>>,
    first: usize,
    read: usize,
    offset: usize,
}

/// The largest preallocation a declared section count can ask for.
const PREALLOCATION_CAP: usize = 256;

impl<'a> Reader<'a> {
    /// A preallocation for a section of `n` declared entries, at most
    /// [`PREALLOCATION_CAP`] (a declared count is the caller's: it must not
    /// size an allocation; a longer section grows its vector, and a count
    /// larger than the text is refused by the section's own count check). A
    /// fixed cap, not the number of lines left: counting those on each call
    /// made reading a many-record text quadratic.
    fn bounded(&self, n: usize) -> usize {
        n.min(PREALLOCATION_CAP)
    }

    fn syntax(&self, line: usize, detail: impl Into<String>) -> FactorizedError {
        FactorizedError::TextSyntax {
            line,
            detail: detail.into(),
        }
    }

    fn raw(&mut self) -> Result<Line<'a>> {
        let Some((index, text)) = self.lines.next() else {
            return Err(self.syntax(self.offset + 1, "unexpected end of input"));
        };
        self.read += 1;
        let number = self.offset + index + 1;
        let (key, rest) = text.split_once(['\t', ' ']).unwrap_or((text, ""));
        let fields = if rest.is_empty() {
            Vec::new()
        } else {
            rest.split('\t').collect()
        };
        Ok(Line {
            number,
            key,
            fields,
        })
    }

    fn expect(&mut self, key: &str) -> Result<Line<'a>> {
        let line = self.raw()?;
        if line.key != key {
            return Err(self.syntax(
                line.number,
                format!("expected `{key}`, found `{}`", line.key),
            ));
        }
        Ok(line)
    }

    fn peek_key(&mut self) -> Option<&'a str> {
        self.lines
            .peek()
            .map(|(_, text)| text.split_once(['\t', ' ']).map_or(*text, |(k, _)| k))
    }
}

fn field<'a>(line: &Line<'a>, name: &str) -> Result<&'a str> {
    line.fields
        .iter()
        .find_map(|f| f.strip_prefix(name).and_then(|r| r.strip_prefix('=')))
        .ok_or_else(|| FactorizedError::TextSyntax {
            line: line.number,
            detail: format!("`{}` without `{name}=`", line.key),
        })
}

fn number<T: std::str::FromStr>(line: usize, text: &str, what: &str) -> Result<T> {
    text.trim()
        .parse()
        .map_err(|_| FactorizedError::TextSyntax {
            line,
            detail: format!("bad {what} {text:?}"),
        })
}

fn atom(line: usize, text: &str) -> Result<Atom> {
    crate::utils::vakint_macros::vk_parse!(text).map_err(|error| FactorizedError::TextSyntax {
        line,
        detail: format!("cannot parse {:.80}: {error}", text),
    })
}

fn last_atom(line: &Line<'_>, position: usize) -> Result<Atom> {
    if line.fields.len() != position + 1 {
        return Err(FactorizedError::TextSyntax {
            line: line.number,
            detail: format!(
                "`{}` has {} fields, expected {}",
                line.key,
                line.fields.len(),
                position + 1
            ),
        });
    }
    atom(line.number, line.fields[position])
}

fn count_section(reader: &mut Reader<'_>, key: &str) -> Result<(usize, usize)> {
    let line = reader.expect(key)?;
    let n = number(line.number, line.fields.first().copied().unwrap_or(""), key)?;
    Ok((line.number, n))
}

fn end_point(line: usize, text: &str) -> Result<Endpoint> {
    if text == "ext" {
        Ok(Endpoint::External)
    } else {
        Ok(Endpoint::Vertex(number(line, text, "endpoint")?))
    }
}

fn hedge_of(line: usize, text: &str) -> Result<Option<usize>> {
    if text == "-" {
        Ok(None)
    } else {
        Ok(Some(number(line, text, "half-edge")?))
    }
}

fn edge_ref(line: usize, text: &str) -> Result<Option<usize>> {
    if text == "none" {
        Ok(None)
    } else {
        Ok(Some(number(line, text, "edge")?))
    }
}

fn integers<T: std::str::FromStr>(line: usize, text: &str, what: &str) -> Result<Vec<T>> {
    if text.is_empty() {
        return Ok(Vec::new());
    }
    text.split(',').map(|c| number(line, c, what)).collect()
}

fn read_one(reader: &mut Reader<'_>) -> Result<FactorizedNumerator> {
    // schema v1: the base's reader, verbatim (a v1 record has no new
    // features; insertions and the closure kind are Reserved there)
    {
        let mut probe = reader.clone();
        if let Ok(schema) = probe.raw()
            && schema.key == "schema"
            && schema.fields.join("\t") == FACTORIZED_SCHEMA
        {
            return read_one_v1(reader);
        }
    }
    let schema = reader.raw()?;
    if schema.key != "schema" {
        return Err(reader.syntax(schema.number, "a record starts with `schema`"));
    }
    let found = schema.fields.join("\t");
    let v3 = found == FACTORIZED_SCHEMA_V3;
    if !v3 && found.trim_end() == FACTORIZED_SCHEMA_V3 {
        return Err(reader.syntax(
            schema.number,
            &format!(
                "the schema {found:?} has trailing whitespace: write exactly {FACTORIZED_SCHEMA_V3}"
            ),
        ));
    }
    if found != FACTORIZED_SCHEMA && !v3 {
        return Err(FactorizedError::TextSchema {
            line: schema.number,
            found,
        });
    }
    let label = reader.expect("label")?.fields.join("\t");
    let kind = reader.expect("kind")?;
    // v3: kind closure (this reader); v1: kind direct
    let direct_v3 = v3
        && matches!(
            kind.fields.as_slice(),
            ["direct"] | ["direct", "externals=conserved"]
        );
    let independent_externals = match (v3, kind.fields.as_slice()) {
        (false, ["direct"]) => false,
        (true, ["direct"]) | (true, ["direct", "externals=conserved"]) => false,
        (false, ["closure", ..]) => {
            return Err(FactorizedError::Reserved {
                feature: "the closure kind in schema v1 (it is written in schema v3)".into(),
            });
        }
        (true, ["closure"]) | (true, ["closure", "externals=conserved"]) => false,
        (true, ["closure", "externals=independent"]) => true,
        (true, [other, ..]) if *other == "direct" || *other == "vacuum" => {
            return Err(FactorizedError::Reserved {
                feature: format!(
                    "kind {other} in schema v3 with {:?} (direct is conserved; vacuum is not a v3 kind)",
                    kind.fields
                ),
            });
        }
        _ => return Err(reader.syntax(kind.number, "unknown kind")),
    };
    let dimension_line = reader.expect("dimension")?;
    let dimension = last_atom(&dimension_line, 0)?;
    let dimension = match dimension.as_view() {
        AtomView::Var(v) => v.get_symbol(),
        _ => return Err(reader.syntax(dimension_line.number, "the dimension is not a symbol")),
    };
    let momenta = reader.expect("momenta")?;
    let loops: usize = number(momenta.number, field(&momenta, "loops")?, "loops")?;
    let externals: usize = number(momenta.number, field(&momenta, "externals")?, "externals")?;
    let mut loop_edges = Vec::with_capacity(reader.bounded(loops));
    for i in 1..=loops {
        let line = reader.expect("loop")?;
        let index: usize = number(
            line.number,
            line.fields.first().copied().unwrap_or(""),
            "loop",
        )?;
        if index != i {
            return Err(reader.syntax(line.number, format!("loop {index}, expected {i}")));
        }
        loop_edges.push(edge_ref(line.number, field(&line, "edge")?)?);
    }
    let mut external_edges = Vec::with_capacity(reader.bounded(externals));
    for j in 1..=externals {
        let line = reader.expect("external")?;
        let index: usize = number(
            line.number,
            line.fields.first().copied().unwrap_or(""),
            "external",
        )?;
        if index != j {
            return Err(reader.syntax(line.number, format!("external {index}, expected {j}")));
        }
        external_edges.push(edge_ref(line.number, field(&line, "edge")?)?);
    }
    let (edges_line, n_edges) = count_section(reader, "edges")?;
    let mut edges = Vec::with_capacity(reader.bounded(n_edges));
    while reader.peek_key() == Some("edge") {
        let line = reader.expect("edge")?;
        let id: usize = number(
            line.number,
            line.fields.first().copied().unwrap_or(""),
            "edge id",
        )?;
        let hedges = field(&line, "hedges")?;
        let (hs, ht) = hedges
            .split_once(',')
            .ok_or_else(|| reader.syntax(line.number, "hedges=<source>,<target>"))?;
        let line_kind = field(&line, "line")?;
        let edge = FactorizedEdge {
            id,
            source: end_point(line.number, field(&line, "source")?)?,
            target: end_point(line.number, field(&line, "target")?)?,
            source_hedge: hedge_of(line.number, hs)?,
            target_hedge: hedge_of(line.number, ht)?,
            particle: field(&line, "particle")?.to_string(),
            line: LineKind::parse(line_kind)
                .ok_or_else(|| reader.syntax(line.number, format!("unknown line {line_kind}")))?,
            power: number(line.number, field(&line, "power")?, "power")?,
            momentum: integers(line.number, field(&line, "momentum")?, "momentum")?,
            mass_squared: last_atom(&line, 8)?,
            factor: Atom::num(1),
        };
        let factor = reader.expect("edge_factor")?;
        let fid: usize = number(
            factor.number,
            factor.fields.first().copied().unwrap_or(""),
            "edge id",
        )?;
        if fid != id {
            return Err(reader.syntax(factor.number, format!("edge_factor {fid} after edge {id}")));
        }
        let mut edge = edge;
        edge.factor = last_atom(&factor, 1)?;
        edges.push(edge);
    }
    if edges.len() != n_edges {
        return Err(FactorizedError::TextSectionCount {
            line: edges_line,
            section: "edges".into(),
            declared: n_edges,
            found: edges.len(),
        });
    }
    let (vertices_line, n_vertices) = count_section(reader, "vertices")?;
    let mut vertices = Vec::with_capacity(reader.bounded(n_vertices));
    while reader.peek_key() == Some("vertex") {
        let line = reader.expect("vertex")?;
        vertices.push(FactorizedVertex {
            id: number(
                line.number,
                line.fields.first().copied().unwrap_or(""),
                "vertex id",
            )?,
            rule: field(&line, "rule")?.to_string(),
            factor: last_atom(&line, 2)?,
        });
    }
    if vertices.len() != n_vertices {
        return Err(FactorizedError::TextSectionCount {
            line: vertices_line,
            section: "vertices".into(),
            declared: n_vertices,
            found: vertices.len(),
        });
    }
    let (ports_line, n_ports) = count_section(reader, "ports")?;
    let mut ports = Vec::with_capacity(reader.bounded(n_ports));
    while reader.peek_key() == Some("port") {
        let line = reader.expect("port")?;
        ports.push(last_atom(&line, 0)?);
    }
    if ports.len() != n_ports {
        return Err(FactorizedError::TextSectionCount {
            line: ports_line,
            section: "ports".into(),
            declared: n_ports,
            found: ports.len(),
        });
    }
    let projector_line = reader.expect("projector")?;
    let projector = match projector_line.fields.as_slice() {
        ["none"] => None,
        _ => Some(Projector {
            p2_power: number(
                projector_line.number,
                field(&projector_line, "p2_power")?,
                "p2_power",
            )?,
            tensor: last_atom(&projector_line, 1)?,
        }),
    };
    let mut open = None;
    if v3 {
        let (line, n) = count_section(reader, "open")?;
        let mut slots = Vec::with_capacity(reader.bounded(n));
        for _ in 0..n {
            let slot = reader.expect("open_slot")?;
            slots.push(slot.fields.join("\t"));
        }
        open = Some((line, slots));
    }
    let signs = reader.expect("prefactor")?;
    let loop_sign = number(signs.number, field(&signs, "loop_sign")?, "loop_sign")?;
    let nf_power = number(signs.number, field(&signs, "nf_power")?, "nf_power")?;
    let overall = last_atom(&reader.expect("prefactor_overall")?, 0)?;
    let scalars = last_atom(&reader.expect("prefactor_scalars")?, 0)?;
    let colour = last_atom(&reader.expect("prefactor_colour")?, 0)?;
    let branch = reader.expect("colour_branch")?;
    let colour_branch = match branch.fields.as_slice() {
        ["none"] => None,
        _ => Some(ColourBranchInfo {
            index: number(
                branch.number,
                branch.fields.first().copied().unwrap_or(""),
                "colour branch",
            )?,
            count: number(branch.number, field(&branch, "of")?, "colour branch count")?,
            rows: integers(branch.number, field(&branch, "rows")?, "colour rows")?,
        }),
    };
    let mut record = None;
    let closure = if direct_v3 {
        record = Some(read_direct_section(reader)?);
        None
    } else if v3 {
        let (spec, ignored) = read_closure_section(reader, independent_externals)?;
        if !ignored.is_empty() {
            eprintln!(
                "warning: closure record {label:?}: ignored keys {ignored:?} (vakint derives the cut)"
            );
        }
        Some(spec)
    } else {
        None
    };
    let (insertions_line, n_insertions) = count_section(reader, "insertions")?;
    let mut insertions = Vec::with_capacity(reader.bounded(n_insertions));
    while reader.peek_key() == Some("insertion") {
        let line = reader.expect("insertion")?;
        let vertex: usize = number(
            line.number,
            line.fields.first().copied().unwrap_or(""),
            "insertion vertex",
        )?;
        let hedges = integers(line.number, field(&line, "hedges")?, "half-edge")?;
        let cycles = field(&line, "cycles")?;
        let (fermion, ghost) = cycles
            .split_once(',')
            .ok_or_else(|| reader.syntax(line.number, "cycles=<fermion>,<ghost>"))?;
        let coefficient = last_atom(&line, 3)?;
        let factor_line = reader.expect("insertion_factor")?;
        let fv: usize = number(
            factor_line.number,
            factor_line.fields.first().copied().unwrap_or(""),
            "insertion vertex",
        )?;
        if fv != vertex {
            return Err(reader.syntax(
                factor_line.number,
                format!("insertion_factor {fv} after insertion {vertex}"),
            ));
        }
        insertions.push(Insertion {
            vertex,
            hedges,
            cycles: (
                number(line.number, fermion, "fermion cycles")?,
                number(line.number, ghost, "ghost cycles")?,
            ),
            coefficient,
            // a compact slash gamma(a, b, Q(x, mink(D))) is opened on read
            factor: super::expand_compact_slash(&last_atom(&factor_line, 1)?, vertex as i64),
        });
    }
    if insertions.len() != n_insertions {
        return Err(FactorizedError::TextSectionCount {
            line: insertions_line,
            section: "insertions".into(),
            declared: n_insertions,
            found: insertions.len(),
        });
    }
    // The order of the insertion lines is free (the writer sorts by vertex).
    insertions.sort_by_key(|i| i.vertex);
    if let Some(key) = reader.peek_key()
        && key != "lines"
    {
        let line = reader.raw()?;
        return Err(reader.syntax(line.number, format!("unknown or misplaced key `{key}`")));
    }
    let lines = reader.expect("lines")?;
    let asserted: usize = number(
        lines.number,
        lines.fields.first().copied().unwrap_or(""),
        "line count",
    )?;
    let found = reader.read - reader.first;
    if asserted != found {
        return Err(FactorizedError::TextLineCount {
            line: lines.number,
            asserted,
            found,
        });
    }
    let end = reader.raw()?;
    if end.key != "end" || !end.fields.is_empty() {
        return Err(reader.syntax(end.number, "a record ends with `end`"));
    }
    if !insertions.is_empty()
        && record.as_ref().is_some_and(|r: &DirectRecord| {
            r.original_cycles.is_none() && r.depth.is_none() && r.weights.is_empty()
        })
    {
        record = None;
    }
    let mut value = FactorizedNumerator {
        label,
        kind: match closure {
            Some(spec) => FactorizedKind::Closure(spec),
            None => FactorizedKind::Direct,
        },
        dimension,
        loop_edges,
        external_edges,
        edges,
        vertices,
        ports,
        projector,
        prefactor: Prefactor {
            overall,
            loop_sign,
            nf_power,
            scalars,
            colour,
        },
        colour_branch,
        insertions,
        record,
    };
    // rule R4: compact slashes in the projector are written out on read,
    // before the `open` lines are checked against the value
    if value.has_new_features() {
        value.expand_compact_slashes()?;
    }
    if let Some((line, slots)) = open {
        let derived = open_slots(&value)?;
        // compare as atoms' canonical strings
        let mut given = slots
            .iter()
            .map(|s| atom(line, s).map(|a| a.to_canonical_string()))
            .collect::<Result<Vec<_>>>()?;
        given.sort();
        if given != derived {
            return Err(FactorizedError::TextSyntax {
                line,
                detail: format!(
                    "open slots {given:?} differ from the ports the projector leaves open {derived:?}"
                ),
            });
        }
    }
    Ok(value)
}

/// The base's `read_one`, verbatim (the rules of a value without new
/// features: [`FactorizedNumerator::has_new_features`]).
fn read_one_v1(reader: &mut Reader<'_>) -> Result<FactorizedNumerator> {
    let schema = reader.raw()?;
    if schema.key != "schema" {
        return Err(reader.syntax(schema.number, "a record starts with `schema`"));
    }
    let found = schema.fields.join("\t");
    if found != FACTORIZED_SCHEMA {
        return Err(FactorizedError::TextSchema {
            line: schema.number,
            found,
        });
    }
    let label = reader.expect("label")?.fields.join("\t");
    let kind = reader.expect("kind")?;
    match kind.fields.as_slice() {
        ["direct"] => {}
        ["closure", ..] => {
            return Err(FactorizedError::Reserved {
                feature: "the closure kind (R* remainder with a closing line)".into(),
            });
        }
        _ => return Err(reader.syntax(kind.number, "unknown kind")),
    }
    let dimension_line = reader.expect("dimension")?;
    let dimension = last_atom(&dimension_line, 0)?;
    let dimension = match dimension.as_view() {
        AtomView::Var(v) => v.get_symbol(),
        _ => return Err(reader.syntax(dimension_line.number, "the dimension is not a symbol")),
    };
    let momenta = reader.expect("momenta")?;
    let loops: usize = number(momenta.number, field(&momenta, "loops")?, "loops")?;
    let externals: usize = number(momenta.number, field(&momenta, "externals")?, "externals")?;
    let mut loop_edges = Vec::with_capacity(reader.bounded(loops));
    for i in 1..=loops {
        let line = reader.expect("loop")?;
        let index: usize = number(
            line.number,
            line.fields.first().copied().unwrap_or(""),
            "loop",
        )?;
        if index != i {
            return Err(reader.syntax(line.number, format!("loop {index}, expected {i}")));
        }
        loop_edges.push(edge_ref(line.number, field(&line, "edge")?)?);
    }
    let mut external_edges = Vec::with_capacity(reader.bounded(externals));
    for j in 1..=externals {
        let line = reader.expect("external")?;
        let index: usize = number(
            line.number,
            line.fields.first().copied().unwrap_or(""),
            "external",
        )?;
        if index != j {
            return Err(reader.syntax(line.number, format!("external {index}, expected {j}")));
        }
        external_edges.push(edge_ref(line.number, field(&line, "edge")?)?);
    }
    let (edges_line, n_edges) = count_section(reader, "edges")?;
    let mut edges = Vec::with_capacity(reader.bounded(n_edges));
    while reader.peek_key() == Some("edge") {
        let line = reader.expect("edge")?;
        let id: usize = number(
            line.number,
            line.fields.first().copied().unwrap_or(""),
            "edge id",
        )?;
        let hedges = field(&line, "hedges")?;
        let (hs, ht) = hedges
            .split_once(',')
            .ok_or_else(|| reader.syntax(line.number, "hedges=<source>,<target>"))?;
        let line_kind = field(&line, "line")?;
        let edge = FactorizedEdge {
            id,
            source: end_point(line.number, field(&line, "source")?)?,
            target: end_point(line.number, field(&line, "target")?)?,
            source_hedge: hedge_of(line.number, hs)?,
            target_hedge: hedge_of(line.number, ht)?,
            particle: field(&line, "particle")?.to_string(),
            line: LineKind::parse(line_kind)
                .ok_or_else(|| reader.syntax(line.number, format!("unknown line {line_kind}")))?,
            power: number(line.number, field(&line, "power")?, "power")?,
            momentum: integers(line.number, field(&line, "momentum")?, "momentum")?,
            mass_squared: last_atom(&line, 8)?,
            factor: Atom::num(1),
        };
        let factor = reader.expect("edge_factor")?;
        let fid: usize = number(
            factor.number,
            factor.fields.first().copied().unwrap_or(""),
            "edge id",
        )?;
        if fid != id {
            return Err(reader.syntax(factor.number, format!("edge_factor {fid} after edge {id}")));
        }
        let mut edge = edge;
        edge.factor = last_atom(&factor, 1)?;
        edges.push(edge);
    }
    if edges.len() != n_edges {
        return Err(FactorizedError::TextSectionCount {
            line: edges_line,
            section: "edges".into(),
            declared: n_edges,
            found: edges.len(),
        });
    }
    let (vertices_line, n_vertices) = count_section(reader, "vertices")?;
    let mut vertices = Vec::with_capacity(reader.bounded(n_vertices));
    while reader.peek_key() == Some("vertex") {
        let line = reader.expect("vertex")?;
        vertices.push(FactorizedVertex {
            id: number(
                line.number,
                line.fields.first().copied().unwrap_or(""),
                "vertex id",
            )?,
            rule: field(&line, "rule")?.to_string(),
            factor: last_atom(&line, 2)?,
        });
    }
    if vertices.len() != n_vertices {
        return Err(FactorizedError::TextSectionCount {
            line: vertices_line,
            section: "vertices".into(),
            declared: n_vertices,
            found: vertices.len(),
        });
    }
    let (ports_line, n_ports) = count_section(reader, "ports")?;
    let mut ports = Vec::with_capacity(reader.bounded(n_ports));
    while reader.peek_key() == Some("port") {
        let line = reader.expect("port")?;
        ports.push(last_atom(&line, 0)?);
    }
    if ports.len() != n_ports {
        return Err(FactorizedError::TextSectionCount {
            line: ports_line,
            section: "ports".into(),
            declared: n_ports,
            found: ports.len(),
        });
    }
    let projector_line = reader.expect("projector")?;
    let projector = match projector_line.fields.as_slice() {
        ["none"] => None,
        _ => Some(Projector {
            p2_power: number(
                projector_line.number,
                field(&projector_line, "p2_power")?,
                "p2_power",
            )?,
            tensor: last_atom(&projector_line, 1)?,
        }),
    };
    let signs = reader.expect("prefactor")?;
    let loop_sign = number(signs.number, field(&signs, "loop_sign")?, "loop_sign")?;
    let nf_power = number(signs.number, field(&signs, "nf_power")?, "nf_power")?;
    let overall = last_atom(&reader.expect("prefactor_overall")?, 0)?;
    let scalars = last_atom(&reader.expect("prefactor_scalars")?, 0)?;
    let colour = last_atom(&reader.expect("prefactor_colour")?, 0)?;
    let branch = reader.expect("colour_branch")?;
    let colour_branch = match branch.fields.as_slice() {
        ["none"] => None,
        _ => Some(ColourBranchInfo {
            index: number(
                branch.number,
                branch.fields.first().copied().unwrap_or(""),
                "colour branch",
            )?,
            count: number(branch.number, field(&branch, "of")?, "colour branch count")?,
            rows: integers(branch.number, field(&branch, "rows")?, "colour rows")?,
        }),
    };
    let insertions = reader.expect("insertions")?;
    if insertions.fields.as_slice() != ["0"] {
        return Err(FactorizedError::Reserved {
            feature: "insertion vertices (child counterterms)".into(),
        });
    }
    if let Some(key) = reader.peek_key()
        && key != "lines"
    {
        let line = reader.raw()?;
        return Err(reader.syntax(line.number, format!("unknown or misplaced key `{key}`")));
    }
    let lines = reader.expect("lines")?;
    let asserted: usize = number(
        lines.number,
        lines.fields.first().copied().unwrap_or(""),
        "line count",
    )?;
    let found = reader.read - reader.first;
    if asserted != found {
        return Err(FactorizedError::TextLineCount {
            line: lines.number,
            asserted,
            found,
        });
    }
    let end = reader.raw()?;
    if end.key != "end" || !end.fields.is_empty() {
        return Err(reader.syntax(end.number, "a record ends with `end`"));
    }
    Ok(FactorizedNumerator {
        label,
        kind: FactorizedKind::Direct,
        dimension,
        loop_edges,
        external_edges,
        edges,
        vertices,
        ports,
        projector,
        prefactor: Prefactor {
            overall,
            loop_sign,
            nf_power,
            scalars,
            colour,
        },
        colour_branch,
        insertions: Vec::new(),
        record: None,
    })
}

/// Parse exactly one record (no validation).
pub(super) fn parse_record(input: &str) -> Result<FactorizedNumerator> {
    crate::Vakint::initialize_vakint_symbols();
    let mut reader = Reader {
        lines: input.lines().enumerate().peekable(),
        first: 0,
        read: 0,
        offset: 0,
    };
    let value = read_one(&mut reader)?;
    if let Some((index, _)) = reader.lines.next() {
        return Err(FactorizedError::TextSyntax {
            line: index + 1,
            detail: "text after `end`".into(),
        });
    }
    Ok(value)
}

/// Parse a concatenation of records (e.g. all colour branches of a graph in
/// one file), without validation.
pub fn parse_factorized_records(input: &str) -> Result<Vec<FactorizedNumerator>> {
    crate::Vakint::initialize_vakint_symbols();
    let mut reader = Reader {
        lines: input.lines().enumerate().peekable(),
        first: 0,
        read: 0,
        offset: 0,
    };
    let mut values = Vec::new();
    while reader.lines.peek().is_some() {
        reader.first = reader.read;
        values.push(read_one(&mut reader)?);
    }
    Ok(values)
}

/// [`parse_factorized_records`] of a file, each record validated.
pub fn read_factorized_records(path: &std::path::Path) -> Result<Vec<FactorizedNumerator>> {
    let text = std::fs::read_to_string(path).map_err(|e| FactorizedError::TextSyntax {
        line: 0,
        detail: format!("{}: {e}", path.display()),
    })?;
    let values = parse_factorized_records(&text)?;
    for value in &values {
        value.validate()?;
    }
    Ok(values)
}
