//! Parsing, validation and exact normalization of p-integral inputs
//! (`docs/p_integrals.md`, section 1).

use std::collections::{BTreeMap, BTreeSet};

use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::domains::rational::Rational;

use super::{P_FAMILIES, PFamily, families, invalid};
use crate::routing::{LinearMomentum, linear_momentum};
use crate::symbols::S;
use crate::utils::vakint_macros::vk_symbol;
use crate::{VakintError, get_integer_from_atom};

/// One propagator `1/(q²)^power` of a p-integral, written `edge(from, to)`.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct PLine {
    pub id: i64,
    pub ends: (i64, i64),
    pub momentum: Atom,
    pub power: i64,
}

/// A validated p-integral graph: momentum conservation holds with the written
/// orientations, `external` enters at `entry` and leaves at `exit`.
#[derive(Clone, Debug)]
pub(crate) struct PGraph {
    pub lines: Vec<PLine>,
    /// The user's external head `p(i)`.
    pub external: Atom,
    pub entry: i64,
    pub exit: i64,
    /// Number of distinct loop momenta, equal to the cycle rank.
    pub loops: usize,
}

/// Why an integral with an external momentum is not a p-integral vakint can
/// identify (`docs/p_integrals.md`, section 1, "Out of domain").
#[derive(Debug)]
pub(crate) enum PRejection {
    /// A well-formed integral outside the p-integral domain: a massive line, a
    /// non-integer power or several external heads (and, in the catalogue,
    /// more than four loops or a graph outside the 1L-4L registry). Under
    /// `allow_unknown_integrals` it is matched as an unknown integral, as
    /// before p-integrals existed; otherwise it is an error.
    OutOfDomain(String),
    /// Malformed input (not a `prop`, a bad edge or vertex, broken momentum
    /// conservation, a loop-count mismatch, ...): always this error.
    Malformed(VakintError),
}

impl PRejection {
    pub(crate) fn is_out_of_domain(&self) -> bool {
        matches!(self, Self::OutOfDomain(_))
    }
}

impl From<VakintError> for PRejection {
    fn from(error: VakintError) -> Self {
        Self::Malformed(error)
    }
}

impl From<PRejection> for VakintError {
    fn from(rejection: PRejection) -> Self {
        match rejection {
            PRejection::OutOfDomain(detail) => invalid(detail),
            PRejection::Malformed(error) => error,
        }
    }
}

fn out_of_domain(detail: impl std::fmt::Display) -> PRejection {
    PRejection::OutOfDomain(detail.to_string())
}

/// The exact normalization of a p-integral (section 1, items 5 and 6).
#[derive(Clone, Debug)]
pub(crate) enum Normalized {
    /// The integral is scaleless and vanishes exactly.
    Zero,
    /// Tree lines, merged quadrics and pinched lines are removed; `factor`
    /// (tree factors and negative powers, in the user's momenta) multiplies the
    /// numerator. Every remaining line has a positive power.
    Graph { graph: PGraph, factor: Atom },
}

/// The `topo(...)` argument of an integral, or `None` for another expression.
fn topology_argument(integral: AtomView<'_>) -> Option<AtomView<'_>> {
    match integral {
        AtomView::Fun(function) if function.get_symbol() == S.topo && function.get_nargs() == 1 => {
            function.iter().next()
        }
        _ => None,
    }
}

fn factors(view: AtomView<'_>) -> Vec<AtomView<'_>> {
    match view {
        AtomView::Mul(product) => product.iter().collect(),
        view => vec![view],
    }
}

/// The family whose short-form head is `view`'s head, if any.
pub(super) fn short_form_family(view: AtomView) -> Option<&'static PFamily> {
    let AtomView::Fun(function) = view else {
        return None;
    };
    let symbol = function.get_symbol();
    P_FAMILIES
        .iter()
        .position(|family| symbol == vk_symbol!(family.name))
        .map(|index| &families()[index])
}

/// `prop(...)` momenta of a long form, or the head of a short form.
pub(super) fn mentions_external_momentum(integral: AtomView) -> bool {
    let Some(argument) = topology_argument(integral) else {
        return false;
    };
    if short_form_family(argument).is_some() {
        return true;
    }
    factors(argument).into_iter().any(|factor| match factor {
        AtomView::Fun(prop) if prop.get_symbol() == S.prop && prop.get_nargs() == 5 => prop
            .iter()
            .nth(2)
            .is_some_and(|momentum| momentum.contains_symbol(S.p)),
        _ => false,
    })
}

/// A vacuum long form whose every mass is literally `0`: scaleless (item 6).
pub(super) fn is_massless_vacuum(integral: AtomView) -> bool {
    let Some(argument) = topology_argument(integral) else {
        return false;
    };
    factors(argument).into_iter().all(|factor| {
        matches!(factor, AtomView::Fun(prop) if prop.get_symbol() == S.prop
            && prop.get_nargs() == 5 && prop.iter().nth(3).is_some_and(|mass| mass.is_zero()))
    })
}

fn integer(view: AtomView, what: impl Fn() -> String) -> Result<i64, VakintError> {
    get_integer_from_atom(view).ok_or_else(|| invalid(what()))
}

/// The integer coefficient of every head, rejecting anything but `k(j)`/`p(i)`.
fn line_momentum(id: i64, momentum: AtomView) -> Result<LinearMomentum, VakintError> {
    let reject = || {
        invalid(format!(
            "prop {id}: momentum {momentum} is not an integer combination of k(j) and p(i)"
        ))
    };
    let linear = linear_momentum(momentum).map_err(|_| reject())?;
    for (head, coefficient) in &linear {
        let AtomView::Fun(function) = head.as_view() else {
            return Err(reject());
        };
        if ![S.k, S.p].contains(&function.get_symbol())
            || function.get_nargs() != 1
            || get_integer_from_atom(function.iter().next().unwrap()).is_none()
            || !coefficient.is_integer()
        {
            return Err(reject());
        }
    }
    Ok(linear)
}

fn is_external(head: &Atom) -> bool {
    matches!(head.as_view(), AtomView::Fun(function) if function.get_symbol() == S.p)
}

/// `momentum` with every loop head removed.
fn loop_part(momentum: &LinearMomentum) -> LinearMomentum {
    momentum
        .iter()
        .filter(|(head, _)| !is_external(head))
        .map(|(head, coefficient)| (head.clone(), coefficient.clone()))
        .collect()
}

fn to_atom(momentum: &LinearMomentum) -> Atom {
    momentum
        .iter()
        .fold(Atom::Zero, |sum, (head, coefficient)| {
            sum + Atom::num(coefficient.clone()) * head
        })
}

/// Union-find over vertex ids; the representative of a class is its smallest id.
#[derive(Default)]
pub(super) struct VertexClasses(BTreeMap<i64, i64>);

impl VertexClasses {
    pub(super) fn find(&mut self, vertex: i64) -> i64 {
        let parent = *self.0.entry(vertex).or_insert(vertex);
        if parent == vertex {
            return vertex;
        }
        let root = self.find(parent);
        self.0.insert(vertex, root);
        root
    }

    pub(super) fn merge(&mut self, (a, b): (i64, i64)) {
        let (a, b) = (self.find(a), self.find(b));
        self.0.insert(a.max(b), a.min(b));
    }
}

/// Connected components of `vertices` under `lines`, skipping `removed`.
fn components(
    vertices: &BTreeSet<i64>,
    lines: &[(i64, i64)],
    removed: Option<i64>,
) -> Vec<BTreeSet<i64>> {
    let mut classes = VertexClasses::default();
    for &(a, b) in lines {
        if Some(a) != removed && Some(b) != removed {
            classes.merge((a, b));
        }
    }
    let mut grouped = BTreeMap::<i64, BTreeSet<i64>>::new();
    for &vertex in vertices {
        if Some(vertex) != removed {
            grouped
                .entry(classes.find(vertex))
                .or_default()
                .insert(vertex);
        }
    }
    grouped.into_values().collect()
}

/// Section 1, item 6: a normalized massless p-integral vanishes iff its two
/// external vertices coincide, a self-loop remains, a loop momentum no longer
/// flows through any line, or a vertex (possibly external) separates a part
/// that contains neither external vertex.
pub(super) fn is_scaleless(lines: &[(i64, i64)], entry: i64, exit: i64, loops: usize) -> bool {
    if entry == exit || lines.iter().any(|(a, b)| a == b) {
        return true;
    }
    let vertices = lines
        .iter()
        .flat_map(|&(a, b)| [a, b])
        .chain([entry, exit])
        .collect::<BTreeSet<_>>();
    let parts = components(&vertices, lines, None);
    if lines.len() + parts.len() < vertices.len() + loops {
        return true;
    }
    vertices
        .iter()
        .map(|&vertex| Some(vertex))
        .chain([None])
        .any(|removed| {
            components(&vertices, lines, removed)
                .iter()
                .any(|part| !part.contains(&entry) && !part.contains(&exit))
        })
}

impl PGraph {
    /// Parse and validate a long or short form p-integral (section 1, items 1-4).
    pub(crate) fn parse(integral: AtomView) -> Result<Self, VakintError> {
        Self::parse_checked(integral).map_err(VakintError::from)
    }

    /// `integral` mentions an external momentum but is outside the p-integral
    /// domain ([`PRejection::OutOfDomain`]) before any catalogue lookup.
    pub(crate) fn is_out_of_domain(integral: AtomView) -> bool {
        mentions_external_momentum(integral)
            && Self::parse_checked(integral).is_err_and(|rejection| rejection.is_out_of_domain())
    }

    /// [`Self::parse`] with the reason of a rejection: out of the p-integral
    /// domain, or malformed.
    pub(crate) fn parse_checked(integral: AtomView) -> Result<Self, PRejection> {
        let argument = topology_argument(integral)
            .ok_or_else(|| invalid(format!("expected topo(...), found {integral}")))?;
        if let Some(family) = short_form_family(argument) {
            return Self::from_lines(family.short_form_lines(argument)?);
        }
        let mut lines = Vec::new();
        for factor in factors(argument) {
            let AtomView::Fun(prop) = factor else {
                return Err(invalid(format!("{factor} is not a prop(...)")).into());
            };
            let arguments = prop.iter().collect::<Vec<_>>();
            if prop.get_symbol() != S.prop || arguments.len() != 5 {
                return Err(
                    invalid(format!("{factor} is not a prop(id,edge(a,b),q,0,power)")).into(),
                );
            }
            let id = integer(arguments[0], || {
                format!("prop ids must be integers: {factor}")
            })?;
            let ends = match arguments[1] {
                AtomView::Fun(edge) if edge.get_symbol() == S.edge && edge.get_nargs() == 2 => {
                    let mut ends = edge.iter().map(|vertex| {
                        integer(vertex, || format!("prop {id}: vertex ids must be integers"))
                    });
                    (ends.next().unwrap()?, ends.next().unwrap()?)
                }
                edge => return Err(invalid(format!("prop {id}: {edge} is not edge(a,b)")).into()),
            };
            if ends.0 < 0 || ends.1 < 0 {
                return Err(invalid(format!("prop {id}: vertex ids must be nonnegative")).into());
            }
            if !arguments[3].is_zero() {
                return Err(out_of_domain(format!(
                    "p-integrals are massless; prop {id} has mass {}",
                    arguments[3]
                )));
            }
            lines.push((id, ends, arguments[2].to_owned(), arguments[4].to_owned()));
        }
        lines.sort_by_key(|(id, ..)| *id);
        if lines.windows(2).any(|pair| pair[0].0 == pair[1].0) {
            return Err(invalid("prop ids must be distinct").into());
        }
        let mut momenta = Vec::with_capacity(lines.len());
        for (id, _, momentum, _) in &lines {
            momenta.push(line_momentum(*id, momentum.as_view())?);
        }
        let heads = momenta
            .iter()
            .flat_map(|momentum| momentum.keys().filter(|head| is_external(head)).cloned())
            .collect::<BTreeSet<_>>();
        if heads.len() != 1 {
            let detail = format!(
                "a p-integral has one external momentum; found {}",
                if heads.is_empty() {
                    "none".to_owned()
                } else {
                    heads
                        .iter()
                        .map(Atom::to_string)
                        .collect::<Vec<_>>()
                        .join(", ")
                }
            );
            // Several heads: a multi-scale integral; none: not a p-integral.
            return Err(if heads.is_empty() {
                invalid(detail).into()
            } else {
                out_of_domain(detail)
            });
        }
        // Conservation is checked before the powers, as in the spec.
        let graph = Self::validated(
            lines
                .iter()
                .map(|(id, ends, momentum, _)| (*id, *ends, momentum.clone()))
                .collect(),
            heads.into_iter().next().unwrap(),
        )?;
        let lines = lines
            .into_iter()
            .map(|(id, ends, momentum, power)| {
                let power = get_integer_from_atom(power.as_view()).ok_or_else(|| {
                    out_of_domain(format!(
                        "powers must be integers; prop {id} has power {power}"
                    ))
                })?;
                Ok(PLine {
                    id,
                    ends,
                    momentum,
                    power,
                })
            })
            .collect::<Result<Vec<_>, PRejection>>()?;
        let graph = Self { lines, ..graph };
        graph.check_loop_count()?;
        Ok(graph)
    }

    /// Lines with integer powers, as produced by a short form.
    fn from_lines(lines: Vec<PLine>) -> Result<Self, PRejection> {
        let external = lines
            .iter()
            .flat_map(|line| {
                linear_momentum(line.momentum.as_view())
                    .unwrap_or_default()
                    .into_keys()
                    .filter(is_external)
            })
            .next()
            .ok_or_else(|| invalid("a short form needs ±p(i) as its first argument"))?;
        let graph = Self::validated(
            lines
                .iter()
                .map(|line| (line.id, line.ends, line.momentum.clone()))
                .collect(),
            external,
        )?;
        let graph = Self { lines, ..graph };
        graph.check_loop_count()?;
        Ok(graph)
    }

    /// Momentum conservation with the written orientations (item 3).
    fn validated(lines: Vec<(i64, (i64, i64), Atom)>, external: Atom) -> Result<Self, VakintError> {
        let mut balances = BTreeMap::<i64, LinearMomentum>::new();
        let mut loop_heads = BTreeSet::new();
        for (id, (from, to), momentum) in &lines {
            let momentum = line_momentum(*id, momentum.as_view())?;
            for (vertex, sign) in [(*from, 1), (*to, -1)] {
                let balance = balances.entry(vertex).or_default();
                for (head, coefficient) in &momentum {
                    let entry = balance.entry(head.clone()).or_insert_with(Rational::zero);
                    *entry = &*entry + &(coefficient * &Rational::from(sign));
                }
                balance.retain(|_, coefficient| !coefficient.is_zero());
            }
            loop_heads.extend(loop_part(&momentum).into_keys());
        }
        let unit = |sign: i64| LinearMomentum::from([(external.clone(), Rational::from(sign))]);
        let (mut entry, mut exit) = (None, None);
        for (vertex, balance) in &balances {
            let slot = if balance.is_empty() {
                continue;
            } else if *balance == unit(1) {
                &mut entry
            } else if *balance == unit(-1) {
                &mut exit
            } else {
                return Err(invalid(format!(
                    "momentum conservation: {vertex} -> {}",
                    to_atom(balance)
                )));
            };
            if slot.replace(*vertex).is_some() {
                return Err(invalid(format!(
                    "momentum conservation: {vertex} -> {}",
                    to_atom(balance)
                )));
            }
        }
        let (Some(entry), Some(exit)) = (entry, exit) else {
            return Err(invalid(format!(
                "momentum conservation: {external} enters and leaves at no vertex"
            )));
        };
        Ok(Self {
            lines: lines
                .into_iter()
                .map(|(id, ends, momentum)| PLine {
                    id,
                    ends,
                    momentum,
                    power: 1,
                })
                .collect(),
            external,
            entry,
            exit,
            loops: loop_heads.len(),
        })
    }

    /// The cycle rank must equal the number of loop momenta (item 4).
    fn check_loop_count(&self) -> Result<(), VakintError> {
        let ends = self.lines.iter().map(|line| line.ends).collect::<Vec<_>>();
        let vertices = ends
            .iter()
            .flat_map(|&(a, b)| [a, b])
            .collect::<BTreeSet<_>>();
        let rank = ends.len() + components(&vertices, &ends, None).len() - vertices.len();
        if rank != self.loops {
            return Err(invalid(format!(
                "the routing does not match the graph's loop count ({} loop momenta, cycle rank {rank})",
                self.loops
            )));
        }
        Ok(())
    }

    /// Section 1, items 5 and 6: tree lines become `dot(p,p)^(-a)`, lines with
    /// identical quadrics merge into the lowest id, lines with `power <= 0` are
    /// contracted (a negative power becomes `dot(q,q)^(-a)`), and a scaleless
    /// result is an exact zero.
    pub(crate) fn normalize(self) -> Result<Normalized, VakintError> {
        let mut factor = Atom::one();
        let mut classes = VertexClasses::default();
        let mut kept: Vec<PLine> = Vec::new();
        let mut quadrics = BTreeMap::<LinearMomentum, usize>::new();
        for line in self.lines {
            let momentum = line_momentum(line.id, line.momentum.as_view())?;
            if loop_part(&momentum).is_empty() {
                if line.power != 0 {
                    if momentum.is_empty() {
                        // A massless line without momentum: scaleless.
                        return Ok(Normalized::Zero);
                    }
                    factor *= S
                        .dot(&self.external, &self.external)
                        .pow(Atom::num(-line.power));
                }
                classes.merge(line.ends);
                continue;
            }
            // Quadrics are equal iff the momenta agree up to an overall sign.
            let negative = momentum.values().next().is_some_and(|c| c.is_negative());
            let key = momentum
                .into_iter()
                .map(|(head, c)| (head, if negative { -c } else { c }))
                .collect();
            if let Some(&index) = quadrics.get(&key) {
                kept[index].power += line.power;
                classes.merge(line.ends);
            } else {
                quadrics.insert(key, kept.len());
                kept.push(line);
            }
        }
        let mut remaining = Vec::with_capacity(kept.len());
        for line in kept {
            if line.power > 0 {
                remaining.push(line);
                continue;
            }
            if line.power < 0 {
                factor *= S
                    .dot(&line.momentum, &line.momentum)
                    .pow(Atom::num(-line.power));
            }
            classes.merge(line.ends);
        }
        for line in &mut remaining {
            line.ends = (classes.find(line.ends.0), classes.find(line.ends.1));
        }
        let (entry, exit) = (classes.find(self.entry), classes.find(self.exit));
        let ends = remaining.iter().map(|line| line.ends).collect::<Vec<_>>();
        if is_scaleless(&ends, entry, exit, self.loops) {
            return Ok(Normalized::Zero);
        }
        Ok(Normalized::Graph {
            graph: Self {
                lines: remaining,
                entry,
                exit,
                ..self
            },
            factor,
        })
    }

    /// The `topo(prop(...)*...)` long form of these lines.
    pub(crate) fn to_atom(&self) -> Atom {
        S.topo
            .call(self.lines.iter().fold(Atom::one(), |product, line| {
                product
                    * S.prop.call((
                        line.id,
                        S.edge.call(line.ends),
                        &line.momentum,
                        Atom::Zero,
                        line.power,
                    ))
            }))
    }
}

impl PFamily {
    /// The lines of a short form `NAME(±p(i), a_1, ..., a_N)`.
    fn short_form_lines(&self, view: AtomView) -> Result<Vec<PLine>, PRejection> {
        let AtomView::Fun(function) = view else {
            unreachable!("short forms are functions")
        };
        let arguments = function.iter().collect::<Vec<_>>();
        let short = || PRejection::Malformed(VakintError::InvalidShortExpression(view.to_string()));
        if arguments.len() != self.lines.len() + 1 {
            return Err(short());
        }
        let external = arguments[0].to_owned();
        let momentum = linear_momentum(external.as_view()).map_err(|_| short())?;
        if momentum.len() != 1
            || !momentum.keys().all(is_external)
            || !momentum.values().all(|c| c.abs().is_one())
        {
            return Err(short());
        }
        self.lines
            .iter()
            .zip(&arguments[1..])
            .map(|(line, power)| {
                Ok(PLine {
                    momentum: line
                        .momentum
                        .replace(S.p.call(1).to_pattern())
                        .with(external.to_pattern()),
                    power: get_integer_from_atom(*power)
                        .ok_or_else(|| out_of_domain(format!("powers must be integers: {view}")))?,
                    ..line.clone()
                })
            })
            .collect()
    }
}
