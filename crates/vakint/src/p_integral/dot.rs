//! DOT export of p-integrals in FeynKit's compact convention
//! (`docs/p_integrals.md`, section 3).

use std::fmt::Write;

use symbolica::atom::{Atom, AtomView};

use super::input::{Normalized, PGraph};
use super::invalid;
use crate::VakintError;
use crate::routing::linear_momentum;
use crate::symbols::S;

/// Write a p-integral `topo(...)` as a compact FeynKit DOT graph of
/// `particle` lines. Lines with power 0 are contracted and repeated quadrics
/// merged first (section 1, item 5); an input that needs a numerator factor
/// (a tree line or a negative power) or that is scaleless is rejected.
///
/// The external momentum enters through `exte0 -> v{in}` (edge id 0) and
/// leaves through `v{out} -> exte1` (edge id 1); the lines follow in prop-id
/// order with edge ids from 2. When every loop momentum `k(j)` is carried
/// alone by a line, those lines get `lmb_id` slots in the order of `j` and are
/// oriented along `k(j)`, so that `from_diagram` reproduces the routing.
/// Numerators are not exported.
pub fn to_dot(integral: AtomView, particle: &str) -> Result<String, VakintError> {
    let normalized = match PGraph::parse(integral)?.normalize()? {
        Normalized::Zero => {
            return Err(invalid(format!(
                "a scaleless p-integral has no DOT form: {integral}"
            )));
        }
        Normalized::Graph { graph, factor } => {
            if factor != Atom::one() {
                return Err(invalid(format!(
                    "DOT export needs nonnegative powers and no tree lines: {integral}"
                )));
            }
            graph
        }
    };
    // For each loop momentum, the line carrying it alone and whether it is
    // written against that line's orientation.
    let mut loop_lines = vec![None::<(usize, bool)>; normalized.loops];
    let mut loop_heads = normalized
        .lines
        .iter()
        .flat_map(|line| {
            linear_momentum(line.momentum.as_view())
                .unwrap_or_default()
                .into_keys()
                .filter(|head| matches!(head.as_view(), AtomView::Fun(f) if f.get_symbol() == S.k))
        })
        .collect::<Vec<_>>();
    loop_heads.sort();
    loop_heads.dedup();
    for (index, line) in normalized.lines.iter().enumerate() {
        let momentum = linear_momentum(line.momentum.as_view())?;
        if momentum.len() != 1 {
            continue;
        }
        let (head, coefficient) = momentum.into_iter().next().unwrap();
        if let Some(slot) = loop_heads.iter().position(|h| *h == head)
            && coefficient.abs().is_one()
            && loop_lines[slot].is_none()
        {
            loop_lines[slot] = Some((index, coefficient.is_negative()));
        }
    }
    let basis = loop_lines.iter().all(Option::is_some);
    let mut dot = String::new();
    let mut emit = || -> std::fmt::Result {
        writeln!(dot, "digraph p_integral {{")?;
        writeln!(dot, "  exte0 [style=invis]; exte1 [style=invis];")?;
        writeln!(
            dot,
            "  exte0 -> v{} [id=0 particle=\"{particle}\"];",
            normalized.entry
        )?;
        writeln!(
            dot,
            "  v{} -> exte1 [id=1 particle=\"{particle}\"];",
            normalized.exit
        )?;
        for (index, line) in normalized.lines.iter().enumerate() {
            let slot = loop_lines
                .iter()
                .position(|entry| entry.is_some_and(|(line, _)| line == index))
                .filter(|_| basis);
            let reversed = slot.is_some_and(|slot| loop_lines[slot].unwrap().1);
            let (from, to) = if reversed {
                (line.ends.1, line.ends.0)
            } else {
                line.ends
            };
            write!(
                dot,
                "  v{from} -> v{to} [id={} particle=\"{particle}\"",
                index + 2
            )?;
            if let Some(slot) = slot {
                write!(dot, " lmb_id={slot}")?;
            }
            if line.power != 1 {
                write!(dot, " power={}", line.power)?;
            }
            writeln!(dot, "];")?;
        }
        writeln!(dot, "}}")
    };
    emit().expect("writing to a String cannot fail");
    Ok(dot)
}
