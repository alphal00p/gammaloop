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
//! Version 1 (this file) reads each pattern on its own. Block reuse across
//! patterns (the loop-polynomial reading and slot substitution depend only on
//! the presentation, the routing and the block, the line powers only enter the
//! final keys) and chunked reading of a block with budget checkpoints are the
//! next steps; they do not change this type.
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
    /// declined (an error instead when strict).
    pub fn lower(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        options: BoundaryOptions,
    ) -> Result<(Vec<LoweredClosure>, Vec<usize>), VakintError> {
        let mut classes = Vec::<LoweredClosure>::new();
        let mut declined = Vec::new();
        for (index, pattern) in self.patterns.iter().enumerate() {
            // A checkpoint between patterns: the resource budget and the
            // execution limits of the settings apply to every pattern.
            settings.validate_execution()?;
            crate::execution::checkpoint().map_err(VakintError::from)?;
            let closure = self.closure(pattern);
            let input = self.direct_input(pattern);
            match closure.lower_direct(vakint, settings, &input)? {
                Some(lowered) => {
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
                None if options.strict() => {
                    return Err(super::invalid(format!(
                        "remainder boundary {:?}: the direct reader declined pattern {index} (strict mode)",
                        self.label
                    )));
                }
                None => declined.push(index),
            }
        }
        Ok((classes, declined))
    }

    /// The value of the remainder: Σ over patterns of the closed p-integral,
    /// one reduction per closing line.
    pub fn evaluate(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        options: BoundaryOptions,
    ) -> Result<Atom, VakintError> {
        settings.with_resource_budget(|| {
            let (classes, declined) = self.lower(vakint, settings, options)?;
            let mut result = Atom::Zero;
            for class in classes {
                result += class.evaluate(vakint, settings)?;
            }
            for index in declined {
                let pattern = &self.patterns[index];
                result += MassiveClosure::evaluate_batch(
                    vakint,
                    settings,
                    &[(self.closure(pattern), self.direct_input(pattern).ordinary())],
                )?
                .pop()
                .unwrap_or_default();
            }
            Ok(result.expand())
        })
    }
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
    fn remainder_boundary_from_direct_round_trips() {
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
