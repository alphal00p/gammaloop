//! The Dirac-trace seam of the engine: "the trace of this ordered list of
//! slashed momenta and open Lorentz indices → a Lorentz tensor polynomial".
//!
//! A closed fermion loop is a cyclic word over two letters: `Slash(e)` (the
//! slashed momentum of edge `e`) and `Index(μ)` (a gamma carrying an open
//! Lorentz index `μ`), in matrix-product order. A [`TraceEngine`] returns
//! `Tr(word)` (`Tr 1 = 4`, `D`-dimensional Clifford algebra, no γ5) as an atom
//! in source notation: metrics `g(mink(D,a), mink(D,b))`, momenta
//! `Q(e, mink(D,a))` (contracted pairs allowed) and the dimension symbol.
//!
//! * [`IdensoTrace`] (the first version, the default): idenso's production
//!   gamma simplification (`GammaSimplifySettings::canonical()`, the R*
//!   remainder chain of the beta thread) on the word alone.
//! * [`NaiveTrace`]: the textbook recursion over the first gamma's partners
//!   ((2n-1)!! terms): a cross-check of the seam, not a production engine.
//!
//! A memoised word engine (idea (e) of the plan: decompose the slashes on the
//! canonical vectors, memoise words over a small alphabet across graphs)
//! implements the same trait (a later lane).

use idenso::dirac::GammaSimplifySettings;
use idenso::tensor::{AlgebraSettings, SymbolicTensor};
use idenso::{CookMode, CookSettings};
use symbolica::atom::{Atom, AtomCore, Symbol};
use symbolica::{function, symbol};

use crate::p_integral::factorized::names;

/// One letter of a trace word.
#[derive(Clone, Debug, PartialEq)]
pub enum Letter {
    /// The slashed momentum `Q(e)` of edge `e`.
    Slash(i64),
    /// A gamma with an open Lorentz slot (`mink(D, ..)` atom).
    Index(Atom),
}

/// `Tr(word)` as a Lorentz tensor polynomial (module documentation).
pub trait TraceEngine {
    fn trace(&mut self, word: &[Letter], dimension: Symbol) -> Result<Atom, String>;
}

/// Which trace engine the lowering uses.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum TraceBackend {
    /// idenso's production gamma algebra on each trace ([`IdensoTrace`]).
    #[default]
    Idenso,
    /// The textbook recursion ([`NaiveTrace`]).
    Naive,
}

impl TraceBackend {
    pub fn engine(self) -> Box<dyn TraceEngine> {
        match self {
            Self::Idenso => Box::new(IdensoTrace),
            Self::Naive => Box::new(NaiveTrace::default()),
        }
    }
}

/// A fresh index payload `hedge(base + n, 1)` (outside every FeynGen
/// half-edge id: indices are identified by their full slot atom).
fn dummy(base: i64, n: usize) -> Atom {
    function!(symbol!(names::HEDGE), base + n as i64, 1)
}

/// The letters of `word` as a product of gamma matrices `γ(s_i, s_{i+1}, μ_i)`
/// (cyclic spinor dummies) times the slashed momenta `Q(e, μ_i)`.
fn word_atom(word: &[Letter], dimension: Symbol) -> Atom {
    let (mink, bis, gamma, q) = (
        symbol!(names::MINK),
        symbol!(names::BIS),
        symbol!(names::GAMMA),
        symbol!(names::MOMENTUM),
    );
    let d = Atom::var(dimension);
    let n = word.len();
    let spinor = |i: usize| function!(bis, 4, dummy(9_000_000, i));
    let mut product = Atom::num(1);
    for (i, letter) in word.iter().enumerate() {
        let lorentz = match letter {
            Letter::Slash(edge) => {
                let slot = function!(mink, &d, dummy(8_000_000, i));
                product *= function!(q, *edge, &slot);
                slot
            }
            Letter::Index(slot) => slot.clone(),
        };
        product *= function!(gamma, spinor(i), spinor((i + 1) % n), &lorentz);
    }
    product
}

/// idenso's production gamma algebra on one trace word.
pub struct IdensoTrace;

impl TraceEngine for IdensoTrace {
    fn trace(&mut self, word: &[Letter], dimension: Symbol) -> Result<Atom, String> {
        if word.len() % 2 == 1 {
            return Ok(Atom::Zero);
        }
        let product = word_atom(word, dimension);
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let cooked = cooking
            .try_cook(product.as_view())
            .map_err(|error| format!("idenso cooking: {error}"))?;
        let settings = AlgebraSettings {
            gamma: Some(GammaSimplifySettings::canonical()),
            epsilon: true,
            ..Default::default()
        };
        let reduced = SymbolicTensor::infer_indexed(cooked)
            .and_then(|tensor| tensor.simplify_algebra(&settings))
            .map_err(|error| format!("idenso trace: {error}"))?;
        Ok(cooking.uncook(reduced.expression().as_view()).expand())
    }
}

/// The textbook recursion `Tr(a_1 … a_n) = Σ_k (−1)^k (a_1·a_k) Tr(…)`.
#[derive(Default)]
pub struct NaiveTrace {
    dummies: usize,
}

impl NaiveTrace {
    fn pair(&mut self, a: &Letter, b: &Letter, dimension: Symbol) -> Atom {
        let (metric, mink, q) = (
            symbol!(names::METRIC),
            symbol!(names::MINK),
            symbol!(names::MOMENTUM),
        );
        match (a, b) {
            (Letter::Index(x), Letter::Index(y)) => function!(metric, x, y),
            (Letter::Index(x), Letter::Slash(e)) | (Letter::Slash(e), Letter::Index(x)) => {
                function!(q, *e, x)
            }
            (Letter::Slash(e), Letter::Slash(f)) => {
                self.dummies += 1;
                let slot = function!(mink, Atom::var(dimension), dummy(7_000_000, self.dummies));
                function!(q, *e, &slot) * function!(q, *f, &slot)
            }
        }
    }

    fn recurse(&mut self, word: &[Letter], dimension: Symbol) -> Atom {
        if word.is_empty() {
            return Atom::num(4);
        }
        if word.len() % 2 == 1 {
            return Atom::Zero;
        }
        let mut total = Atom::Zero;
        for k in 1..word.len() {
            let sign = if k % 2 == 1 { 1 } else { -1 };
            let rest = word[1..k]
                .iter()
                .chain(&word[k + 1..])
                .cloned()
                .collect::<Vec<_>>();
            let pair = self.pair(&word[0], &word[k], dimension);
            total += Atom::num(sign) * pair * self.recurse(&rest, dimension);
        }
        total
    }
}

impl TraceEngine for NaiveTrace {
    fn trace(&mut self, word: &[Letter], dimension: Symbol) -> Result<Atom, String> {
        Ok(self.recurse(word, dimension))
    }
}
