use std::sync::LazyLock;
use symbolica::{atom::Symbol, symbol};

pub struct OneLoopSymbols {
    /// Regulator, `d = 4 - 2*ep`
    pub ep: Symbol,
    /// Spacetime dimension `d`.
    pub d: Symbol,
    /// The bubble invariant `p^2`.
    pub psq: Symbol,
    /// The loop momentum.
    pub k: Symbol,
    /// External momenta
    pub q1: Symbol,
    pub q2: Symbol,
    pub q3: Symbol,
    /// Symmetric, linear dot product for numerators
    pub dot: Symbol,
    /// The master-integral heads A0/B0/C0/D0.
    pub a0: Symbol,
    pub b0: Symbol,
    pub c0: Symbol,
    pub d0: Symbol,
}

pub static S: LazyLock<OneLoopSymbols> = LazyLock::new(|| OneLoopSymbols {
    ep: symbol!("oneloop::ep"),
    d: symbol!("oneloop::d"),
    psq: symbol!("oneloop::psq"),
    k: symbol!("oneloop::k"),
    q1: symbol!("oneloop::q1"),
    q2: symbol!("oneloop::q2"),
    q3: symbol!("oneloop::q3"),
    dot: symbol!("oneloop::dot"; Symmetric, Linear),
    a0: symbol!("oneloop::A0"),
    b0: symbol!("oneloop::B0"),
    c0: symbol!("oneloop::C0"),
    d0: symbol!("oneloop::D0"),
});
