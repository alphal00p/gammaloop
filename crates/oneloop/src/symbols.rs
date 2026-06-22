use std::sync::LazyLock;
use symbolica::{atom::Symbol, symbol};

pub struct OneLoopSymbols {
    pub ep: Symbol,
    pub d: Symbol,
    pub psq: Symbol,
    pub k: Symbol,
    pub q1: Symbol,
    pub dot: Symbol,
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
    dot: symbol!("oneloop::dot"; Symmetric, Linear),
    a0: symbol!("oneloop::A0"),
    b0: symbol!("oneloop::B0"),
    c0: symbol!("oneloop::C0"),
    d0: symbol!("oneloop::D0"),
});
