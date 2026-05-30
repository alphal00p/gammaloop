use std::sync::LazyLock;
use symbolica::{atom::Symbol, symbol};

pub struct OneLoopSymbols {
    pub ep: Symbol,
    pub psq: Symbol,
    pub a0: Symbol,
    pub b0: Symbol,
    pub c0: Symbol,
    pub d0: Symbol,
}

pub static S: LazyLock<OneLoopSymbols> = LazyLock::new(|| OneLoopSymbols {
    ep: symbol!("oneloop::ep"),
    psq: symbol!("oneloop::psq"),
    a0: symbol!("oneloop::A0"),
    b0: symbol!("oneloop::B0"),
    c0: symbol!("oneloop::C0"),
    d0: symbol!("oneloop::D0"),
});
