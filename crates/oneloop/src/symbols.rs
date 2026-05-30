use std::sync::LazyLock;
use symbolica::{atom::Symbol, symbol};

pub struct OneLoopSymbols {
    pub ep: Symbol,
    pub psq: Symbol,
    pub b0: Symbol,
}

pub static S: LazyLock<OneLoopSymbols> = LazyLock::new(|| OneLoopSymbols {
    ep: symbol!("oneloop::ep"),
    psq: symbol!("oneloop::psq"),
    b0: symbol!("oneloop::B0"),
});
