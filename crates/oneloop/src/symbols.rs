use std::sync::LazyLock;
use symbolica::{atom::Symbol, symbol};

pub struct OneLoopSymbols {
    pub ep: Symbol,
    pub psq: Symbol,
    pub musq: Symbol,
    pub log: Symbol,
}

pub static S: LazyLock<OneLoopSymbols> = LazyLock::new(|| OneLoopSymbols {
    ep: symbol!("oneloop::ep"),
    psq: symbol!("oneloop::psq"),
    musq: symbol!("oneloop::musq"),
    log: symbol!("oneloop::log"),
});
