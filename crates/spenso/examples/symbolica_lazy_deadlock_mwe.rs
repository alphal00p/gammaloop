use std::sync::LazyLock;

use symbolica::{atom::Symbol, initialize, symbol};

struct MweSymbols {
    metric: Symbol,
}

static DEADLOCKING: LazyLock<MweSymbols> = LazyLock::new(|| MweSymbols {
    metric: symbol!("deadlock_mwe::metric"),
});

spenso::symbolica_init_lazy_static! {
    static SAFE, SAFE_INNER: MweSymbols = || MweSymbols {
        metric: symbol!("safe_mwe::metric"),
    };
}

initialize!(|| {
    if safe_mode() {
        spenso::symbolica_init::in_symbolica_initializer(|| {
            let _ = SAFE.force_in_initializer().metric;
        });
    } else {
        let _ = DEADLOCKING.metric;
    }
});

fn safe_mode() -> bool {
    std::env::args().any(|arg| arg == "--safe")
}

fn main() {
    if safe_mode() {
        let _ = SAFE.metric;
        println!("initializer-safe lazy access completed");
    } else {
        eprintln!("expected to deadlock: run with `-- --safe` for the fixed pattern");
        let _ = DEADLOCKING.metric;
    }
}
