use std::{cell::Cell, ops::Deref, sync::LazyLock};

thread_local! {
    static IN_SYMBOLICA_INITIALIZER: Cell<usize> = const { Cell::new(0) };
}

/// Runs `f` while marking the current thread as executing a Symbolica
/// initializer.
///
/// Public [`SymbolicaInitLazy`] handles skip their Symbolica initialization
/// probe while this marker is active. This lets Symbolica `initialize!`
/// callbacks force symbol bundles without recursively entering Symbolica's
/// global state initializer.
pub fn in_symbolica_initializer<T>(f: impl FnOnce() -> T) -> T {
    struct Guard;

    impl Drop for Guard {
        fn drop(&mut self) {
            IN_SYMBOLICA_INITIALIZER.with(|depth| depth.set(depth.get() - 1));
        }
    }

    IN_SYMBOLICA_INITIALIZER.with(|depth| depth.set(depth.get() + 1));
    let _guard = Guard;
    f()
}

fn is_in_symbolica_initializer() -> bool {
    IN_SYMBOLICA_INITIALIZER.with(|depth| depth.get() > 0)
}

fn ensure_symbolica_initialized() {
    if !is_in_symbolica_initializer() {
        let _ = symbolica::state::State::is_builtin("__spenso_symbolica_init_probe__");
    }
}

/// A lazy static handle whose public access initializes Symbolica before
/// forcing the wrapped lazy value.
///
/// Use this for symbol bundles that are also warmed from Symbolica
/// `initialize!` callbacks. The public handle avoids acquiring the bundle's
/// lazy lock before Symbolica callbacks have had a chance to run.
pub struct SymbolicaInitLazy<T: 'static> {
    inner: &'static LazyLock<T>,
}

impl<T: 'static> Clone for SymbolicaInitLazy<T> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<T: 'static> Copy for SymbolicaInitLazy<T> {}

impl<T: 'static> SymbolicaInitLazy<T> {
    pub const fn new(inner: &'static LazyLock<T>) -> Self {
        Self { inner }
    }

    /// Forces the wrapped lazy value from inside a Symbolica initializer.
    ///
    /// This intentionally does not probe Symbolica first; callers should use it
    /// only from registered `initialize!` callbacks or from
    /// [`in_symbolica_initializer`].
    pub fn force_in_initializer(&self) -> &T {
        in_symbolica_initializer(|| LazyLock::force(self.inner))
    }
}

impl<T: 'static> Deref for SymbolicaInitLazy<T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        ensure_symbolica_initialized();
        LazyLock::force(self.inner)
    }
}

#[macro_export]
macro_rules! symbolica_init_lazy_static {
    (
        $(#[$meta:meta])*
        $vis:vis static $public_name:ident, $inner_name:ident: $ty:ty = $init:expr;
    ) => {
        static $inner_name: ::std::sync::LazyLock<$ty> = ::std::sync::LazyLock::new($init);

        $(#[$meta])*
        $vis static $public_name: $crate::symbolica_init::SymbolicaInitLazy<$ty> =
            $crate::symbolica_init::SymbolicaInitLazy::new(&$inner_name);
    };
}
