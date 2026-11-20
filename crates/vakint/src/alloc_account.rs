//! Per-job allocation accounting for a process that runs several jobs (one per
//! thread) under one RSS guard.
//!
//! [`Counting`] is a global allocator that adds the size of every allocation
//! made on a registered thread to that thread's live-byte counter and subtracts
//! every deallocation made there. vakint does not install it (its own global
//! allocator is plain mimalloc): a binary that needs the counters installs it
//! itself and enables vakint's feature `external-global-allocator`, e.g.
//!
//! ```ignore
//! #[global_allocator]
//! static GLOBAL: vakint::alloc_account::Counting<vakint::mimalloc::MiMalloc> =
//!     vakint::alloc_account::Counting(vakint::mimalloc::MiMalloc);
//! ```
//!
//! A driver registers each job thread ([`register_current_thread`]), resets
//! its counter when a job starts ([`reset_current`]) and, when the process
//! exceeds its RSS guard, compares the counters ([`current_live_bytes`],
//! [`max_other_live_bytes`]) to fail the job that holds the memory instead of
//! whichever job reaches a checkpoint first. Without [`Counting`] as the global
//! allocator nothing is counted, and [`current_live_bytes`] says so (`None`)
//! instead of reading 0 ([`counting_installed`]).
//!
//! The counters measure requested bytes, not resident pages: a relative measure
//! between jobs. Memory freed on another thread than it was allocated on is
//! subtracted from that thread (rare in the single-threaded jobs of beta).

use std::{
    alloc::{GlobalAlloc, Layout},
    cell::Cell,
    sync::{
        Mutex,
        atomic::{AtomicI64, AtomicU8, Ordering},
    },
};

thread_local! {
    static COUNTER: Cell<*const AtomicI64> = const { Cell::new(std::ptr::null()) };
}

static REGISTRY: Mutex<Vec<&'static AtomicI64>> = Mutex::new(Vec::new());

#[inline]
fn account(delta: i64) {
    // try_with: no access during thread teardown; the cell has no destructor
    // and its const initializer does not allocate.
    let _ = COUNTER.try_with(|counter| {
        let pointer = counter.get();
        if !pointer.is_null() {
            // SAFETY: a registered counter is leaked ('static) and never freed;
            // the probe counter of `probe_counting` is a local that is
            // uninstalled before it goes out of scope.
            unsafe { &*pointer }.fetch_add(delta, Ordering::Relaxed);
        }
    });
}

/// A global allocator that delegates to `A` and accounts live bytes per
/// registered thread.
pub struct Counting<A>(pub A);

unsafe impl<A: GlobalAlloc> GlobalAlloc for Counting<A> {
    #[inline]
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let pointer = unsafe { self.0.alloc(layout) };
        if !pointer.is_null() {
            account(layout.size() as i64);
        }
        pointer
    }

    #[inline]
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let pointer = unsafe { self.0.alloc_zeroed(layout) };
        if !pointer.is_null() {
            account(layout.size() as i64);
        }
        pointer
    }

    #[inline]
    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        unsafe { self.0.dealloc(pointer, layout) };
        account(-(layout.size() as i64));
    }

    #[inline]
    unsafe fn realloc(&self, pointer: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        let new = unsafe { self.0.realloc(pointer, layout, new_size) };
        if !new.is_null() {
            account(new_size as i64 - layout.size() as i64);
        }
        new
    }
}

/// Whether [`Counting`] is the global allocator of this process. Probed once
/// (one heap allocation of the calling thread, counted by a temporary counter)
/// and cached; costs one relaxed load afterwards.
pub fn counting_installed() -> bool {
    const UNKNOWN: u8 = 0;
    const INSTALLED: u8 = 1;
    const ABSENT: u8 = 2;
    static STATE: AtomicU8 = AtomicU8::new(UNKNOWN);
    match STATE.load(Ordering::Relaxed) {
        INSTALLED => true,
        ABSENT => false,
        _ => {
            let installed = probe_counting();
            STATE.store(
                if installed { INSTALLED } else { ABSENT },
                Ordering::Relaxed,
            );
            installed
        }
    }
}

/// Allocates one box on the calling thread with a local counter in place of
/// the thread's own, and reports whether the global allocator counted it.
fn probe_counting() -> bool {
    const PROBE_BYTES: usize = 64;
    let probe = AtomicI64::new(0);
    COUNTER.with(|counter| {
        let previous = counter.replace(&probe as *const AtomicI64);
        // black_box: the allocation must happen (and go through the global
        // allocator) although its contents are never read.
        let boxed = std::hint::black_box(Box::new([0u8; PROBE_BYTES]));
        let counted = probe.load(Ordering::Relaxed);
        drop(boxed);
        counter.set(previous);
        counted >= PROBE_BYTES as i64
    })
}

/// Register the calling thread as a job thread with a counter of its own (once
/// per thread; later calls keep the first counter).
pub fn register_current_thread() {
    COUNTER.with(|counter| {
        if counter.get().is_null() {
            let leaked: &'static AtomicI64 = Box::leak(Box::new(AtomicI64::new(0)));
            REGISTRY
                .lock()
                .unwrap_or_else(|e| e.into_inner())
                .push(leaked);
            counter.set(leaked as *const AtomicI64);
        }
    });
}

/// Reset the calling job thread's counter (a new job starts).
pub fn reset_current() {
    COUNTER.with(|counter| {
        let pointer = counter.get();
        if !pointer.is_null() {
            // SAFETY: see `account`.
            unsafe { &*pointer }.store(0, Ordering::Relaxed);
        }
    });
}

/// The calling thread's live bytes since its last reset; `None` on a thread
/// that is not registered, and on every thread when [`Counting`] is not the
/// global allocator (then nothing is counted, and a guard must not read 0).
pub fn current_live_bytes() -> Option<i64> {
    if !counting_installed() {
        return None;
    }
    registered_live_bytes()
}

/// The calling thread's counter, whether or not [`Counting`] counts into it.
fn registered_live_bytes() -> Option<i64> {
    COUNTER.with(|counter| {
        let pointer = counter.get();
        // SAFETY: see `account`.
        (!pointer.is_null()).then(|| unsafe { &*pointer }.load(Ordering::Relaxed))
    })
}

/// The largest live-byte count among the other registered threads (0 if none,
/// and always 0 when [`Counting`] is not the global allocator: check
/// [`current_live_bytes`] first).
pub fn max_other_live_bytes() -> i64 {
    let mine = COUNTER.with(|counter| counter.get());
    REGISTRY
        .lock()
        .unwrap_or_else(|e| e.into_inner())
        .iter()
        .filter(|counter| !std::ptr::eq(**counter as *const AtomicI64, mine))
        .map(|counter| counter.load(Ordering::Relaxed))
        .max()
        .unwrap_or(0)
        .max(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// vakint's own binaries (this test binary among them) run on plain
    /// mimalloc: a heap allocation lies in mimalloc's heap, and Counting is not
    /// installed, so the counters report `None`, not 0.
    #[cfg(all(feature = "mimalloc", not(feature = "external-global-allocator")))]
    #[test]
    fn vakint_binaries_run_on_plain_mimalloc() {
        unsafe extern "C" {
            fn mi_is_in_heap_region(pointer: *const std::ffi::c_void) -> bool;
        }
        let boxed = std::hint::black_box(Box::new([7u8; 64]));
        // SAFETY: mi_is_in_heap_region only compares the address with
        // mimalloc's heap regions.
        let in_mimalloc = unsafe { mi_is_in_heap_region(boxed.as_ptr().cast()) };
        let installed = counting_installed();
        let live = std::thread::spawn(|| {
            register_current_thread();
            reset_current();
            let _held = std::hint::black_box(vec![0u8; 1 << 20]);
            current_live_bytes()
        })
        .join()
        .unwrap();
        eprintln!(
            "ALLOCATOR_PROBE binary=vakint-lib-tests global_allocator=mimalloc::MiMalloc mi_is_in_heap_region={in_mimalloc} counting_installed={installed} current_live_bytes={live:?}"
        );
        assert!(in_mimalloc, "the global allocator is not mimalloc");
        assert!(!installed, "vakint must not install Counting");
        assert_eq!(live, None);
    }

    /// The accounting itself, with Counting called directly (it is not the
    /// global allocator of this binary, so only these calls are counted).
    #[test]
    fn counting_accounts_the_registered_thread() {
        std::thread::spawn(|| {
            let counting = Counting(std::alloc::System);
            let layout = Layout::from_size_align(1000, 8).unwrap();
            // Unregistered: nothing to count into.
            let pointer = unsafe { counting.alloc(layout) };
            assert!(!pointer.is_null());
            assert_eq!(registered_live_bytes(), None);
            unsafe { counting.dealloc(pointer, layout) };

            register_current_thread();
            reset_current();
            assert_eq!(registered_live_bytes(), Some(0));
            let pointer = unsafe { counting.alloc(layout) };
            assert_eq!(registered_live_bytes(), Some(1000));
            let pointer = unsafe { counting.realloc(pointer, layout, 3000) };
            assert!(!pointer.is_null());
            assert_eq!(registered_live_bytes(), Some(3000));
            let grown = Layout::from_size_align(3000, 8).unwrap();
            let zeroed = unsafe { counting.alloc_zeroed(layout) };
            assert_eq!(registered_live_bytes(), Some(4000));
            unsafe { counting.dealloc(zeroed, layout) };
            unsafe { counting.dealloc(pointer, grown) };
            assert_eq!(registered_live_bytes(), Some(0));
            let pointer = unsafe { counting.alloc(layout) };
            reset_current();
            assert_eq!(registered_live_bytes(), Some(0));
            unsafe { counting.dealloc(pointer, layout) };
            assert_eq!(registered_live_bytes(), Some(-1000));
        })
        .join()
        .unwrap();
    }
}
