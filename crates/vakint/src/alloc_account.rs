//! Per-job allocation accounting for a process that runs several jobs (one per
//! thread) under one RSS guard.
//!
//! The global allocator ([`Counting`]) adds the size of every allocation made on
//! a registered thread to that thread's live-byte counter and subtracts every
//! deallocation made there. A driver registers each job thread
//! ([`register_current_thread`]), resets its counter when a job starts
//! ([`reset_current`]) and, when the process exceeds its RSS guard, compares the
//! counters ([`current_live_bytes`], [`max_other_live_bytes`]) to fail the job
//! that holds the memory instead of whichever job reaches a checkpoint first.
//!
//! The counters measure requested bytes, not resident pages: a relative measure
//! between jobs. Memory freed on another thread than it was allocated on is
//! subtracted from that thread (rare in the single-threaded jobs of beta).

use std::{
    alloc::{GlobalAlloc, Layout},
    cell::Cell,
    sync::{
        Mutex,
        atomic::{AtomicI64, Ordering},
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
            // SAFETY: registered counters are leaked ('static) and never freed.
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
/// that is not registered.
pub fn current_live_bytes() -> Option<i64> {
    COUNTER.with(|counter| {
        let pointer = counter.get();
        // SAFETY: see `account`.
        (!pointer.is_null()).then(|| unsafe { &*pointer }.load(Ordering::Relaxed))
    })
}

/// The largest live-byte count among the other registered threads (0 if none).
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
