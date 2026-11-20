//! Per-job allocation accounting for a process that runs several jobs (one per
//! thread) under one RSS guard.
//!
//! [`Counting`] is a global allocator that adds the size of every allocation
//! made on a registered thread to that thread's live-byte count and subtracts
//! every deallocation made there. vakint does not install it (its own global
//! allocator is plain mimalloc): a binary that needs the counts installs it
//! itself and enables vakint's feature `external-global-allocator`, e.g.
//!
//! ```ignore
//! #[global_allocator]
//! static GLOBAL: vakint::alloc_account::Counting<vakint::mimalloc::MiMalloc> =
//!     vakint::alloc_account::Counting(vakint::mimalloc::MiMalloc);
//! ```
//!
//! A driver registers each job thread ([`register_current_thread_with_key`]),
//! resets its count when a job starts ([`reset_current`]) and, when the process
//! exceeds its RSS guard, compares the counts ([`live_bytes_rank`]) to fail the
//! job that holds the memory instead of whichever job reaches a checkpoint
//! first. Without [`Counting`] as the global allocator nothing is counted, and
//! [`current_live_bytes`] and [`live_bytes_rank`] say so (`None`) instead of
//! reading 0 ([`counting_installed`]).
//!
//! The counts measure requested bytes, not resident pages: a relative measure
//! between jobs.
//!
//! # Cost
//!
//! Every allocation, reallocation and deallocation through [`Counting`] makes
//! one thread-local access, to a const-initialised pointer without a destructor
//! (no lazy initialisation, no state check, no allocation). On a thread that is
//! not registered (null) that read and its branch are all. On a registered
//! thread the size change goes to the thread's counter with a relaxed load and a
//! relaxed store, not a locked add: each counter has one writer, its own
//! thread, so no update can be lost, and nothing waits for the thread's
//! outstanding stores (which made the locked add of every call expensive next
//! to code that writes the memory it allocates). Each counter sits on cache
//! lines of its own.
//!
//! # What other threads see
//!
//! Every count is exact at every moment at which its thread is not inside a
//! [`Counting`] call: readers of other threads' counters
//! ([`max_other_live_bytes`], [`live_bytes_rank`]) see the live bytes as of the
//! thread's last allocator call, as with a locked add. A registered thread's
//! counter stays in the registry after the thread exits, with its final value
//! (frees by thread-local destructors during its teardown included).
//!
//! # The RSS guard's choice
//!
//! [`live_bytes_rank`] compares the caller's count with the other threads'
//! counts, all exact, as `mine >= max(0, others)`. Equal counts are broken
//! deterministically: the thread with the smaller key
//! ([`register_current_thread_with_key`]; the registration order for
//! [`register_current_thread`]) comes first, then the earlier registration. So
//! of several threads with exactly equal counts exactly one comes first (with
//! `>=` alone all of them did).
//!
//! # Frees on another thread
//!
//! Bytes are attributed to the thread that makes the call: memory freed on
//! another thread than it was allocated on is subtracted from the freeing
//! thread's count (if that thread is registered; otherwise from none) and stays
//! counted on the allocating thread. Each thread only ever writes its own
//! counter, so a cross-thread free cannot corrupt any count (nor make an update
//! of another thread get lost); it can make the freeing thread's count
//! negative. This is rare in beta's single-threaded jobs.

use std::{
    alloc::{GlobalAlloc, Layout},
    cell::Cell,
    sync::{
        Mutex,
        atomic::{AtomicI64, AtomicU8, Ordering},
    },
};

thread_local! {
    // const-initialised and without a destructor: every access is a plain
    // thread-local access, valid for the whole life of the thread (also during
    // the teardown of other thread-locals), and it never allocates.
    static COUNTER: Cell<*const AtomicI64> = const { Cell::new(std::ptr::null()) };
}

/// A registered thread's counter, alone on its cache lines: it is written at
/// every allocator call of that thread.
#[repr(align(128))]
struct OwnLine(AtomicI64);

struct Registered {
    counter: &'static AtomicI64,
    /// Tie-break key of the RSS guard (smaller first).
    key: u64,
    /// Registration order, the second tie-break.
    sequence: u64,
}

static REGISTRY: Mutex<Vec<Registered>> = Mutex::new(Vec::new());

// `#[inline]` and `try_with`, as pr-18: with `#[inline(always)]` and `with`,
// fat LTO inlined this into the ~150k allocation sites of qcd-background-beta
// but left the thread-local accessor out of line (a call per allocation).
#[inline]
fn account(delta: i64) {
    let _ = COUNTER.try_with(|counter| {
        let pointer = counter.get();
        if !pointer.is_null() {
            // SAFETY: a registered counter is leaked ('static) and never freed;
            // the probe counter of `probe_counting` is a local that is
            // uninstalled before it goes out of scope.
            let counter = unsafe { &*pointer };
            // The calling thread is the only writer of this counter.
            counter.store(
                counter.load(Ordering::Relaxed).wrapping_add(delta),
                Ordering::Relaxed,
            );
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
        // Accounted first: the inner deallocation can then be a tail call.
        account(-(layout.size() as i64));
        unsafe { self.0.dealloc(pointer, layout) };
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
/// the thread's own, and reports whether the global allocator counted it. The
/// thread's own counter is left unchanged.
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
/// per thread; later calls keep the first registration). Its tie-break key
/// ([`live_bytes_rank`]) is its registration order.
pub fn register_current_thread() {
    register(None);
}

/// [`register_current_thread`] with the key that breaks ties between equal
/// counts in [`live_bytes_rank`] (smaller first), e.g. the driver's worker
/// index, so that the choice does not depend on which thread registered first.
pub fn register_current_thread_with_key(key: u64) {
    register(Some(key));
}

fn register(key: Option<u64>) {
    if COUNTER.with(|counter| !counter.get().is_null()) {
        return;
    }
    let counter: &'static AtomicI64 = &Box::leak(Box::new(OwnLine(AtomicI64::new(0)))).0;
    {
        let mut registry = REGISTRY.lock().unwrap_or_else(|e| e.into_inner());
        let sequence = registry.len() as u64;
        registry.push(Registered {
            counter,
            key: key.unwrap_or(sequence),
            sequence,
        });
    }
    COUNTER.with(|cell| cell.set(counter as *const AtomicI64));
}

/// Reset the calling job thread's count to 0 (a new job starts).
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

/// The calling thread's registered counter, whether or not [`Counting`] counts
/// into it.
fn own_counter() -> Option<&'static AtomicI64> {
    // SAFETY: see `account`; outside `probe_counting` the cell holds a leaked
    // registered counter or null.
    COUNTER.with(|counter| unsafe { counter.get().as_ref() })
}

/// The calling thread's count, whether or not [`Counting`] counts into it.
fn registered_live_bytes() -> Option<i64> {
    own_counter().map(|counter| counter.load(Ordering::Relaxed))
}

/// The largest live-byte count among the other registered threads (0 if none,
/// never below 0, and always 0 when [`Counting`] is not the global allocator:
/// check [`current_live_bytes`] first).
pub fn max_other_live_bytes() -> i64 {
    let mine = COUNTER.with(|counter| counter.get());
    REGISTRY
        .lock()
        .unwrap_or_else(|e| e.into_inner())
        .iter()
        .filter(|registered| !std::ptr::eq(registered.counter as *const AtomicI64, mine))
        .map(|registered| registered.counter.load(Ordering::Relaxed))
        .max()
        .unwrap_or(0)
        .max(0)
}

/// What the RSS guard compares at a checkpoint of a job thread.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LiveBytesRank {
    /// The calling thread's live bytes since its last reset.
    pub mine: i64,
    /// [`max_other_live_bytes`] at the same moment.
    pub largest_other: i64,
    /// Whether the calling thread comes first: `mine >= largest_other` (so
    /// `mine >= 0`), and every other thread with exactly `mine` bytes has a
    /// larger key (then a later registration). Of several threads with equal
    /// counts at most one comes first.
    pub largest: bool,
}

/// The calling thread's count against the other registered threads' counts;
/// `None` on a thread that is not registered, and on every thread when
/// [`Counting`] is not the global allocator.
pub fn live_bytes_rank() -> Option<LiveBytesRank> {
    if !counting_installed() {
        return None;
    }
    registered_rank()
}

fn registered_rank() -> Option<LiveBytesRank> {
    let counter = own_counter()?;
    let mine = counter.load(Ordering::Relaxed);
    let registry = REGISTRY.lock().unwrap_or_else(|e| e.into_inner());
    let me = registry
        .iter()
        .find(|registered| std::ptr::eq(registered.counter, counter))
        .map_or((u64::MAX, u64::MAX), |registered| {
            (registered.key, registered.sequence)
        });
    let others = registry
        .iter()
        .filter(|registered| !std::ptr::eq(registered.counter, counter))
        .map(|registered| {
            (
                registered.counter.load(Ordering::Relaxed),
                (registered.key, registered.sequence),
            )
        });
    let (largest_other, largest) = rank_against(mine, me, others);
    Some(LiveBytesRank {
        mine,
        largest_other,
        largest,
    })
}

/// The largest other count (at least 0) and whether `mine` comes first: iff
/// `mine >= max(0, every other count)` (the pr-18 comparison) and no other
/// thread with exactly `mine` bytes has a smaller `(key, sequence)`.
fn rank_against(
    mine: i64,
    me: (u64, u64),
    others: impl Iterator<Item = (i64, (u64, u64))>,
) -> (i64, bool) {
    let mut largest_other = 0i64;
    let mut beaten = false;
    for (bytes, order) in others {
        largest_other = largest_other.max(bytes);
        beaten |= bytes == mine && order < me;
    }
    (largest_other, mine >= largest_other && !beaten)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{Arc, Barrier};

    /// A size above any threshold an earlier, batched version of this module
    /// used (1 MiB).
    const T: i64 = 1 << 20;

    /// The calling thread's count as another thread reads it.
    fn published() -> i64 {
        own_counter().unwrap().load(Ordering::Relaxed)
    }

    /// An inner allocator without memory: it hands out dangling (aligned,
    /// non-null) pointers that are never dereferenced, so the accounting can be
    /// driven to any size.
    struct NoMemory;

    unsafe impl GlobalAlloc for NoMemory {
        unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
            std::ptr::without_provenance_mut(layout.align())
        }
        unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
            unsafe { self.alloc(layout) }
        }
        unsafe fn dealloc(&self, _pointer: *mut u8, _layout: Layout) {}
        unsafe fn realloc(&self, pointer: *mut u8, _layout: Layout, _new_size: usize) -> *mut u8 {
            pointer
        }
    }

    /// An inner allocator that always fails.
    struct Failing;

    unsafe impl GlobalAlloc for Failing {
        unsafe fn alloc(&self, _layout: Layout) -> *mut u8 {
            std::ptr::null_mut()
        }
        unsafe fn dealloc(&self, _pointer: *mut u8, _layout: Layout) {}
        unsafe fn realloc(&self, _pointer: *mut u8, _layout: Layout, _new: usize) -> *mut u8 {
            std::ptr::null_mut()
        }
    }

    fn bytes(size: usize) -> Layout {
        Layout::from_size_align(size, 1).unwrap()
    }

    /// Accounts `size` bytes (allocation) on the calling thread.
    fn grow(size: usize) {
        let pointer = unsafe { Counting(NoMemory).alloc(bytes(size)) };
        assert!(!pointer.is_null());
    }

    /// Accounts `-size` bytes (deallocation) on the calling thread.
    fn shrink(size: usize) {
        unsafe { Counting(NoMemory).dealloc(std::ptr::without_provenance_mut(1), bytes(size)) };
    }

    fn on_registered_thread<R: Send + 'static>(f: impl FnOnce() -> R + Send + 'static) -> R {
        std::thread::spawn(move || {
            register_current_thread();
            reset_current();
            f()
        })
        .join()
        .unwrap()
    }

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
        let (live, rank) = std::thread::spawn(|| {
            register_current_thread();
            reset_current();
            let _held = std::hint::black_box(vec![0u8; 1 << 20]);
            (current_live_bytes(), live_bytes_rank())
        })
        .join()
        .unwrap();
        eprintln!(
            "ALLOCATOR_PROBE binary=vakint-lib-tests global_allocator=mimalloc::MiMalloc mi_is_in_heap_region={in_mimalloc} counting_installed={installed} current_live_bytes={live:?} live_bytes_rank={rank:?}"
        );
        assert!(in_mimalloc, "the global allocator is not mimalloc");
        assert!(!installed, "vakint must not install Counting");
        assert_eq!(live, None);
        assert_eq!(rank, None);
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

    /// What another thread reads is the exact count after every call: small
    /// and large changes, upwards and downwards, below zero.
    #[test]
    fn every_change_is_visible_at_once() {
        on_registered_thread(|| {
            let t = T as usize;
            let mut live = 0i64;
            for (size, up) in [
                (t - 1, true),
                (1, true),
                (t - 1, false),
                (1, false),
                (t - 1, false),
                (1, false),
                (t + 1, true),
                (10 * t, true),
                (7, false),
                (1, true),
            ] {
                if up {
                    grow(size);
                    live += size as i64;
                } else {
                    shrink(size);
                    live -= size as i64;
                }
                assert_eq!(published(), live);
            }
            for _ in 0..1000 {
                grow(t - 1);
                shrink(t - 1);
            }
            assert_eq!(published(), live);
        });
    }

    /// The counter of a thread that has exited keeps its exact final value,
    /// including frees made by thread-local destructors during its teardown
    /// (registered before and after the thread's registration).
    #[test]
    fn exited_threads_keep_their_exact_count() {
        struct FreeOnExit(usize);
        impl Drop for FreeOnExit {
            fn drop(&mut self) {
                shrink(self.0);
            }
        }
        thread_local! {
            static BEFORE: FreeOnExit = const { FreeOnExit(30) };
            static AFTER: FreeOnExit = const { FreeOnExit(20) };
        }
        let counter = std::thread::spawn(|| {
            BEFORE.with(|_| {});
            register_current_thread();
            AFTER.with(|_| {});
            grow(5000);
            assert_eq!(published(), 5000);
            reset_current();
            assert_eq!(published(), 0);
            grow(1100);
            assert_eq!(registered_live_bytes(), Some(1100));
            assert_eq!(registered_rank().map(|rank| rank.mine), Some(1100));
            own_counter().unwrap()
        })
        .join()
        .unwrap();
        assert_eq!(counter.load(Ordering::Relaxed), 1100 - 30 - 20);
    }

    /// Readers of other threads' counters see every count exactly at the
    /// points where those threads stop, while 8 threads count concurrently; no
    /// update of any thread is lost.
    #[test]
    fn counts_are_exact_under_concurrent_threads() {
        const THREADS: usize = 8;
        const ROUNDS: usize = 200;
        let barrier = Arc::new(Barrier::new(THREADS + 1));
        let exact: Arc<Vec<AtomicI64>> =
            Arc::new((0..THREADS).map(|_| AtomicI64::new(0)).collect());
        let (sender, receiver) = std::sync::mpsc::channel();
        let workers: Vec<_> = (0..THREADS)
            .map(|worker| {
                let barrier = barrier.clone();
                let exact = exact.clone();
                let sender = sender.clone();
                std::thread::spawn(move || {
                    register_current_thread();
                    reset_current();
                    sender.send((worker, own_counter().unwrap())).unwrap();
                    let mut live = 0i64;
                    let mut state = 0x9e37_79b9_7f4a_7c15u64 ^ worker as u64;
                    for _ in 0..ROUNDS {
                        for _ in 0..2000 {
                            state = state
                                .wrapping_mul(6364136223846793005)
                                .wrapping_add(1442695040888963407);
                            let r = state >> 33;
                            let size = if r % 50 == 0 {
                                (r % (3 * T as u64)) as usize
                            } else {
                                (r % 65536) as usize
                            };
                            if (r >> 20) % 2 == 0 || live < size as i64 {
                                grow(size);
                                live += size as i64;
                            } else {
                                shrink(size);
                                live -= size as i64;
                            }
                        }
                        exact[worker].store(live, Ordering::Relaxed);
                        barrier.wait(); // round done: the checker reads
                        barrier.wait(); // checked: next round
                    }
                    assert_eq!(registered_live_bytes(), Some(live));
                })
            })
            .collect();
        drop(sender);
        let mut counters: Vec<_> = receiver.iter().take(THREADS).collect();
        counters.sort_by_key(|(worker, _)| *worker);
        for _ in 0..ROUNDS {
            barrier.wait();
            let mut largest = 0;
            for (worker, counter) in &counters {
                let exact = exact[*worker].load(Ordering::Relaxed);
                let seen = counter.load(Ordering::Relaxed);
                assert_eq!(seen, exact, "worker {worker}");
                largest = largest.max(seen);
            }
            // This thread is not registered: every registered counter is
            // another's (other tests may have registered more).
            assert!(max_other_live_bytes() >= largest);
            barrier.wait();
        }
        for worker in workers {
            worker.join().unwrap();
        }
    }

    /// A free on another thread is subtracted from the freeing thread (if it is
    /// registered) and leaves the allocating thread's count alone.
    #[test]
    fn cross_thread_frees_count_on_the_freeing_thread() {
        let small = Layout::from_size_align(5000, 8).unwrap();
        let large = Layout::from_size_align(2 * T as usize, 8).unwrap();
        // A (registered) allocates real memory: small, large, small.
        let (pointers, live_a, counter_a) = on_registered_thread(move || {
            let pointers: Vec<usize> = [small, large, small]
                .into_iter()
                .map(|layout| unsafe { Counting(std::alloc::System).alloc(layout) } as usize)
                .collect();
            assert!(pointers.iter().all(|pointer| *pointer != 0));
            (pointers, registered_live_bytes(), own_counter().unwrap())
        });
        assert_eq!(live_a, Some(2 * 5000 + 2 * T));
        // B (registered) frees the large one and one small one ...
        let (large_pointer, small_pointer) = (pointers[1], pointers[2]);
        let live_b = on_registered_thread(move || {
            let counting = Counting(std::alloc::System);
            unsafe { counting.dealloc(large_pointer as *mut u8, large) };
            assert_eq!(published(), -2 * T);
            unsafe { counting.dealloc(small_pointer as *mut u8, small) };
            registered_live_bytes()
        });
        assert_eq!(live_b, Some(-2 * T - 5000));
        // ... and C (not registered) the other small one: counted nowhere.
        let first = pointers[0];
        std::thread::spawn(move || {
            unsafe { Counting(std::alloc::System).dealloc(first as *mut u8, small) };
            assert_eq!(registered_live_bytes(), None);
        })
        .join()
        .unwrap();
        assert_eq!(counter_a.load(Ordering::Relaxed), 2 * 5000 + 2 * T);
    }

    /// realloc accounts the size change, growing and shrinking, small and
    /// large; a failed allocation or reallocation changes nothing.
    #[test]
    fn realloc_accounts_the_size_change() {
        on_registered_thread(|| {
            let counting = Counting(NoMemory);
            let pointer = unsafe { counting.alloc(bytes(100)) };
            let pointer = unsafe { counting.realloc(pointer, bytes(100), 1000) };
            assert_eq!(published(), 1000);
            let pointer = unsafe { counting.realloc(pointer, bytes(1000), T as usize + 500) };
            assert_eq!(published(), T + 500);
            let pointer = unsafe { counting.realloc(pointer, bytes(T as usize + 500), 600) };
            assert_eq!(published(), 600);
            let pointer = unsafe { counting.realloc(pointer, bytes(600), 400) };
            assert_eq!(published(), 400);
            unsafe { counting.dealloc(pointer, bytes(400)) };
            assert_eq!(registered_live_bytes(), Some(0));

            let failing = Counting(Failing);
            assert!(unsafe { failing.alloc(bytes(T as usize)) }.is_null());
            assert!(unsafe { failing.alloc_zeroed(bytes(10)) }.is_null());
            let dangling = std::ptr::without_provenance_mut(1);
            assert!(unsafe { failing.realloc(dangling, bytes(10), 3 * T as usize) }.is_null());
            assert_eq!(published(), 0);
        });
    }

    /// Zero-size changes (a realloc to the same size; a zero-size layout, which
    /// GlobalAlloc's callers never pass to a real allocator and Rust's
    /// collections never allocate) leave the count alone.
    #[test]
    fn zero_size_changes_nothing() {
        on_registered_thread(|| {
            let counting = Counting(NoMemory);
            for start in [T - 1, -(T - 1), 0] {
                reset_current();
                let pointer = unsafe { counting.alloc(bytes(64)) };
                if start >= 64 {
                    grow((start - 64) as usize);
                } else {
                    shrink((64 - start) as usize);
                }
                assert_eq!(published(), start);
                let zero = unsafe { counting.alloc(bytes(0)) };
                unsafe { counting.dealloc(zero, bytes(0)) };
                let zero = unsafe { counting.alloc_zeroed(bytes(0)) };
                unsafe { counting.dealloc(zero, bytes(0)) };
                let pointer = unsafe { counting.realloc(pointer, bytes(64), 64) };
                assert_eq!(published(), start);
                unsafe { counting.dealloc(pointer, bytes(64)) };
            }
        });
    }

    /// On a thread that is not registered nothing is counted anywhere.
    #[test]
    fn unregistered_threads_count_nothing() {
        std::thread::spawn(|| {
            for _ in 0..100 {
                grow(3 * T as usize);
                shrink(5);
            }
            assert!(own_counter().is_none());
            assert_eq!(registered_live_bytes(), None);
            assert_eq!(registered_rank(), None);
        })
        .join()
        .unwrap();
    }

    /// The probe of counting_installed leaves a registered thread's counter
    /// as it was.
    #[test]
    fn the_probe_leaves_the_thread_count_alone() {
        on_registered_thread(|| {
            grow(T as usize + 3);
            let (counter, before) = (own_counter().unwrap() as *const AtomicI64, published());
            // Counting is not this binary's global allocator.
            assert!(!probe_counting());
            assert_eq!(own_counter().unwrap() as *const AtomicI64, counter);
            assert_eq!(published(), before);
        });
    }

    /// Each registered counter lies alone on its cache lines.
    #[test]
    fn counters_have_lines_of_their_own() {
        let counter = on_registered_thread(|| own_counter().unwrap() as *const AtomicI64 as usize);
        assert_eq!(counter % 128, 0);
        assert_eq!(std::mem::size_of::<OwnLine>(), 128);
    }

    /// The guard's order: the larger count first; equal counts by key, then by
    /// registration; the pr-18 comparison `mine >= max(0, others)` otherwise.
    #[test]
    fn rank_breaks_ties_by_key_then_registration() {
        let rank =
            |mine, me, others: &[(i64, (u64, u64))]| rank_against(mine, me, others.iter().copied());
        assert_eq!(rank(10, (5, 5), &[(9, (0, 0)), (3, (9, 9))]), (9, true));
        assert_eq!(rank(10, (5, 5), &[(11, (9, 9))]), (11, false));
        assert_eq!(rank(10, (5, 5), &[(10, (6, 0))]), (10, true));
        assert_eq!(rank(10, (5, 5), &[(10, (4, 9))]), (10, false));
        assert_eq!(rank(10, (5, 5), &[(10, (5, 6))]), (10, true));
        assert_eq!(rank(10, (5, 5), &[(10, (5, 4))]), (10, false));
        // As pr-18: a count below zero never comes first, 0 does without others.
        assert_eq!(rank(-5, (0, 0), &[(-9, (1, 1))]), (0, false));
        assert_eq!(rank(0, (0, 0), &[]), (0, true));
        assert_eq!(rank(0, (3, 3), &[(0, (1, 1))]), (0, false));
        assert_eq!(rank(0, (3, 3), &[(-1, (1, 1))]), (0, true));

        // End to end, with counts far above anything else in this process.
        let huge = 1usize << 50;
        let barrier = Arc::new(Barrier::new(2));
        let run = |key: u64, extra: usize, barrier: Arc<Barrier>| {
            std::thread::spawn(move || {
                register_current_thread_with_key(key);
                reset_current();
                grow(huge + extra);
                barrier.wait(); // both counted
                let rank = registered_rank().unwrap();
                barrier.wait(); // both ranked
                shrink(huge + extra);
                rank
            })
        };
        for (extra_3, extra_7, first_3) in [(0, 0, true), (0, 1, false), (1, 0, true), (5, 5, true)]
        {
            let three = run(3, extra_3, barrier.clone());
            let seven = run(7, extra_7, barrier.clone());
            let (three, seven) = (three.join().unwrap(), seven.join().unwrap());
            assert_eq!((three.largest, seven.largest), (first_3, !first_3));
            assert_eq!(three.mine, (huge + extra_3) as i64);
            assert_eq!(three.largest_other, (huge + extra_7) as i64);
        }
    }

    /// The pr-18 hot path, for the microbenchmark: one relaxed atomic add per
    /// call on a registered thread.
    struct Pr18Counting<A>(A);

    thread_local! {
        static PR18_COUNTER: Cell<*const AtomicI64> = const { Cell::new(std::ptr::null()) };
    }

    #[inline]
    fn pr18_account(delta: i64) {
        let _ = PR18_COUNTER.try_with(|counter| {
            let pointer = counter.get();
            if !pointer.is_null() {
                unsafe { &*pointer }.fetch_add(delta, Ordering::Relaxed);
            }
        });
    }

    unsafe impl<A: GlobalAlloc> GlobalAlloc for Pr18Counting<A> {
        #[inline]
        unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
            let pointer = unsafe { self.0.alloc(layout) };
            if !pointer.is_null() {
                pr18_account(layout.size() as i64);
            }
            pointer
        }
        #[inline]
        unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
            unsafe { self.0.dealloc(pointer, layout) };
            pr18_account(-(layout.size() as i64));
        }
    }

    /// Times N allocation/free pairs through plain MiMalloc, Counting<MiMalloc>
    /// on an unregistered and on a registered thread, and the pr-18 hot path on
    /// a registered thread: one pair at a time and in batches of 64, each with
    /// the allocated bytes left alone and written after each allocation (`touch`,
    /// as code does that fills what it allocates: the stores a locked add waits
    /// for). Each variant runs on a fresh thread; best and median of several
    /// interleaved repetitions. Run on one core with
    /// `-- --ignored --nocapture --exact alloc_account::tests::microbenchmark`.
    #[cfg(feature = "mimalloc")]
    #[test]
    #[ignore = "microbenchmark"]
    fn microbenchmark() {
        use std::time::Instant;
        const PAIRS: usize = 1 << 24;
        const REPEATS: usize = 7;
        const SIZES: [usize; 8] = [16, 24, 32, 48, 64, 128, 256, 1024];

        fn run<A: GlobalAlloc>(allocator: &A, batch: usize, touch: bool) -> f64 {
            let layouts: Vec<Layout> = (0..batch)
                .map(|i| Layout::from_size_align(SIZES[i % SIZES.len()], 8).unwrap())
                .collect();
            let mut pointers = vec![std::ptr::null_mut::<u8>(); batch];
            let start = Instant::now();
            for _ in 0..PAIRS / batch {
                for (pointer, layout) in pointers.iter_mut().zip(&layouts) {
                    *pointer = std::hint::black_box(unsafe { allocator.alloc(*layout) });
                    if touch {
                        unsafe { std::ptr::write_bytes(*pointer, 0x5a, layout.size()) };
                        std::hint::black_box(*pointer);
                    }
                }
                for (pointer, layout) in pointers.iter().zip(&layouts) {
                    unsafe { allocator.dealloc(*pointer, *layout) };
                }
            }
            start.elapsed().as_secs_f64() * 1e9 / PAIRS as f64
        }

        fn on_thread(variant: &'static str, batch: usize, touch: bool) -> f64 {
            std::thread::spawn(move || match variant {
                "plain" => run(&mimalloc::MiMalloc, batch, touch),
                "unregistered" => run(&Counting(mimalloc::MiMalloc), batch, touch),
                "registered" => {
                    register_current_thread();
                    reset_current();
                    let ns = run(&Counting(mimalloc::MiMalloc), batch, touch);
                    assert_eq!(registered_live_bytes(), Some(0));
                    ns
                }
                "pr18-registered" => {
                    let counter = Box::leak(Box::new(AtomicI64::new(0)));
                    PR18_COUNTER.with(|cell| cell.set(counter));
                    let ns = run(&Pr18Counting(mimalloc::MiMalloc), batch, touch);
                    assert_eq!(counter.load(Ordering::Relaxed), 0);
                    ns
                }
                _ => unreachable!(),
            })
            .join()
            .unwrap()
        }

        let variants = ["plain", "unregistered", "registered", "pr18-registered"];
        for touch in [false, true] {
            for batch in [1, 64] {
                let mut times = vec![Vec::new(); variants.len()];
                for _ in 0..REPEATS {
                    for (variant, times) in variants.iter().zip(&mut times) {
                        times.push(on_thread(variant, batch, touch));
                    }
                }
                let plain_best = times[0].iter().copied().fold(f64::INFINITY, f64::min);
                for (variant, times) in variants.iter().zip(&mut times) {
                    times.sort_by(f64::total_cmp);
                    let (best, median) = (times[0], times[REPEATS / 2]);
                    eprintln!(
                        "ALLOC_BENCH pairs={PAIRS} touch={touch} batch={batch} variant={variant} best_ns_per_pair={best:.3} median_ns_per_pair={median:.3} best_vs_plain={:+.1}%",
                        100.0 * (best - plain_best) / plain_best
                    );
                }
            }
        }
    }
}
