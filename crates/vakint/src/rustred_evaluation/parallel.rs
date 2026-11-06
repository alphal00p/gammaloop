//! Order-preserving data parallelism for independent frontend work items.

use std::sync::atomic::{AtomicUsize, Ordering};

use symbolica::utils::resources::ResourceBudget;

/// Map `items` with `task` on up to `threads` scoped workers and return the
/// results in input order, so every caller-visible effect stays independent
/// of scheduling. With `threads <= 1` or fewer than two items the map runs
/// inline on the calling thread: no pool and no spawned thread. Workers inherit
/// the caller's Symbolica resource budgets and request scope
/// ([`crate::Cancellation`]); a worker panic is re-raised here.
pub(crate) fn map_ordered<T: Sync, R: Send>(
    threads: usize,
    items: &[T],
    task: impl Fn(usize, &T) -> R + Sync,
) -> Vec<R> {
    map_ordered_by_cost(threads, items, |_| 0, task)
}

/// [`map_ordered`], with workers taking the items of largest estimated `cost`
/// first (stable for equal costs). This only shortens the makespan of uneven
/// items; results are still returned in input order. The inline path ignores
/// the estimate and runs in input order.
pub(crate) fn map_ordered_by_cost<T: Sync, R: Send>(
    threads: usize,
    items: &[T],
    cost: impl Fn(&T) -> usize,
    task: impl Fn(usize, &T) -> R + Sync,
) -> Vec<R> {
    let workers = threads.min(items.len());
    if workers <= 1 {
        return items
            .iter()
            .enumerate()
            .map(|(index, item)| task(index, item))
            .collect();
    }
    let mut order = (0..items.len()).collect::<Vec<_>>();
    order.sort_by_key(|index| std::cmp::Reverse(cost(&items[*index])));
    let budget = ResourceBudget::capture_current();
    let request = crate::Cancellation::current();
    let next = AtomicUsize::new(0);
    let mut results = (0..items.len()).map(|_| None).collect::<Vec<_>>();
    std::thread::scope(|scope| {
        let handles = (0..workers)
            .map(|_| {
                scope.spawn(|| {
                    budget.run_captured(|| {
                        let work = || {
                            let mut done = Vec::new();
                            loop {
                                let Some(&index) = order.get(next.fetch_add(1, Ordering::Relaxed))
                                else {
                                    break done;
                                };
                                done.push((index, task(index, &items[index])));
                            }
                        };
                        match &request {
                            Some(request) => request.run_coarse(work),
                            None => work(),
                        }
                    })
                })
            })
            .collect::<Vec<_>>();
        for handle in handles {
            match handle.join() {
                Ok(done) => {
                    for (index, result) in done {
                        results[index] = Some(result);
                    }
                }
                Err(panic) => std::panic::resume_unwind(panic),
            }
        }
    });
    results
        .into_iter()
        .map(|result| result.expect("every item is mapped exactly once"))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn results_keep_input_order_for_every_thread_count() {
        let items = (0..37).collect::<Vec<u64>>();
        let expected = items.iter().map(|item| item * item).collect::<Vec<_>>();
        for threads in [0, 1, 2, 5, 64] {
            assert_eq!(
                map_ordered(threads, &items, |index, item| {
                    assert_eq!(index as u64, *item);
                    item * item
                }),
                expected
            );
        }
    }

    #[test]
    fn cost_order_changes_scheduling_but_not_results() {
        let items = (0..23).collect::<Vec<usize>>();
        let started = std::sync::Mutex::new(Vec::new());
        let results = map_ordered_by_cost(
            3,
            &items,
            |item| *item % 5,
            |index, item| {
                started.lock().unwrap().push(index);
                item + 1
            },
        );
        assert_eq!(results, (1..24).collect::<Vec<_>>());
        let started = started.into_inner().unwrap();
        let mut sorted = started.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, items);
        // A worker must finish an item before claiming another, so the first
        // started item is one of the first three claims: all of cost 4.
        assert_eq!(started[0] % 5, 4);
    }

    #[test]
    fn one_thread_runs_inline_on_the_caller() {
        let caller = std::thread::current().id();
        let seen = map_ordered(1, &[(), (), ()], |_, _| std::thread::current().id());
        assert!(seen.iter().all(|id| *id == caller));
    }
}
