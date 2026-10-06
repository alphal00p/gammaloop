use std::cell::Cell;

use spenso_macros::track_usage;

thread_local! {
    static CALLS: Cell<usize> = const { Cell::new(0) };
}

fn record() {
    CALLS.set(CALLS.get() + 1);
}

#[track_usage(record, on_success)]
fn fallible(early: bool, input: Result<(), ()>) -> Result<u8, ()> {
    if early {
        return Ok(1);
    }
    input?;
    Ok(2)
}

struct Value(u8);

#[track_usage(record, on_success)]
impl Value {
    fn increment(&mut self) -> u8 {
        self.0 += 1;
        self.0
    }

    fn into_value(self) -> u8 {
        self.0
    }

    fn __repr__(&self) -> String {
        self.0.to_string()
    }
}

#[test]
fn records_successful_returns_and_preserves_errors_and_ownership() {
    CALLS.set(0);
    assert_eq!(fallible(false, Err(())), Err(()));
    assert_eq!(CALLS.get(), 0);
    assert_eq!(fallible(true, Err(())), Ok(1));
    assert_eq!(fallible(false, Ok(())), Ok(2));
    assert_eq!(CALLS.get(), 2);
    let mut value = Value(4);
    assert_eq!(value.__repr__(), "4");
    assert_eq!(CALLS.get(), 2);
    assert_eq!(value.increment(), 5);
    assert_eq!(value.into_value(), 5);
    assert_eq!(CALLS.get(), 4);
}

#[track_usage(record)]
fn on_entry() -> Result<(), ()> {
    Err(())
}

#[test]
fn default_mode_retains_entry_tracking() {
    CALLS.set(0);
    assert_eq!(on_entry(), Err(()));
    assert_eq!(CALLS.get(), 1);
}
