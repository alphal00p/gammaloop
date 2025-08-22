use std::sync::atomic::{AtomicBool, Ordering};

pub trait IsDefault {
    fn is_default(&self) -> bool;
}

pub static SHOWDEFAULTS: AtomicBool = AtomicBool::new(false);

impl<T: Default + PartialEq> IsDefault for T {
    fn is_default(&self) -> bool {
        // println!("HHHHIII {}", std::any::type_name::<T>());
        if SHOWDEFAULTS.load(Ordering::Relaxed) {
            false
        } else {
            self == &T::default()
        }
    }
}
