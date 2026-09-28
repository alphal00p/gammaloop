pub use linnet::num_traits::{RefOne, RefZero};

pub trait One {
    fn one() -> Self;
}

pub trait Zero {
    fn zero() -> Self;
}

pub trait IsZero {
    fn is_zero(&self) -> bool;

    fn is_non_zero(&self) -> bool {
        !self.is_zero()
    }
}

impl<T: RefZero + PartialEq> IsZero for T {
    fn is_zero(&self) -> bool {
        self.ref_zero() == *self
    }
}
