use symbolica::representations::Num;

use super::Expr;

pub trait Dupe {
    fn dupe(&self) -> Self;
}

macro_rules! impl_dupe {
    ( $($t:ty),* ) => {
    $( impl Dupe for $t
    {
        fn dupe(&self) -> Self
        {
        self.clone()
        }
    }) *
    }
}

impl_dupe!(f64, f32, i64, i32, i16, i8, u64, u32, u16, u8, usize, isize);

impl<'a> Dupe for Expr<'a> {
    fn dupe(&self) -> Self {
        self.clone()
    }
}
