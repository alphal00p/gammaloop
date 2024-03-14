use std::{borrow::Cow, ops::Add};

trait LCM<T> {
    type Output: Clone;
    fn lcm(&self) -> Cow<Self::Output>;
}

#[derive(Clone)]
struct A {}
#[derive(Clone)]
struct B {}

impl<'a> LCM<&B> for &'a A {
    type Output = A;
    fn lcm(&self) -> Cow<Self::Output> {
        Cow::Owned(A {})
    }
}

impl<'a> Add<&A> for &'a A {
    type Output = A;
    fn add(self, _: &A) -> Self::Output {
        A {}
    }
}

impl<'a> LCM<&A> for &'a A {
    type Output = A;
    fn lcm(&self) -> Cow<Self::Output> {
        Cow::Owned(A {})
    }
}
trait MyAdd<T> {
    type Output;
    fn my_add(self, other: T) -> Self::Output;
}

impl<T, U, Out> MyAdd<U> for T
where
    T: LCM<U, Output = Out>,
    U: LCM<T, Output = Out>,
    for<'a, 'b> &'a Out: Add<&'b Out, Output = Out>,
    Out: Clone,
{
    type Output = Out;
    fn my_add(self, other: U) -> Self::Output {
        self.lcm().as_ref() + other.lcm().as_ref()
    }
}

// struct Contain<T> {
//     value: T,
// }

// impl<'a, T> Add<&Contain<T>> for &'a Contain<T>
// where
//     for<'b, 'c> &'b T: Add<&'c T, Output = T>,
// {
//     type Output = Contain<T>;
//     fn add(self, other: &Contain<T>) -> Contain<T> {
//         Contain {
//             value: &self.value + &other.value,
//         }
//     }
// }

fn main() {
    let a = &A {};
    let b = B {};
    let c = a.my_add(a);
}
