use std::any::Any;

struct A {}

struct B {}

struct C<T> {
    t: T,
}

trait Contract<T> {
    type Out;
    fn contract(&self, other: &T) -> Option<Self::Out>;
}

// impl<T, U, Out> Contract<C<T>> for C<U>
// where
//     Out: Default,
// {
//     type Out = C<Out>;
//     fn contract(&self, other: &C<T>) -> Option<Self::Out> {
//         Some(C {
//             t: Default::default(),
//         })
//     }
// }

// impl DynContract for A {
//     fn dyn_contract(&self, other: &dyn Any) -> Option<Box<dyn Any>> {
//         if let Some(other_b) = other.downcast_ref::<B>() {
//             self.contract(other_b).map(Box::new)
//         } else {
//             None
//         }
//     }
// }

// impl DynContract for B {
//     fn dyn_contract(&self, other: &dyn Any) -> Option<Box<dyn Any>> {
//         if let Some(other_a) = other.downcast_ref::<A>() {
//             self.contract(other_a).map(Box::new)
//         } else {
//             None
//         }
//     }
// }

fn main() {
    // let elements: Vec<Box<dyn DynContract>> = vec![Box::new(A {}), Box::new(B {})];

    // Example of using the dyn_contract method
    // Note: actual usage depends on the logic of how you want to pair elements.
    // if let Some(result) = elements[0].dyn_contract(&*elements[1]) {
    // Do something with the result
    // }
}
