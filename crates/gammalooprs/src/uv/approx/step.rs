// use symbolica::atom::{Atom, AtomCore};

// use crate::uv::approx::{ForestNodeLike, UVCtx};

// pub struct UvStep<'a, S> {
//     pub ctx: &'a UVCtx<'a>,
//     pub current: &'a S,
//     pub given: &'a S,
// }

// pub struct Dim<const N: u8>;

// pub trait Scheme<S: ForestNodeLike> {
//     fn integrate(step: &UvStep<'_, S>, integrand: impl AtomCore) -> Atom;

//     fn full(step: &UvStep<'_, S>, integrand: impl AtomCore) -> Atom;
// }
// pub trait Approximate<S: ForestNodeLike> {
//     fn local(step: &UvStep<'_, S>, integrand: impl AtomCore) -> Atom;
//     fn full(step: &UvStep<'_, S>, integrand: impl AtomCore) -> Atom;
// }

// impl<S: ForestNodeLike> Approximate<S> for Dim<3> {
//     fn local(step: &UvStep<'_, S>, integrand: impl AtomCore) -> Atom {
//         todo!()
//     }

//     fn full(step: &UvStep<'_, S>, integrand: impl AtomCore) -> Atom {
//         todo!()
//     }

//     fn final(step: &UvStep<'_, S>, integrand: impl AtomCore) -> Atom {
//         todo!()
//     }
// }

// impl<S: ForestNodeLike> UvStep<'_, S> {
//     pub fn new(ctx: &UVCtx<'_>, current: &S, given: &S) -> Self {
//         Self { ctx, current, given }
//     }

//     pub fn local<F: Approximate<S>>(step: &UvStep<'_, S>, integrand: impl AtomCore) -> Atom {
//         <F>::local(step, integrand)
//     }

//     pub fn full<F: Approximate<S>>(step: &UvStep<'_, S>, integrand: impl AtomCore) -> Atom {
//         <F>::full(step, integrand)
//     }
// }
