struct A {}
struct B {}
pub trait Upgradable {}

impl Upgradable for A {}
impl Upgradable for B {}

trait SmallestUpgrade<T: Upgradable> {
    type Output: Upgradable;
}

// impl SmallestUpgrade<A> for B {
//     type Output = B;
// }

impl<T: Upgradable> SmallestUpgrade<T> for T {
    type Output = T;
}

impl SmallestUpgrade<B> for A {
    type Output = B;
}

impl SmallestUpgrade<A> for B {
    type Output = B;
}

struct Foo<T> {
    ts: Vec<T>,
}

pub trait Mult<T: Upgradable> {
    fn mult<U: SmallestUpgrade<T>>(&self, other: &Foo<U>) -> Option<Foo<U::Output>>;
}

impl<T: Upgradable> Mult<T> for Foo<T> {
    fn mult<U: SmallestUpgrade<T>>(&self, other: &Foo<U>) -> Option<Foo<U::Output>> {
        None
    }
}

fn main() {
    let a = Foo { ts: vec![A {}] };
    let b = Foo { ts: vec![B {}] };
    let ap = Foo { ts: vec![A {}] };
    let bp = Foo { ts: vec![B {}] };

    let c = b.mult(&a);
    let d = a.mult(&b);
    let e = ap.mult(&a);
    let f = bp.mult(&b);
}
