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

#[allow(dead_code)]
struct Foo<T> {
    ts: Vec<T>,
}

#[allow(private_bounds, private_interfaces)]
pub trait Mult<T: Upgradable> {
    fn mult<U: SmallestUpgrade<T>>(&self, other: &Foo<U>) -> Option<Foo<U::Output>>;
}

impl<T: Upgradable> Mult<T> for Foo<T> {
    fn mult<U: SmallestUpgrade<T>>(&self, _other: &Foo<U>) -> Option<Foo<U::Output>> {
        None
    }
}

fn main() {
    let a = Foo { ts: vec![A {}] };
    let b = Foo { ts: vec![B {}] };
    let ap = Foo { ts: vec![A {}] };
    let bp = Foo { ts: vec![B {}] };

    let _c = b.mult(&a);
    let _d = a.mult(&b);
    let _e = ap.mult(&a);
    let _f = bp.mult(&b);
}
