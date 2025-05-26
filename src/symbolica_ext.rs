use symbolica::atom::{Atom, AtomOrView, FunctionBuilder, Symbol};

pub trait CallSymbol<T> {
    fn f(&self, args: T) -> Atom;
}

impl<'a, I> CallSymbol<&'a [I]> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
{
    fn f(&self, args: &'a [I]) -> Atom {
        FunctionBuilder::new(*self).add_args(args).finish()
    }
}

impl<'a, 'b, I, J> CallSymbol<(&'a [I], &'b [J])> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
    &'b J: Into<AtomOrView<'b>>,
{
    fn f(&self, args: (&'a [I], &'b [J])) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.0)
            .add_args(args.1)
            .finish()
    }
}

impl<'a, 'b, I, J, const N: usize> CallSymbol<(&'a [I; N], &'b [J])> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
    &'b J: Into<AtomOrView<'b>>,
{
    fn f(&self, args: (&'a [I; N], &'b [J])) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.0.as_slice())
            .add_args(args.1)
            .finish()
    }
}

impl<I, const N: usize> CallSymbol<[I; N]> for Symbol
where
    for<'a> &'a I: Into<AtomOrView<'a>>,
{
    fn f(&self, args: [I; N]) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.as_slice())
            .finish()
    }
}

impl<'b, I, J, const N: usize> CallSymbol<([I; N], &'b [J])> for Symbol
where
    for<'a> &'a I: Into<AtomOrView<'a>>,
    &'b J: Into<AtomOrView<'b>>,
{
    fn f(&self, args: ([I; N], &'b [J])) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.0.as_slice())
            .add_args(args.1)
            .finish()
    }
}

impl<'a, I, J, K> CallSymbol<(&'a [I], &'a [J], &'a [K])> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
    &'a J: Into<AtomOrView<'a>>,
    &'a K: Into<AtomOrView<'a>>,
{
    fn f(&self, args: (&'a [I], &'a [J], &'a [K])) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.0)
            .add_args(args.1)
            .add_args(args.2)
            .finish()
    }
}
