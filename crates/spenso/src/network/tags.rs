use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    id::{Condition, PatternRestriction},
    printer::PrintState,
    symbol, tag,
};

use crate::{
    shadowing::symbolica_utils::SpensoPrintSettings, structure::abstract_index::AIND_SYMBOLS,
};

pub struct SpensoTags {
    pub tag: String,
    pub rank1: String,
    pub rank1_: Symbol,
    pub upper: String,
    pub lower: String,
    pub bracket: Symbol,
    pub pure_scalar: Symbol,
    pub tensor: String,
    pub tensor_: Symbol,
    pub index: String,
    pub representation: String,
    pub i_: Symbol,
    pub dot: Symbol,
    pub rep_: Symbol,
    pub self_dual: String,
    pub self_dual_: Symbol,
    pub dualizable: String,
    pub dualizable_: Symbol,
}

pub static SPENSO_TAG: std::sync::LazyLock<SpensoTags> = std::sync::LazyLock::new(SpensoTags::new);

impl SpensoTags {
    fn new() -> Self {
        let tag = tag!("broadcast");
        let upper = tag!("upper");
        let lower = tag!("lower");
        let rank1 = tag!("rank1");
        let tensor = tag!("tensor");
        let index = tag!("index");
        let representation = tag!("representation");
        let self_dual = tag!("self_dual");
        let dualizable = tag!("dualizable");
        Self {
            rank1_: symbol!("rank1_", tags = [&tensor, &rank1]),
            bracket: symbol!("bracket"),
            pure_scalar: symbol!("pure_scalar"),
            dot: symbol!("dot";Symmetric,Linear; print = |a,opt|{
                    match opt.custom_print_mode {
                        Some(("spenso",i))=>{
                            let SpensoPrintSettings{
                                parens,
                                with_dim,..
                            } = SpensoPrintSettings::from(i);


                    let AtomView::Fun(f) = a else {
                        return None;
                    };

                    if f.get_nargs() != 3 {
                        return None;
                    }
                    let mut args = f.iter();

                    let a = args.next().unwrap();
                    let b = args.next().unwrap();
                    let c = args.next().unwrap();

                    fn is_rep(view:AtomView<'_>)->bool{
                        match view {
                            AtomView::Fun(f) if f.get_symbol().has_tag(&SPENSO_TAG.upper) => true,
                            AtomView::Var(s) if s.get_symbol().has_tag(&SPENSO_TAG.upper) => true,
                            _=>false
                        }
                    }

                    let (a,b,c) = if is_rep(a) && !is_rep(b) && !is_rep(c) {
                        (a,b,c)
                    } else if is_rep(b) && !is_rep(a) && !is_rep(c) {
                        (b,c,a)
                    } else if is_rep(c) && !is_rep(a) && !is_rep(b) {
                        (c,a,b)
                    } else { return None};

                    let mut s = String::new();
                    if parens {
                        s.push('(');
                    }
                    b.format(&mut s, opt,PrintState::new()).unwrap();
                    s.push('.');
                    if with_dim {a.format(&mut s, opt, PrintState::new()).unwrap();
                        s.push('.');
                    }
                    c.format(&mut s, opt,PrintState::new()).unwrap();
                    if parens {
                        s.push(')');
                    }
                    Some(s)

                },
                _=>None
            }



            }),
            tensor_: symbol!("tensor_", tag = tensor),
            i_: symbol!("i_", tag = &index),
            rep_: symbol!("rep_", tag = &representation),
            self_dual_: symbol!("self_dual_", tags = [&representation, &self_dual]),
            dualizable_: symbol!("dualizable_", tags = [&representation, &dualizable]),
            tag,
            upper,
            lower,
            rank1,
            tensor,
            index,
            representation,
            self_dual,
            dualizable,
        }
    }

    pub fn rank1_<'a, const N: usize, A: Into<AtomOrView<'a>>>(
        &self,
        args: impl IntoIterator<Item = A>,
    ) -> Atom {
        let symbol = if N == 0 {
            self.rank1_
        } else {
            self.rank_one_tensor_symbol(&format!("rank1{}_", N))
        };
        let mut f = FunctionBuilder::new(symbol);
        let mut has_args = false;
        for a in args.into_iter() {
            has_args = true;
            f = f.add_arg(a);
        }
        if has_args {
            f.finish()
        } else {
            Atom::var(symbol)
        }
    }

    pub fn index_fiter(&self, symbol: Symbol) -> Condition<PatternRestriction> {
        symbol.filter_tag(self.index.clone()) | symbol.filter_single(|a| a.is_integer())
    }
    pub fn rep_<'a, const N: usize, A: Into<AtomOrView<'a>>>(
        &self,
        args: impl IntoIterator<Item = A>,
    ) -> Atom {
        let symbol = if N == 0 {
            self.rep_
        } else {
            self.representation_symbol(&format!("rep{}_", N))
        };
        let mut f = FunctionBuilder::new(symbol);
        let mut has_args = false;
        for a in args.into_iter() {
            has_args = true;
            f = f.add_arg(a);
        }
        if has_args {
            f.finish()
        } else {
            Atom::var(symbol)
        }
    }

    pub fn tensor_symbol(&self, name: &str) -> Symbol {
        symbol!(name, tag = &self.tensor)
    }

    pub fn representation_symbol(&self, name: &str) -> Symbol {
        symbol!(name, tag = &self.representation)
    }

    pub fn self_dual_symbol(&self, name: &str) -> Symbol {
        symbol!(name, tags = [&self.representation, &self.self_dual])
    }

    pub fn dualizable_symbol(&self, name: &str) -> Symbol {
        symbol!(name, tags = [&self.representation, &self.dualizable])
    }

    pub fn rank_one_tensor_symbol(&self, name: &str) -> Symbol {
        symbol!(name, tags = [&self.tensor, &self.rank1])
    }

    pub fn tensor_<'a, const N: usize, A: Into<AtomOrView<'a>>>(
        &self,
        args: impl IntoIterator<Item = A>,
    ) -> Atom {
        let symbol = if N == 0 {
            self.tensor_
        } else {
            self.tensor_symbol(&format!("tensor{}_", N))
        };
        let mut f = FunctionBuilder::new(symbol);
        let mut has_args = false;
        for a in args.into_iter() {
            has_args = true;
            f = f.add_arg(a);
        }
        if has_args {
            f.finish()
        } else {
            Atom::var(symbol)
        }
    }

    pub fn self_dual_<'a, const N: usize, A: Into<AtomOrView<'a>>>(
        &self,
        args: impl IntoIterator<Item = A>,
    ) -> Atom {
        let symbol = if N == 0 {
            self.self_dual_
        } else {
            self.self_dual_symbol(&format!("self_dual{}_", N))
        };
        let mut f = FunctionBuilder::new(symbol);
        let mut has_args = false;
        for a in args.into_iter() {
            has_args = true;
            f = f.add_arg(a);
        }
        if has_args {
            f.finish()
        } else {
            Atom::var(symbol)
        }
    }

    pub fn dualizable_<'a, const N: usize, A: Into<AtomOrView<'a>>>(
        &self,
        args: impl IntoIterator<Item = A>,
    ) -> Atom {
        let symbol = if N == 0 {
            self.dualizable_
        } else {
            self.dualizable_symbol(&format!("self_dual{}_", N))
        };
        let mut f = FunctionBuilder::new(symbol);
        let mut has_args = false;
        for a in args.into_iter() {
            has_args = true;
            f = f.add_arg(a);
        }
        if has_args {
            f.finish()
        } else {
            Atom::var(symbol)
        }
    }

    pub fn dualizable_dual_<'a, const N: usize, A: Into<AtomOrView<'a>>>(
        &self,
        args: impl IntoIterator<Item = A>,
    ) -> Atom {
        AIND_SYMBOLS.dual(self.dualizable_::<N, A>(args))
    }
}
