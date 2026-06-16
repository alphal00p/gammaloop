use crate::structure::concrete_index::ConcreteIndex;
use ::symbolica_utils::IntoSymbol;
use ahash::HashMap;
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder, Symbol},
    printer::{PrintOptions, PrintUserData},
    symbol,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct SpensoPrintSettings {
    pub with_dim: bool,
    pub parens: bool,
    pub commas: bool,
    pub index_subscripts: bool,
    pub symbol_scripts: bool,
}

impl From<SpensoPrintSettings> for HashMap<String, PrintUserData> {
    fn from(settings: SpensoPrintSettings) -> Self {
        HashMap::from_iter([(
            "spenso".to_string(),
            PrintUserData::Integer(settings.to_usize() as i64),
        )])
    }
}

impl From<&SpensoPrintSettings> for HashMap<String, PrintUserData> {
    fn from(settings: &SpensoPrintSettings) -> Self {
        (*settings).into()
    }
}

impl From<usize> for SpensoPrintSettings {
    fn from(flag: usize) -> Self {
        Self::from_usize(flag)
    }
}

impl SpensoPrintSettings {
    fn from_usize(x: usize) -> Self {
        Self {
            parens: (x & 0b00001) != 0,
            commas: (x & 0b00010) != 0,
            with_dim: (x & 0b00100) != 0,
            symbol_scripts: (x & 0b01000) != 0,
            index_subscripts: (x & 0b10000) != 0,
        }
    }

    fn to_usize(self) -> usize {
        (self.parens as usize)
            | ((self.commas as usize) << 1)
            | ((self.with_dim as usize) << 2)
            | ((self.symbol_scripts as usize) << 3)
            | ((self.index_subscripts as usize) << 4)
    }

    pub fn typst() -> Self {
        Self {
            parens: true,
            commas: false,
            with_dim: false,
            symbol_scripts: true,
            index_subscripts: true,
        }
    }

    pub fn is_typst(&self) -> bool {
        self == &Self::typst()
    }

    pub fn compact() -> Self {
        Self {
            parens: true,
            commas: false,
            with_dim: false,
            symbol_scripts: false,
            index_subscripts: false,
        }
    }

    pub fn nice_symbolica(&self) -> PrintOptions {
        PrintOptions {
            custom_print_mode: self.into(),
            color_builtin_symbols: true,
            terms_on_new_line: true,
            max_line_length: Some(120),
            color_namespace: false,
            multiplication_operator: '·',
            hide_all_namespaces: true,
            color_top_level_sum: true,
            num_exp_as_superscript: true,
            ..Default::default()
        }
    }

    pub fn typst_symbolica(&self) -> PrintOptions {
        PrintOptions {
            custom_print_mode: self.into(),
            color_builtin_symbols: false,
            terms_on_new_line: false,
            color_namespace: false,
            multiplication_operator: ' ',
            hide_all_namespaces: true,
            color_top_level_sum: false,
            num_exp_as_superscript: true,
            ..Default::default()
        }
    }
}

pub trait LogPrint {
    fn log_print(&self, max_line_length: Option<usize>) -> String;
}

impl<A: AtomCore> LogPrint for A {
    fn log_print(&self, max: Option<usize>) -> String {
        let mut settings = SpensoPrintSettings::compact().nice_symbolica();
        settings.max_line_length = max;
        self.printer(settings).to_string()
    }
}

pub fn atomic_expanded_label<I: IntoSymbol>(indices: &[ConcreteIndex], name: I) -> Atom {
    let id = name.ref_into_symbol();
    atomic_expanded_label_id(indices, id, &[])
}

pub fn atomic_flat_label<I: IntoSymbol>(index: usize, name: I) -> Atom {
    let id = name.ref_into_symbol();
    atomic_flat_label_id(index, id)
}

#[allow(clippy::cast_possible_wrap)]
pub fn atomic_flat_label_id(index: usize, id: Symbol) -> Atom {
    FunctionBuilder::new(id)
        .add_arg(Atom::num(index as i64).as_view())
        .finish()
}

#[allow(clippy::cast_possible_wrap)]
pub fn atomic_expanded_label_id(indices: &[ConcreteIndex], name: Symbol, args: &[Atom]) -> Atom {
    let mut value_builder = FunctionBuilder::new(name);
    let mut index_func = FunctionBuilder::new(symbol!("cind"));
    for arg in args {
        value_builder = value_builder.add_arg(arg);
    }
    for &index in indices {
        index_func = index_func.add_arg(Atom::num(index as i64).as_view());
    }

    value_builder.add_arg(index_func.finish()).finish()
}

#[cfg(test)]
mod test {
    use crate::network::tags::SPENSO_TAG;
    use symbolica::{
        atom::AtomCore,
        parse,
        printer::{PrintOptions, PrintUserData},
        symbol, tag,
    };

    #[test]
    fn print() {
        let _lower = symbol!(
            "lower",
            tag = SPENSO_TAG.lower,
            print = |_, opt, _state| {
                if matches!(
                    opt.custom_print_mode.get("typst"),
                    Some(PrintUserData::Integer(1))
                ) {
                    let body = r#"{
let args = arg.pos().map(to-eq).join("")
(content: args,lower:true)
}"#;
                    Some(body.into())
                } else {
                    None
                }
            }
        );

        let _upper = symbol!(
            "upper",
            tags = [SPENSO_TAG.upper.clone(), tag!("Real")],
            print = |_, opt, _state| {
                if matches!(
                    opt.custom_print_mode.get("typst"),
                    Some(PrintUserData::Integer(1))
                ) {
                    let body = r#"{
let args = arg.pos().map(to-eq).join("")
(content: args,upper:true)
}"#;
                    Some(body.into())
                } else {
                    None
                }
            }
        );

        let expr = parse!(
            "a*f(lower(f(lower(upper(x),lower(a)),lower(y,c))))^(sin(x)*cos(x))*g(x,lower(y),upper(x+1),lower(1))/(x+1)/sin(g(y))*smth(3*x)^(-m)"
        );

        println!("{}", expr.printer(PrintOptions::typst()))
    }
}
