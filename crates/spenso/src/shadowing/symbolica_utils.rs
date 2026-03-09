use crate::{
    network::parsing::SPENSO_TAG, structure::concrete_index::ConcreteIndex,
    tensors::parametric::atomcore::PatternReplacement,
};
use derive_more::Display;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{
        AddView, Atom, AtomCore, AtomView, FunctionBuilder, MulView, PowView, Symbol,
        representation::FunView,
    },
    coefficient::CoefficientView,
    domains::{float::Complex, rational::Rational},
    evaluate::{FunctionMap, Instruction, OptimizationSettings, Slot},
    function,
    id::Context,
    printer::{CanonicalOrderingSettings, PrintOptions, PrintState},
    state::State,
    symbol,
    utils::Settable,
};

extern crate derive_more;

use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::{Debug, Display, Error},
};

use eyre::Result;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct SpensoPrintSettings {
    pub with_dim: bool,
    pub parens: bool,
    pub commas: bool,
    pub index_subscripts: bool,
    pub symbol_scripts: bool,
}

impl From<SpensoPrintSettings> for (&'static str, usize) {
    fn from(settings: SpensoPrintSettings) -> Self {
        ("spenso", settings.to_usize())
    }
}

impl From<&SpensoPrintSettings> for (&'static str, usize) {
    fn from(settings: &SpensoPrintSettings) -> Self {
        ("spenso", settings.to_usize())
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

    // x^-2*a^-2*b^-2 -> 1/(x * a * b)^2
    // x^-1*a^-1*b^-1 -> ((1/x)/a)/b
    pub fn compact() -> Self {
        Self {
            parens: false,
            commas: false,
            with_dim: false,
            symbol_scripts: false,
            index_subscripts: false,
        }
    }

    pub fn nice_symbolica(&self) -> PrintOptions {
        PrintOptions {
            custom_print_mode: Some(self.into()),
            color_builtin_symbols: true,
            terms_on_new_line: true,
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
            custom_print_mode: Some(self.into()),
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

pub trait AtomCoreExt {
    fn to_bare_ordered_string(&self) -> String;

    fn typst(&self) -> String;

    fn typst_fmt<W: std::fmt::Write>(
        &self,
        f: &mut W,
        settings: &TypstSettings,
    ) -> Result<(), Error>;

    fn is_upper(&self) -> bool;
    fn is_lower(&self) -> bool;
}
pub struct TypstSettings {
    pub preamble: String,
    pub default_dis: String,
    pub default_sym: String,
}

impl TypstSettings {
    pub fn addln_preamble(&mut self, line: &str) {
        self.preamble.push('\n');
        self.preamble.push_str(line);
    }

    pub fn lowering() -> Self {
        TypstSettings {
            default_dis: r#"
let args = ()
let uppers = ()
let lowers = ()
for a in arg.pos(){
  if type(a)==content{
    args.push(a)
  } else if type(a)==dictionary{
    if a.at("upper",default:false){
      uppers.push(to-eq(a.content))
      lowers.push(hide(to-eq(a.content)))
    }else if a.at("lower",default:false){
      lowers.push(to-eq(a.content))
      uppers.push(hide(to-eq(a.content)))
    } else{
      args.push(to-eq(a.content))
    }
  } else{
    args.push(to-eq(a))
  }
}
let arg = args.join()
if uppers.len()!=0{
  let upper= uppers.join()
  let lower= lowers.join()
  if args.len()==0{
    $attach(op(#name),t:#upper,b:#lower)$
  }else{
    $attach(op(#name),t:#upper,b:#lower)(#arg)$
  }
}else{
  $op(#name)(#arg)$
}
"#
            .into(),
            ..Default::default()
        }
    }
}
impl Default for TypstSettings {
    fn default() -> Self {
        let default_dis = "$ op(#name)(#arg.pos().map(to_eq).join(\", \")) $";
        let default_sym = "$ #name $";

        TypstSettings {
            preamble: "let to-eq(a)=if type(a)==dictionary{a.content}else {$#a$}".into(),
            default_dis: default_dis.into(),
            default_sym: default_sym.into(),
        }
    }
}

impl<A: AtomCore> AtomCoreExt for A {
    fn to_bare_ordered_string(&self) -> String {
        self.to_canonically_ordered_string(CanonicalOrderingSettings {
            include_namespace: false,
            include_attributes: false,
            hide_namespace: None,
        })
    }

    fn typst_fmt<W: std::fmt::Write>(
        &self,
        f: &mut W,
        settings: &TypstSettings,
    ) -> Result<(), Error> {
        let mut params = BTreeMap::new();
        let mut fn_map = FunctionMap::new();
        let mut externals = BTreeSet::new();

        self.visitor(&mut |a| {
            if let AtomView::Var(a) = a {
                params.insert(a.get_symbol(), a.as_view().to_owned());
                false
            } else if let AtomView::Fun(a) = a {
                externals.insert(a.get_symbol());

                let dashed_name = format!(
                    "{}-{}",
                    a.get_symbol().get_namespace(),
                    a.get_symbol().get_stripped_name(),
                );
                let _ = fn_map.add_external_function(a.get_symbol(), dashed_name);
                true
            } else {
                true
            }
        });

        let (params, symbols): (Vec<Atom>, Vec<Symbol>) =
            params.into_iter().map(|(s, a)| (a, s)).collect();
        let eval_tree = self
            .evaluator(
                &fn_map,
                &params,
                OptimizationSettings {
                    horner_iterations: 10,
                    ..Default::default()
                },
            )
            .unwrap();

        writeln!(f, "#{{")?;
        writeln!(f, "{}", &settings.preamble)?;
        fn typst_rat(r: &Rational) -> String {
            if r.is_integer() {
                r.numerator().to_string()
            } else {
                format!(" {}/{}", r.numerator(), r.denominator())
            }
        }

        fn typst_slot(s: Slot, consts: &[Complex<Rational>]) -> String {
            match s {
                Slot::Const(c) => {
                    let Complex { re, im } = &consts[c];

                    match (re.is_zero(), im.is_zero()) {
                        (true, true) => "0".into(),
                        (true, false) => format!("i {}", typst_rat(im)),
                        (false, true) => format!(" {}", typst_rat(re)),
                        _ => format!("({} + i {})", typst_rat(re), typst_rat(im)),
                    }
                }
                Slot::Out(c) => format!("out{c}"),
                Slot::Param(c) => format!("param{c}"),
                Slot::Temp(c) => format!("tmp{c}"),
            }
        }

        let (instr, _, consts) = eval_tree.export_instructions();

        writeln!(
            f,
            "let default_dis(name:\"\",namespace:\"\",..arg)={{{}}}",
            &settings.default_dis
        )?;

        writeln!(
            f,
            "let default_sym(name:\"\",namespace:\"\")={{{}}}",
            &settings.default_sym
        )?;

        for ((i, s), a) in symbols.iter().enumerate().zip(params) {
            write!(f, "let param{i} = ")?;
            if let Some(p) = s.get_print_function()
                && let Some(a) = p(
                    a.as_view(),
                    &PrintOptions {
                        custom_print_mode: Some(("typst", 1)),
                        ..Default::default()
                    },
                )
            {
                writeln!(f, "{a}")?;
            } else {
                let name = s.get_stripped_name();
                let namespace = s.get_namespace();
                writeln!(f, "default_sym(namespace:\"{namespace}\",name:\"{name}\")",)?;
            }
        }

        for s in &externals {
            if s.is_builtin() {
                continue;
            }
            let name = s.get_stripped_name();
            let namespace = s.get_namespace();
            let atom = function!(*s, symbol!("args"));
            write!(f, "let {namespace}-{name}")?;
            if let Some(p) = s.get_print_function()
                && let Some(a) = p(
                    atom.as_view(),
                    &PrintOptions {
                        custom_print_mode: Some(("typst", 1)),
                        ..Default::default()
                    },
                )
            {
                writeln!(f, "{a}")?;
            } else {
                let name = s.get_stripped_name();
                let namespace = s.get_namespace();
                writeln!(
                    f,
                    "(..arg) = default_dis(namespace:\"{namespace}\",name:\"{name}\",..arg)"
                )?;
            }
        }

        for i in instr {
            match i {
                Instruction::Add(s, args, _is_real) => {
                    writeln!(
                        f,
                        "let {} = $({})$",
                        typst_slot(s, &consts),
                        args.into_iter().map(|a| typst_slot(a, &consts)).join(" + ")
                    )?;
                }
                Instruction::Mul(s, args, _is_real) => {
                    writeln!(
                        f,
                        "let {} = $({})$",
                        typst_slot(s, &consts),
                        args.into_iter()
                            .map(|a| typst_slot(a, &consts))
                            .join(" dot ")
                    )?;
                }
                Instruction::ExternalFun(o, name, args) => {
                    writeln!(
                        f,
                        "let {} = {name}({})",
                        typst_slot(o, &consts),
                        args.into_iter().map(|a| typst_slot(a, &consts)).join(",")
                    )?;
                }
                Instruction::Fun(o, builtin, s, _is_real) => {
                    let b = builtin.get_symbol();

                    if b == Symbol::COS {
                        writeln!(
                            f,
                            "let {} = $ cos({})$",
                            typst_slot(o, &consts),
                            typst_slot(s, &consts)
                        )?;
                    } else if b == Symbol::SIN {
                        writeln!(
                            f,
                            "let {} = $ sin({})$",
                            typst_slot(o, &consts),
                            typst_slot(s, &consts)
                        )?;
                    } else if b == Symbol::SQRT {
                        writeln!(
                            f,
                            "let {} = $ sqrt({})$",
                            typst_slot(o, &consts),
                            typst_slot(s, &consts)
                        )?;
                    } else {
                        let name = b.get_stripped_name().to_string();
                        writeln!(
                            f,
                            "let {} = $op(\"{name}\")({})$",
                            typst_slot(o, &consts),
                            typst_slot(s, &consts)
                        )?;
                    }
                }
                Instruction::Powf(o, b, e, _is_real) => {
                    writeln!(
                        f,
                        "let {} = ${}^({})$",
                        typst_slot(o, &consts),
                        typst_slot(b, &consts),
                        typst_slot(e, &consts)
                    )?;
                }
                Instruction::Pow(o, b, e, _is_real) => {
                    writeln!(
                        f,
                        "let {} = ${}^({})$",
                        typst_slot(o, &consts),
                        typst_slot(b, &consts),
                        e,
                    )?;
                }
                _ => {
                    println!("{i:?}")
                }
            }
        }
        writeln!(f, "out0")?;
        writeln!(f, "}}")
    }

    fn typst(&self) -> String {
        let mut out = String::new();
        self.as_atom_view()
            .fmt_output(
                &mut out,
                &PrintOptions {
                    custom_print_mode: Some(("typst", 2)),
                    ..Default::default()
                },
                PrintState::new(),
            )
            .unwrap();
        out
    }

    fn is_upper(&self) -> bool {
        match self.as_atom_view() {
            AtomView::Fun(a) => a.get_symbol().has_tag(&SPENSO_TAG.upper),
            AtomView::Var(a) => a.get_symbol().has_tag(&SPENSO_TAG.upper),
            _ => false,
        }
    }

    fn is_lower(&self) -> bool {
        match self.as_atom_view() {
            AtomView::Fun(a) => a.get_symbol().has_tag(&SPENSO_TAG.lower),
            AtomView::Var(a) => a.get_symbol().has_tag(&SPENSO_TAG.lower),
            _ => false,
        }
    }
}
pub struct Typst;

pub trait FormatWithState {
    fn fmt_output<W: std::fmt::Write>(
        &self,
        f: &mut W,
        opts: &PrintOptions,
        print_state: PrintState,
    ) -> Result<bool, Error>;
}

impl FormatWithState for AtomView<'_> {
    fn fmt_output<W: std::fmt::Write>(
        &self,
        fmt: &mut W,
        opts: &PrintOptions,
        print_state: PrintState,
    ) -> Result<bool, Error> {
        match self {
            AtomView::Num(n) => n.as_view().format(fmt, opts, print_state),
            AtomView::Var(v) => v.as_view().format(fmt, opts, print_state),
            AtomView::Fun(f) => f.fmt_output(fmt, opts, print_state),
            AtomView::Pow(p) => p.fmt_output(fmt, opts, print_state),
            AtomView::Mul(t) => t.fmt_output(fmt, opts, print_state),
            AtomView::Add(e) => e.fmt_output(fmt, opts, print_state),
        }
    }
}

impl FormatWithState for FunView<'_> {
    fn fmt_output<W: std::fmt::Write>(
        &self,
        f: &mut W,
        opts: &PrintOptions,
        print_state: PrintState,
    ) -> Result<bool, Error> {
        if print_state.in_sum {
            f.write_char('+')?;
        }

        // let id = self.get_symbol();

        // if let Some(custom_print) = &id.fo {
        //     if let Some(s) = custom_print(self.as_view(), opts) {
        //         f.write_str(&s)?;
        //         return Ok(false);
        //     }
        // }
        let mut uppers = vec![];
        let mut lowers = vec![];
        // f.write_str("attach(")?;
        for a in self.iter() {
            if a.is_upper() {
                let mut out = String::new();
                a.fmt_output(&mut out, opts, print_state)?;
                let out_low = format!("#hide(${}$)", out);
                uppers.push(out);
                lowers.push(out_low);
            }

            if a.is_lower() {
                let mut out = String::new();
                a.fmt_output(&mut out, opts, print_state)?;
                let out_low = format!("#hide(${}$)", out);
                lowers.push(out);
                uppers.push(out_low);
            }
        }

        if uppers.is_empty() {
            f.write_str("op(\"")?;
            self.get_symbol().format(opts, f)?;
            f.write_str("\")")?;
            let n_args = self.get_nargs();

            if n_args > 0 {
                f.write_char('(')?;
            }
            for (i, a) in self.iter().enumerate() {
                if i + 1 < n_args {
                    f.write_char(',')?;
                }
                a.fmt_output(f, opts, print_state)?;
            }
            if n_args > 0 {
                f.write_char(')')?;
            }
        } else {
            f.write_str("scripts(attach(")?;
            f.write_str("op(\"")?;

            self.get_symbol().format(opts, f)?;
            f.write_str("\")")?;
            f.write_str(", tr: ")?;
            f.write_str(&uppers.join(" "))?;
            f.write_str(", br: ")?;
            f.write_str(&lowers.join(" "))?;
            f.write_char(')')?;
            f.write_char(')')?;
            let n_args = self.get_nargs() - lowers.len();

            if n_args > 0 {
                f.write_char('(')?;
            }

            for (i, a) in self
                .iter()
                .filter(|a| !(a.is_lower() || a.is_upper()))
                .enumerate()
            {
                if i + 1 < n_args {
                    f.write_char(',')?;
                }
                a.fmt_output(f, opts, print_state)?;
            }
            // f.write_char(')');
            if n_args > 0 {
                f.write_char(')')?;
            }
        }
        Ok(false)
    }
}

impl FormatWithState for MulView<'_> {
    fn fmt_output<W: std::fmt::Write>(
        &self,
        f: &mut W,
        opts: &PrintOptions,
        mut print_state: PrintState,
    ) -> Result<bool, Error> {
        let add_paren = print_state.in_exp || print_state.in_exp_base;
        if add_paren {
            if print_state.in_sum {
                print_state.in_sum = false;
                f.write_char('+')?;
            }

            f.write_char('(')?;
            print_state.in_exp = false;
            print_state.in_exp_base = false;
        }

        print_state.in_product = true;

        // write the coefficient first
        let mut first = true;
        let mut skip_num = false;
        if let Some(AtomView::Num(n)) = self.iter().last() {
            print_state.suppress_one = true;
            first = n.as_view().format(f, opts, print_state)?;
            print_state.suppress_one = false;
            skip_num = true;
        } else if print_state.in_sum {
            f.write_char('+')?;
        }

        print_state.top_level_add_child = false;
        print_state.level += 1;
        print_state.in_sum = false;

        for x in self.iter().take(if skip_num {
            self.get_nargs() - 1
        } else {
            self.get_nargs()
        }) {
            if !first {
                f.write_char(' ')?;
            }
            first = false;

            x.fmt_output(f, opts, print_state)?;
        }

        if add_paren {
            f.write_char(')')?;
        }
        Ok(false)
    }
}

impl FormatWithState for PowView<'_> {
    fn fmt_output<W: std::fmt::Write>(
        &self,
        f: &mut W,
        opts: &PrintOptions,
        mut print_state: PrintState,
    ) -> Result<bool, Error> {
        if print_state.in_sum {
            f.write_char('+')?;
        }

        let add_paren = print_state.in_exp_base; // right associative
        if add_paren {
            f.write_char('(')?;
            print_state.in_exp = false;
            print_state.in_exp_base = false;
        }

        let b = self.get_base();
        let e = self.get_exp();

        print_state.top_level_add_child = false;
        print_state.level += 1;
        print_state.in_sum = false;
        print_state.in_product = false;
        print_state.suppress_one = false;

        if let AtomView::Num(n) = e
            && n.get_coeff_view() == CoefficientView::Natural(-1, 1, 0, 1)
        {
            // TODO: construct the numerator
            f.write_str("1/(")?;
            b.fmt_output(f, opts, print_state)?;
            f.write_char(')')?;
            return Ok(false);
        }

        print_state.in_exp_base = true;

        b.fmt_output(f, opts, print_state)?;

        print_state.in_exp_base = false;
        print_state.in_exp = true;

        f.write_char('^')?;

        f.write_char('(')?;
        print_state.in_exp = false;
        e.fmt_output(f, opts, print_state)?;
        f.write_char(')')?;

        if add_paren {
            f.write_char(')')?;
        }

        Ok(false)
    }
}

impl FormatWithState for AddView<'_> {
    fn fmt_output<W: std::fmt::Write>(
        &self,
        f: &mut W,
        opts: &PrintOptions,
        mut print_state: PrintState,
    ) -> Result<bool, Error> {
        let mut first = true;
        print_state.top_level_add_child = print_state.level == 0;
        print_state.level += 1;
        print_state.suppress_one = false;

        let add_paren = print_state.in_product || print_state.in_exp || print_state.in_exp_base;
        if add_paren {
            if print_state.in_sum {
                f.write_char('+')?;
            }

            print_state.in_sum = false;
            print_state.in_product = false;
            print_state.in_exp = false;
            print_state.in_exp_base = false;

            f.write_char('(')?;
        }

        let mut count = 0;
        for x in self.iter() {
            if !first && print_state.top_level_add_child && opts.terms_on_new_line {
                f.write_char('\n')?;
            }
            first = false;

            x.fmt_output(f, opts, print_state)?;
            print_state.in_sum = true;
            count += 1;
        }

        if opts.max_terms.is_some() && count < self.get_nargs() {
            if print_state.top_level_add_child && opts.terms_on_new_line {
                f.write_char('\n')?;
            }

            f.write_str("+...")?;
        }

        if add_paren {
            f.write_char(')')?;
        }
        Ok(false)
    }
}
// // fn print_fun_view(view:FunView<'_>)-?

// use anyhow::Ok;
use serde::ser::SerializeStruct;

#[derive(
    Debug,
    Copy,
    Clone,
    Ord,
    PartialOrd,
    Eq,
    PartialEq,
    Hash,
    Display,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
    bincode_trait_derive::BorrowDecodeFromDecode,
)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct SerializableSymbol {
    symbol: Symbol,
}

impl SerializableSymbol {
    pub fn get_id(&self) -> u32 {
        self.symbol.get_id()
    }

    pub fn get_name(&self) -> &str {
        self.symbol.get_name()
    }
}

impl Serialize for SerializableSymbol {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.symbol.get_name().serialize(serializer)
    }
}

impl<'d> Deserialize<'d> for SerializableSymbol {
    fn deserialize<D>(deserializer: D) -> Result<SerializableSymbol, D::Error>
    where
        D: serde::Deserializer<'d>,
    {
        let value = String::deserialize(deserializer)?;
        Ok(SerializableSymbol {
            symbol: symbol!(&value),
        })
    }
}

impl From<Symbol> for SerializableSymbol {
    fn from(value: Symbol) -> Self {
        Self { symbol: value }
    }
}

impl From<SerializableSymbol> for Symbol {
    fn from(value: SerializableSymbol) -> Self {
        value.symbol
    }
}

impl From<SerializableSymbol> for u32 {
    fn from(value: SerializableSymbol) -> Self {
        value.symbol.get_id()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct SerializableAtom(pub Atom);

impl PatternReplacement for SerializableAtom {
    fn replace_map_mut<F: Fn(AtomView, &Context, &mut Settable<'_, Atom>)>(&mut self, m: &F) {
        self.0.replace_map_mut(m)
    }

    fn replace_multiple_mut<T: symbolica::id::BorrowReplacement>(&mut self, replacements: &[T]) {
        self.0.replace_multiple_mut(replacements)
    }

    fn replace_multiple_repeat<T: symbolica::id::BorrowReplacement>(
        &self,
        replacements: &[T],
    ) -> Self {
        self.0.replace_multiple_repeat(replacements).into()
    }

    fn replace_multiple_repeat_mut<T: symbolica::id::BorrowReplacement>(
        &mut self,
        replacements: &[T],
    ) {
        self.0.replace_multiple_repeat_mut(replacements)
    }
}

impl Display for SerializableAtom {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Serialize for SerializableAtom {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut state = serializer.serialize_struct("SerializableAtom", 3)?;

        let mut serialized_atom: Vec<u8> = Vec::new();

        self.0.as_view().write(&mut serialized_atom).unwrap();

        state.serialize_field("atom", &serialized_atom)?;

        let mut symbolica_state = Vec::new();

        State::export(&mut symbolica_state).unwrap();

        state.serialize_field("state", &symbolica_state)?;
        state.end()
    }
}

impl<'de> Deserialize<'de> for SerializableAtom {
    fn deserialize<D>(deserializer: D) -> std::result::Result<SerializableAtom, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct SerializableAtomHelper {
            atom: Vec<u8>,
            state: Vec<u8>,
        }

        let helper = SerializableAtomHelper::deserialize(deserializer)?;

        let state = helper.state;

        let map = State::import(&mut state.as_slice(), None).unwrap();

        let atom = Atom::import_with_map(&mut helper.atom.as_slice(), &map).unwrap();

        Ok(SerializableAtom(atom))
    }
}

impl<'a> TryFrom<AtomView<'a>> for SerializableAtom {
    type Error = eyre::Error;

    fn try_from(value: AtomView<'a>) -> Result<Self> {
        Ok(SerializableAtom(value.to_owned()))
    }
}

impl From<Atom> for SerializableAtom {
    fn from(atom: Atom) -> Self {
        SerializableAtom(atom)
    }
}

impl From<SerializableAtom> for Atom {
    fn from(atom: SerializableAtom) -> Self {
        atom.0
    }
}

pub fn atomic_expanded_label<I: IntoSymbol>(indices: &[ConcreteIndex], name: I) -> Atom {
    let id = name.ref_into_symbol();
    atomic_expanded_label_id(indices, id, &[])
}
#[cfg(feature = "shadowing")]
pub fn atomic_flat_label<I: IntoSymbol>(index: usize, name: I) -> Atom {
    let id = name.ref_into_symbol();
    atomic_flat_label_id(index, id)
}

#[allow(clippy::cast_possible_wrap)]
#[cfg(feature = "shadowing")]
pub fn atomic_flat_label_id(index: usize, id: Symbol) -> Atom {
    let mut value_builder = FunctionBuilder::new(id);
    value_builder = value_builder.add_arg(Atom::num(index as i64).as_atom_view());
    value_builder.finish()
}
#[cfg(feature = "shadowing")]
#[allow(clippy::cast_possible_wrap)]
pub fn atomic_expanded_label_id(indices: &[ConcreteIndex], name: Symbol, args: &[Atom]) -> Atom {
    let mut value_builder = FunctionBuilder::new(name);
    let mut index_func = FunctionBuilder::new(symbol!("cind"));
    for arg in args {
        value_builder = value_builder.add_arg(arg);
    }
    for &index in indices {
        index_func = index_func.add_arg(Atom::num(index as i64).as_atom_view());
    }

    let indices = index_func.finish();
    value_builder.add_arg(&indices).finish()
}

#[cfg(feature = "shadowing")]
pub trait IntoSymbol {
    fn ref_into_symbol(&self) -> Symbol;

    fn from_str(s: &str) -> Self;
}

#[cfg(feature = "shadowing")]
pub trait IntoArgs {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom>;
    fn args(&self) -> Vec<Atom> {
        self.ref_into_args().collect()
    }
    fn cooked_name(&self) -> std::string::String;
}

#[cfg(feature = "shadowing")]
impl IntoArgs for usize {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        std::iter::once(Atom::num(*self as i64))
    }
    fn cooked_name(&self) -> std::string::String {
        format!("{self}")
    }
}

#[cfg(feature = "shadowing")]
impl IntoArgs for Vec<SerializableAtom> {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        self.iter().map(|x| x.0.clone())
    }
    fn cooked_name(&self) -> std::string::String {
        let init = "".into();
        self.iter()
            .fold(init, |acc, x| acc + x.to_string().as_str())
    }
}

#[cfg(feature = "shadowing")]
impl IntoArgs for () {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        std::iter::empty()
    }
    fn cooked_name(&self) -> std::string::String {
        "".into()
    }
}

#[derive(Debug, Default, Clone, Copy, Serialize, Deserialize)]
pub struct NoArgs;

impl Display for NoArgs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "")
    }
}

#[cfg(feature = "shadowing")]
impl IntoArgs for NoArgs {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        std::iter::empty()
    }
    fn cooked_name(&self) -> std::string::String {
        "".into()
    }
}

#[cfg(feature = "shadowing")]
impl IntoArgs for Atom {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        std::iter::once(self.clone())
    }

    fn cooked_name(&self) -> std::string::String {
        self.to_string()
    }
}

#[cfg(feature = "shadowing")]
impl IntoArgs for Vec<Atom> {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        self.iter().cloned()
    }

    fn cooked_name(&self) -> std::string::String {
        let init = "".into();
        self.iter()
            .fold(init, |acc, x| acc + x.to_string().as_str())
    }
}

#[cfg(feature = "shadowing")]
impl<const N: usize> IntoArgs for [Atom; N] {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        self.iter().cloned()
    }

    fn cooked_name(&self) -> std::string::String {
        let init = "".into();
        self.iter()
            .fold(init, |acc, x| acc + x.to_string().as_str())
    }
}

// #[cfg(feature = "shadowing")]
// impl IntoSymbol for String {
//     fn ref_into_symbol(&self) -> Symbol {
//         symbol!(self)
//     }

//     fn from_str(s: &str) -> Self {
//         s.into()
//     }
// }

#[cfg(feature = "shadowing")]
impl IntoSymbol for Symbol {
    fn ref_into_symbol(&self) -> Symbol {
        *self
    }

    fn from_str(s: &str) -> Self {
        symbol!(s)
    }
}

#[cfg(feature = "shadowing")]
impl IntoSymbol for SerializableSymbol {
    fn ref_into_symbol(&self) -> Symbol {
        self.symbol
    }

    fn from_str(s: &str) -> Self {
        Self { symbol: symbol!(s) }
    }
}

#[cfg(feature = "shadowing")]
impl IntoSymbol for std::string::String {
    fn ref_into_symbol(&self) -> Symbol {
        symbol!(self)
    }
    fn from_str(s: &str) -> Self {
        s.into()
    }
}

#[cfg(test)]
mod test {
    use crate::{
        network::parsing::SPENSO_TAG,
        shadowing::symbolica_utils::{AtomCoreExt, TypstSettings},
    };

    use symbolica::{parse, symbol, tag};
    #[test]
    fn print() {
        let _lower = symbol!(
            "lower",
            tag = SPENSO_TAG.lower,
            print = |_, opt| {
                if let Some(("typst", 1)) = opt.custom_print_mode {
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
            print = |_, opt| {
                if let Some(("typst", 1)) = opt.custom_print_mode {
                    let body = r#"{
let args = arg.pos().map(to-eq).join("")
(content: args,upper:true)
}"#;
                    Some(body.into())
                } else {
                    None
                }
            } // ; Real
        );

        let expr = parse!(
            "a*f(lower(f(lower(upper(x),lower(a)),lower(y,c))))^(sin(x)*cos(x))*g(x,lower(y),upper(x+1),lower(1))/(x+1)/sin(g(y))*smth(3*x)^(-m)"
        );
        // .r ith(parse_lit!(pow(_x, _y)));

        let mut out = String::new();
        expr.typst_fmt(&mut out, &TypstSettings::lowering())
            .unwrap();

        println!("{}", out);
        println!("{}", expr.typst())
    }
}
