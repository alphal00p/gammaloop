#[cfg(test)]
use std::sync::LazyLock;
use symbolica_notation::IndexDisplay;

use spenso::network::library::symbolic::ExplicitKey;
use spenso::network::tags;
use spenso::network::tags::{SPENSO_TAG as T, mathml};
use spenso::shadowing::IntoAtom;
use spenso::shadowing::symbolica_utils::{SpensoPrintBackend, SpensoPrintSettings};
use spenso::structure::abstract_index::{AIND_SYMBOLS, AbstractIndex};
use spenso::structure::dimension::Dimension;
use spenso::structure::representation::{RepName, RepresentationMetadata};
use spenso::structure::slot::AbsInd;
use spenso::structure::{Canonicalized, TensorStructure};
use spenso::tensor_symbol;
use symbolica::atom::{
    Atom, AtomCore, AtomOrView, AtomView, EvaluationInfo, FunctionBuilder, Symbol,
};
use symbolica::coefficient::CoefficientView;
use symbolica::domains::rational::Rational;
use symbolica::printer::{PrintOptions, PrintState, PrintUserData};
use symbolica::symbol;

use crate::color::casimir::CofDimensionInvariantRewriter;
use crate::representations::ColorAntiFundamental;
use crate::shorthands::metric::PermuteWithMetric;

use super::rep_symbols::RS;
use super::representations::{ColorAdjoint, ColorFundamental};

mod casimir;
mod conjugate;
mod macros;
pub(crate) mod simplify;

pub use conjugate::color_conj_impl;

#[derive(Debug)]
pub enum ColorError {
    NotFully(Atom),
}

pub struct ColorSymbols {
    pub nc_: Symbol,
    pub adj_: Symbol,
    /// Symbol backing the color fundamental representation function.
    pub fundamental_rep: Symbol,
    /// Symbol backing the color adjoint representation function.
    pub adjoint_rep: Symbol,
    /// The generator symbol
    pub t: Symbol,
    /// The structure constant symbol i.e. [T^a, T^b] = i f^{abc} T^c
    pub f: Symbol,
    /// The symmetric color invariant symbol.
    pub d: Symbol,
    /// The degree-k Gram invariant symbol for two symmetric traces.
    pub gram: Symbol,
    /// The degree-k Casimir eigenvalue symbol.
    pub cas: Symbol,
    /// The degree-k Dynkin index symbol.
    pub idx: Symbol,
    /// The number of colors symbol (i.e. the dimension of the fundamental representation) usually Nc=3
    pub nc: Symbol,
}

impl ColorSymbols {
    pub fn chain_t<'a>(&self, adjoint_index: impl Into<AtomOrView<'a>>) -> Atom {
        FunctionBuilder::new(self.t)
            .add_arg(adjoint_index)
            .add_arg(Atom::var(T.chain_in))
            .add_arg(Atom::var(T.chain_out))
            .finish()
    }

    pub fn explicit_t<'a, 'b, 'c>(
        &self,
        adjoint_index: impl Into<AtomOrView<'a>>,
        left_fundamental: impl Into<AtomOrView<'b>>,
        right_fundamental: impl Into<AtomOrView<'c>>,
    ) -> Atom {
        FunctionBuilder::new(self.t)
            .add_arg(adjoint_index)
            .add_arg(left_fundamental)
            .add_arg(right_fundamental)
            .finish()
    }

    pub fn structure_f<'a, 'b, 'c>(
        &self,
        a: impl Into<AtomOrView<'a>>,
        b: impl Into<AtomOrView<'b>>,
        c: impl Into<AtomOrView<'c>>,
    ) -> Atom {
        FunctionBuilder::new(self.f)
            .add_arg(a)
            .add_arg(b)
            .add_arg(c)
            .finish()
    }

    pub fn symmetric_d<'a>(&self, rep: impl Into<AtomOrView<'a>>, args: Vec<Atom>) -> Atom {
        args.into_iter()
            .fold(FunctionBuilder::new(self.d).add_arg(rep), |d, arg| {
                d.add_arg(arg)
            })
            .finish()
    }

    pub fn gram<D: IntoAtom, L: IntoAtom, R: IntoAtom>(
        &self,
        degree: D,
        left: L,
        right: R,
    ) -> Atom {
        FunctionBuilder::new(self.gram)
            .add_arg(degree.into_atom())
            .add_arg(left.into_atom())
            .add_arg(right.into_atom())
            .finish()
    }

    pub fn cas<D: IntoAtom, R: IntoAtom>(&self, degree: D, rep: R) -> Atom {
        FunctionBuilder::new(self.cas)
            .add_arg(degree.into_atom())
            .add_arg(rep.into_atom())
            .finish()
    }

    pub fn idx<D: IntoAtom, R: IntoAtom>(&self, degree: D, rep: R) -> Atom {
        FunctionBuilder::new(self.idx)
            .add_arg(degree.into_atom())
            .add_arg(rep.into_atom())
            .finish()
    }

    pub fn symmetric_generator_trace<R: IntoAtom, A: IntoAtom>(
        &self,
        rep: R,
        adjoint_indices: impl IntoIterator<Item = A>,
    ) -> Atom {
        let factors = adjoint_indices
            .into_iter()
            .map(|adjoint| self.chain_t(adjoint.into_atom()))
            .collect::<Vec<_>>();
        spenso::shadowing::trace_sym(rep.into_atom(), factors)
    }

    #[cfg(test)]
    pub(crate) fn initialize_tensor_symbols(&self) {
        let _ = self.t;
        let _ = self.f;
        let _ = self.d;

        let _ = self.gram;
        let _ = self.cas;
        let _ = self.idx;
    }

    // Generator for the adjoint representation of SU(N)
    pub fn t_strct<Aind: AbsInd>(
        &self,
        fundimd: impl Into<Dimension>,
        adim: impl Into<Dimension>,
    ) -> Canonicalized<ExplicitKey<Aind>> {
        let nc = fundimd.into();
        ExplicitKey::from_iter(
            [
                ColorAdjoint {}.new_rep(adim).cast(),
                ColorFundamental {}.new_rep(nc).to_lib(),
                ColorAntiFundamental {}.new_rep(nc).cast(),
            ],
            self.t,
            None,
        )
    }
    pub fn t_pattern(
        &self,
        fundimd: impl Into<Dimension>,
        adim: impl Into<Dimension>,
        a: impl Into<AbstractIndex>,
        i: impl Into<AbstractIndex>,
        j: impl Into<AbstractIndex>,
    ) -> Atom {
        let structure = self.t_strct::<AbstractIndex>(fundimd, adim);
        let logical_indices = [a.into(), i.into(), j.into()];
        let storage_indices = structure.layout().logical_to_canonical(&logical_indices);
        structure
            .into_canonical()
            .reindex_storage(&storage_indices)
            .unwrap()
            .permute_with_metric()
    }

    pub fn f_strct<Aind: AbsInd>(
        &self,
        adim: impl Into<Dimension>,
    ) -> Canonicalized<ExplicitKey<Aind>> {
        let adim = adim.into();
        ExplicitKey::from_iter(
            [
                ColorAdjoint {}.new_rep(adim),
                ColorAdjoint {}.new_rep(adim),
                ColorAdjoint {}.new_rep(adim),
            ],
            self.f,
            None,
        )
    }

    pub fn f_pattern(
        &self,
        adim: impl Into<Dimension>,
        a: impl Into<AbstractIndex>,
        b: impl Into<AbstractIndex>,
        c: impl Into<AbstractIndex>,
    ) -> Atom {
        let structure = self.f_strct::<AbstractIndex>(adim);
        let logical_indices = [a.into(), b.into(), c.into()];
        let storage_indices = structure.layout().logical_to_canonical(&logical_indices);
        structure
            .into_canonical()
            .reindex_storage(&storage_indices)
            .unwrap()
            .permute_with_metric()
    }
}

#[derive(Clone, Copy)]
enum ColorInvariantPrintKind {
    Gram,
    Casimir,
    Index,
}

impl ColorInvariantPrintKind {
    fn print(self, atom: AtomView<'_>, opt: &PrintOptions, state: &PrintState) -> Option<String> {
        if opt.mode.is_mathml() {
            return self.mathml(atom, opt, state);
        }
        let resolved = SpensoPrintSettings::resolve(opt)?;
        let settings = resolved.presentation;
        let AtomView::Fun(f) = atom else { return None };
        let args = f.iter().collect::<Vec<_>>();
        let arity = if matches!(self, Self::Gram) { 3 } else { 2 };
        if args.len() != arity {
            return None;
        }
        if !settings.explicit_invariants && !settings.with_dim && small_integer(args[0]) == Some(2)
        {
            let alias = match self {
                Self::Casimir if is_color_rep(args[1], "spenso::cof") => Some(("C_F", "CF")),
                Self::Casimir if is_color_rep(args[1], "spenso::coad") => Some(("C_A", "CA")),
                Self::Index if is_color_rep(args[1], "spenso::cof") => Some(("T_R", "TR")),
                _ => None,
            };
            if let Some((script, plain)) = alias {
                let label = if settings.symbol_scripts {
                    script.to_owned()
                } else {
                    Self::function_head(plain, opt)
                };
                return Some(Self::colorize(&label, opt));
            }
        }
        let (symbol, function) = match self {
            Self::Gram => ("G", "gram"),
            Self::Casimir => ("C", "cas"),
            Self::Index => ("I", "idx"),
        };
        let mut degree = String::new();
        args[0].format(&mut degree, opt, PrintState::new()).unwrap();
        let head = if settings.symbol_scripts {
            let (open, close) = resolved.script_delimiters();
            // Group symbolic or multi-digit degrees so following powers apply to
            // the complete invariant, never to part of its degree label.
            if degree.chars().count() == 1 {
                format!("{symbol}_{degree}")
            } else {
                format!("{symbol}_{open}{degree}{close}")
            }
        } else {
            Self::function_head(function, opt)
        };
        let mut arguments = args[1..]
            .iter()
            .map(|arg| Self::representation(*arg, opt))
            .collect::<Vec<_>>();
        if !settings.symbol_scripts {
            arguments.insert(0, degree);
        }
        // These are scalar function arguments, not index lists: always delimit
        // them, independently of the tensor parentheses/comma preferences.
        Some(format!(
            "{}({})",
            Self::colorize(&head, opt),
            arguments.join(",")
        ))
    }

    fn mathml_settings(opt: &PrintOptions) -> SpensoPrintSettings {
        match opt.custom_print_mode.get("spenso") {
            Some(PrintUserData::Integer(encoded)) => SpensoPrintSettings::from(*encoded as usize),
            _ => SpensoPrintSettings::typst(),
        }
    }

    fn mathml_fits(atom: AtomView<'_>, depth: usize, remaining: &mut usize) -> bool {
        if depth > 12
            || *remaining == 0
            || atom.get_byte_size() > 2048
            || atom
                .get_symbol()
                .is_some_and(|symbol| symbol.get_name().len() > 1024)
        {
            return false;
        }
        *remaining -= 1;
        atom.children()
            .all(|child| Self::mathml_fits(child, depth + 1, remaining))
    }

    fn mathml_source(atom: AtomView<'_>, opt: &PrintOptions) -> Option<String> {
        if !Self::mathml_fits(atom, 0, &mut 64) {
            return None;
        }
        if let AtomView::Var(variable) = atom
            && let Some(display) = IndexDisplay::from_symbol(variable.get_symbol())
        {
            return display.to_mathml();
        }
        let source = atom
            .printer(PrintOptions {
                max_terms: Some(64),
                ..opt.clone()
            })
            .to_string();
        mathml::fragment(source)
    }

    fn mathml_group(source: &str) -> String {
        format!("<mrow><mo>(</mo>{source}<mo>)</mo></mrow>")
    }

    fn mathml_representation(
        mut atom: AtomView<'_>,
        opt: &PrintOptions,
        settings: SpensoPrintSettings,
    ) -> Option<String> {
        let original = atom;
        let mut dual = false;
        if let AtomView::Fun(wrapper) = atom
            && wrapper.get_symbol() == AIND_SYMBOLS.dind
            && wrapper.get_nargs() == 1
        {
            dual = true;
            atom = wrapper.iter().next()?;
        }
        let AtomView::Fun(rep) = atom else {
            return Self::mathml_source(original, opt);
        };
        if !rep.get_symbol().has_tag(&T.representation) || rep.get_nargs() != 1 {
            return Self::mathml_source(original, opt);
        }
        let mut label = match rep.get_symbol().get_name() {
            "spenso::cof" => IndexDisplay::Symbol("F".to_owned()).to_mathml()?,
            "spenso::coad" => IndexDisplay::Symbol("A".to_owned()).to_mathml()?,
            _ => RepresentationMetadata::from_symbol(rep.get_symbol())?
                .label
                .to_mathml()?,
        };
        if dual {
            label = format!("<mover accent=\"true\">{label}<mo>¯</mo></mover>");
        }
        if settings.with_dim {
            let dimension = Self::mathml_source(rep.iter().next()?, opt)?;
            label = if settings.symbol_scripts {
                format!("<msub>{label}{dimension}</msub>")
            } else {
                format!("<mrow>{label}<mo>[</mo>{dimension}<mo>]</mo></mrow>")
            };
        }
        (label.len() <= 4096).then_some(label)
    }

    fn mathml(self, atom: AtomView<'_>, opt: &PrintOptions, state: &PrintState) -> Option<String> {
        if !Self::mathml_fits(atom, 0, &mut 64) {
            return None;
        }
        let AtomView::Fun(function) = atom else {
            return None;
        };
        let arity = if matches!(self, Self::Gram) { 3 } else { 2 };
        if function.get_nargs() != arity {
            return None;
        }
        let args = function.iter().collect::<Vec<_>>();
        let settings = Self::mathml_settings(opt);
        let alias = if !settings.explicit_invariants
            && !settings.with_dim
            && small_integer(args[0]) == Some(2)
        {
            match self {
                Self::Casimir if is_color_rep(args[1], "spenso::cof") => Some(("C", "F", "CF")),
                Self::Casimir if is_color_rep(args[1], "spenso::coad") => Some(("C", "A", "CA")),
                Self::Index if is_color_rep(args[1], "spenso::cof") => Some(("T", "R", "TR")),
                _ => None,
            }
        } else {
            None
        };
        let mut source = if let Some((head, script, plain)) = alias {
            if settings.symbol_scripts {
                format!("<msub><mi>{head}</mi><mi>{script}</mi></msub>")
            } else {
                format!("<mi mathvariant=\"normal\">{plain}</mi>")
            }
        } else {
            let degree = Self::mathml_source(args[0], opt)?;
            let (symbol, name) = match self {
                Self::Gram => ("G", "gram"),
                Self::Casimir => ("C", "cas"),
                Self::Index => ("I", "idx"),
            };
            let head = if settings.symbol_scripts {
                format!("<msub><mi>{symbol}</mi>{degree}</msub>")
            } else {
                format!("<mi mathvariant=\"normal\">{name}</mi>")
            };
            let mut arguments = args[1..]
                .iter()
                .map(|arg| Self::mathml_representation(*arg, opt, settings))
                .collect::<Option<Vec<_>>>()?;
            if !settings.symbol_scripts {
                arguments.insert(0, degree);
            }
            format!(
                "<mrow>{head}<mo>\u{2061}</mo>{}</mrow>",
                Self::mathml_group(&arguments.join("<mo>,</mo>"))
            )
        };
        if state.in_exp_base && alias.is_none() {
            source = Self::mathml_group(&source);
        }
        (source.len() <= 4096).then_some(source)
    }

    fn representation(mut arg: AtomView<'_>, opt: &PrintOptions) -> String {
        let resolved = SpensoPrintSettings::resolve(opt).unwrap();
        let original = arg;
        let mut dual = false;
        if let AtomView::Fun(f) = arg
            && f.get_symbol() == AIND_SYMBOLS.dind
            && f.get_nargs() == 1
        {
            dual = true;
            arg = f.iter().next().unwrap();
        }
        let mut out = String::new();
        if let AtomView::Fun(rep) = arg
            && rep.get_symbol().has_tag(&T.representation)
            && rep.get_nargs() == 1
        {
            out = match rep.get_symbol().get_name() {
                "spenso::cof" => "F".to_owned(),
                "spenso::coad" => "A".to_owned(),
                _ => {
                    Atom::var(rep.get_symbol())
                        .format(&mut out, opt, PrintState::new())
                        .unwrap();
                    out
                }
            };
            if dual {
                out = match resolved.backend {
                    SpensoPrintBackend::Typst => format!("accent({out},macron)"),
                    SpensoPrintBackend::Latex => format!("\\overline{{{out}}}"),
                    SpensoPrintBackend::Plain => format!("bar({out})"),
                };
            }
            if resolved.presentation.with_dim {
                let mut dimension = String::new();
                rep.iter()
                    .next()
                    .unwrap()
                    .format(&mut dimension, opt, PrintState::new())
                    .unwrap();
                if resolved.presentation.symbol_scripts {
                    let (open, close) = resolved.script_delimiters();
                    out.push_str(&format!("_{open}{dimension}{close}"));
                } else {
                    out.push_str(&format!("[{dimension}]"));
                }
            }
        } else {
            original.format(&mut out, opt, PrintState::new()).unwrap();
        }
        out
    }

    fn colorize(head: &str, opt: &PrintOptions) -> String {
        if opt.color_builtin_symbols {
            nu_ansi_term::Color::Magenta.paint(head).to_string()
        } else {
            head.to_owned()
        }
    }

    fn function_head(name: &str, opt: &PrintOptions) -> String {
        match SpensoPrintSettings::resolve(opt).unwrap().backend {
            SpensoPrintBackend::Typst => format!("op({name:?})"),
            SpensoPrintBackend::Latex => format!("\\operatorname{{{name}}}"),
            SpensoPrintBackend::Plain => name.to_owned(),
        }
    }

    fn color_count(atom: AtomView<'_>, opt: &PrintOptions, state: &PrintState) -> Option<String> {
        if !matches!(atom, AtomView::Var(_)) {
            return None;
        }
        if opt.mode.is_mathml() {
            let source = if Self::mathml_settings(opt).symbol_scripts {
                "<msub><mi>N</mi><mi>c</mi></msub>"
            } else {
                "<mi mathvariant=\"normal\">Nc</mi>"
            };
            return Some(if state.in_exp_base {
                Self::mathml_group(source)
            } else {
                source.to_owned()
            });
        }
        let settings = SpensoPrintSettings::resolve(opt)?.presentation;
        let out = if settings.symbol_scripts {
            "N_c".to_owned()
        } else {
            Self::function_head("Nc", opt)
        };
        Some(Self::colorize(&out, opt))
    }
}

fn small_integer(expr: AtomView<'_>) -> Option<i64> {
    let AtomView::Num(number) = expr else {
        return None;
    };
    let CoefficientView::Natural(value, 1, 0, 1) = number.get_coeff_view() else {
        return None;
    };
    Some(value)
}

fn is_color_rep(expr: AtomView<'_>, name: &str) -> bool {
    matches!(
        expr,
        AtomView::Fun(rep) if rep.get_symbol().get_name() == name && rep.get_nargs() == 1
    )
}

spenso::symbolica_init_lazy_static! {
pub static CS, CS_INNER: ColorSymbols = || {
    fn representation_symbol(rep: Atom) -> Symbol {
        let AtomView::Fun(f) = rep.as_view() else {
            unreachable!("Color representations are symbolic functions")
        };
        f.get_symbol()
    }

    ColorSymbols {
        t: tensor_symbol!("spenso::t"; Real; print = |a, opt, state| {
            if opt.mode.is_mathml() {
                return tags::tensor_print(a, opt, state);
            }

            match SpensoPrintSettings::resolve(opt) {
                Some(resolved)=>{
                    let (script_open, script_close) = resolved.script_delimiters();
                    let SpensoPrintSettings{
                        parens,
                        symbol_scripts,
                        commas,..
                    } = resolved.presentation;


                    let AtomView::Fun(f)=a else {
                        return None;
                    };
                    if f.get_nargs()!=3 {
                        return None;
                    }
                    let mut argitem = f.iter();
                    let a = argitem.next().unwrap();
                    let b = argitem.next().unwrap();
                    let mut c = argitem.next().unwrap();

                    let mut out = "t".to_string();
                    if symbol_scripts {
                        out.push('^');
                    }
                    if opt.color_builtin_symbols {
                        out = nu_ansi_term::Color::Magenta.paint(out).to_string();
                    }

                    if parens{
                        out.push(script_open);
                    }
                    a.format(&mut out, opt, PrintState::new()).unwrap();
                    if commas{
                        out.push(',');
                    } else {
                        out.push(' ');
                    }
                    b.format(&mut out, opt, PrintState::new()).unwrap();
                    if parens{
                        out.push(script_close);
                    }
                    if symbol_scripts{
                        out.push('_');
                    }


                    if parens{
                        out.push(script_open);
                    }else if !symbol_scripts{
                        out.push(' ');
                    }

                    let AtomView::Fun(f)=c else {
                        return None;
                    };
                    if f.get_nargs()!=1 {
                        return None;
                    }
                    if f.get_symbol()!=AIND_SYMBOLS.dind{
                        return None;
                    }
                    c = f.iter().next().unwrap();
                    c.format(&mut out, opt, PrintState::new()).unwrap();
                    if parens{
                        out.push(script_close);
                    }
                    Some(out)
                }
                _=>None}

        }),
        f: tensor_symbol!("spenso::f"; Real, Antisymmetric; print = |a, opt, state| {
            if opt.mode.is_mathml() {
                return tags::tensor_print(a, opt, state);
            }

            match SpensoPrintSettings::resolve(opt) {
                Some(resolved)=>{
                    let (script_open, script_close) = resolved.script_delimiters();
                    let SpensoPrintSettings{
                        parens,
                        symbol_scripts,
                        commas,..
                    } = resolved.presentation;


                    let AtomView::Fun(f)=a else {
                        return None;
                    };
                    if f.get_nargs()!=3 {
                        return None;
                    }
                    let mut argitem = f.iter();
                    let a = argitem.next().unwrap();
                    let b = argitem.next().unwrap();
                    let c = argitem.next().unwrap();

                    let mut out = "f".to_string();
                    if symbol_scripts {
                        out.push('^');
                    }
                    if opt.color_builtin_symbols {
                        out = nu_ansi_term::Color::Magenta.paint(out).to_string();
                    }

                    if parens{
                        out.push(script_open);
                    }
                    a.format(&mut out, opt, PrintState::new()).unwrap();
                    if commas{
                        out.push(',');
                    } else {
                        out.push(' ');
                    }
                    b.format(&mut out, opt, PrintState::new()).unwrap();
                    if commas{
                        out.push(',');
                    } else {
                        out.push(' ');
                    }
                    c.format(&mut out, opt, PrintState::new()).unwrap();
                    if parens{
                        out.push(script_close);
                    }
                    Some(out)
                }
                _=>None}

        }),
        d: tensor_symbol!("spenso::d"),
        gram: symbol!("spenso::gram"; Real, Scalar; print = |a, opt, state| {
            ColorInvariantPrintKind::Gram.print(a, opt, state)
        }),
        cas: symbol!("spenso::cas"; Real, Scalar; print = |a, opt, state| {
            ColorInvariantPrintKind::Casimir.print(a, opt, state)
        }),
        idx: symbol!("spenso::idx"; Real, Scalar; print = |a, opt, state| {
            ColorInvariantPrintKind::Index.print(a, opt, state)
        }),
        fundamental_rep: representation_symbol(
            ColorFundamental {}.to_symbolic(std::iter::empty::<Atom>()),
        ),
        adjoint_rep: representation_symbol(ColorAdjoint {}.to_symbolic(std::iter::empty::<Atom>())),
        adj_: symbol!("adj_"),
        nc_: symbol!("nc_"),
        nc: symbol!("spenso::Nc";Real; print = ColorInvariantPrintKind::color_count,eval = EvaluationInfo::constant(|_tags, prec| Ok(Rational::new(3,1).to_multi_prec_float(prec).into()))),
    }
};
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ColorSimplifySettings {
    /// Whether closed color chains should be evaluated as traces.
    pub evaluate_traces: bool,
    /// Whether contractions between generators on different open chains or traces should
    /// be expanded with the fundamental Fierz identity.
    pub expand_cross_chain_fierz: bool,
    /// Whether invariant factors for `cof(N)` should be written directly in
    /// terms of the fundamental dimension.
    pub substitute_cof_dimension_invariants: bool,
    /// Whether a trace line that nothing else in its term can reach is
    /// decomposed completely in one kernel call. The kernel then also
    /// certifies colour fixed points, so the planner needs no confirming
    /// colour round, and distributes a colour sum with ports only when a
    /// rule can span it. Without it, every insertion returns to the planner.
    pub one_shot_traces: bool,
}

impl Default for ColorSimplifySettings {
    fn default() -> Self {
        Self {
            evaluate_traces: true,
            expand_cross_chain_fierz: true,
            substitute_cof_dimension_invariants: false,
            one_shot_traces: true,
        }
    }
}

impl ColorSimplifySettings {
    /// Leaves collected `trace(...)` nodes unevaluated; Fierz contractions can still join them.
    pub fn without_trace_evaluation(mut self) -> Self {
        self.evaluate_traces = false;
        self
    }

    /// Keeps separate open chains and traces instead of applying cross-chain Fierz
    /// expansion.
    pub fn without_cross_chain_fierz_expansion(mut self) -> Self {
        self.expand_cross_chain_fierz = false;
        self
    }

    /// Rewrites supported `cof(N)` invariants to explicit dimension formulas.
    pub fn with_cof_dimension_invariants(mut self) -> Self {
        self.substitute_cof_dimension_invariants = true;
        self
    }

    /// Returns to the planner after each trace insertion and confirms every
    /// colour fixed point with a no-op round.
    pub fn without_one_shot_traces(mut self) -> Self {
        self.one_shot_traces = false;
        self
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ColorCasimirSettings {
    /// Rewrite the explicit fundamental dimension using `d_F = C_A`.
    pub rewrite_fundamental_dimension: bool,
    /// Substitute the common fundamental normalization `T_F = 1/2`.
    pub substitute_fundamental_index: bool,
}

impl Default for ColorCasimirSettings {
    fn default() -> Self {
        Self {
            rewrite_fundamental_dimension: true,
            substitute_fundamental_index: false,
        }
    }
}

impl ColorCasimirSettings {
    /// Keep `d_F` explicit instead of applying the SU(N) relation `d_F = C_A`.
    pub fn without_fundamental_dimension_rewrite(mut self) -> Self {
        self.rewrite_fundamental_dimension = false;
        self
    }

    /// Apply the standard fundamental index normalization `T_F = 1/2`.
    pub fn with_fundamental_index_normalization(mut self) -> Self {
        self.substitute_fundamental_index = true;
        self
    }
}

/// Rewrite color dimensions and scalar invariants between SU(N) conventions.
/// Indexed generator, structure-constant and trace reduction belongs to
/// `simplify_algebra` method on [`SymbolicTensor`](crate::tensor::SymbolicTensor).
pub trait ColorSimplifier {
    /// Replace concrete QCD fundamental and adjoint dimensions by their
    /// parametric SU(Nc) expressions while preserving all color indices.
    ///
    /// This is useful for exact color-algebra comparisons: public expressions
    /// may use the physical `cof(3)` and `coad(8)` representations, while an
    /// intermediate symbolic calculation should retain its `Nc` dependence.
    fn to_parametric_color(&self) -> Atom;

    /// Rewrites the explicit representation dimensions into a Casimir basis.
    fn to_color_casimir(&self, fundamental_rep: AtomView<'_>, adjoint_rep: AtomView<'_>) -> Atom;

    /// Rewrites explicit dimensions with control over SU(N) normalization choices.
    fn to_color_casimir_with(
        &self,
        fundamental_rep: AtomView<'_>,
        adjoint_rep: AtomView<'_>,
        settings: ColorCasimirSettings,
    ) -> Atom;

    /// Rewrites supported `cof(N)` invariant factors into explicit dimension formulas.
    fn to_cof_dimension_invariants(&self) -> Atom;
}
impl ColorSimplifier for Atom {
    fn to_parametric_color(&self) -> Atom {
        self.as_view().to_parametric_color()
    }

    fn to_color_casimir(&self, fundamental_rep: AtomView<'_>, adjoint_rep: AtomView<'_>) -> Atom {
        self.to_color_casimir_with(
            fundamental_rep,
            adjoint_rep,
            ColorCasimirSettings::default(),
        )
    }

    fn to_color_casimir_with(
        &self,
        fundamental_rep: AtomView<'_>,
        adjoint_rep: AtomView<'_>,
        settings: ColorCasimirSettings,
    ) -> Atom {
        casimir::color_casimir_basis_impl(
            self.as_atom_view(),
            fundamental_rep,
            adjoint_rep,
            settings,
        )
    }

    fn to_cof_dimension_invariants(&self) -> Atom {
        CofDimensionInvariantRewriter.run(self.as_atom_view())
    }
}

impl ColorSimplifier for AtomView<'_> {
    fn to_parametric_color(&self) -> Atom {
        let adjoint = ColorAdjoint {};
        let fundamental = ColorFundamental {};
        let nc = Atom::var(CS.nc);
        self.replace(adjoint.to_symbolic([RS.d_, RS.a_]))
            .with(
                adjoint.to_symbolic([nc.clone().pow(Atom::num(2)) - Atom::one(), Atom::var(RS.a_)]),
            )
            .replace(fundamental.to_symbolic([RS.d_, RS.a_]))
            .with(fundamental.to_symbolic([nc, Atom::var(RS.a_)]))
    }

    fn to_color_casimir(&self, fundamental_rep: AtomView<'_>, adjoint_rep: AtomView<'_>) -> Atom {
        self.to_color_casimir_with(
            fundamental_rep,
            adjoint_rep,
            ColorCasimirSettings::default(),
        )
    }

    fn to_color_casimir_with(
        &self,
        fundamental_rep: AtomView<'_>,
        adjoint_rep: AtomView<'_>,
        settings: ColorCasimirSettings,
    ) -> Atom {
        casimir::color_casimir_basis_impl(
            self.as_atom_view(),
            fundamental_rep,
            adjoint_rep,
            settings,
        )
    }

    fn to_cof_dimension_invariants(&self) -> Atom {
        CofDimensionInvariantRewriter.run(self.as_atom_view())
    }
}

#[cfg(test)]
mod test;

#[cfg(test)]
mod mathml_tests {
    use super::*;
    use crate::test_support;
    use spenso::structure::representation::{IndexPalette, IndexRow, RepresentationClass};
    use symbolica::atom::{NamespacedSymbol, SymbolBuilder};
    use symbolica::function;

    fn options(settings: SpensoPrintSettings) -> PrintOptions {
        PrintOptions {
            custom_print_mode: settings.into(),
            ..PrintOptions::mathml()
        }
    }

    #[test]
    fn mathml_quadratic_aliases_and_color_count_honor_symbol_scripts() {
        test_support::test_initialize();
        let fundamental = ColorFundamental {}.to_symbolic([Atom::num(3)]);
        let adjoint = ColorAdjoint {}.to_symbolic([Atom::num(8)]);
        for (kind, atom, script, plain) in [
            (
                ColorInvariantPrintKind::Casimir,
                CS.cas(2, fundamental.clone()),
                "<msub><mi>C</mi><mi>F</mi></msub>",
                "CF",
            ),
            (
                ColorInvariantPrintKind::Casimir,
                CS.cas(2, adjoint),
                "<msub><mi>C</mi><mi>A</mi></msub>",
                "CA",
            ),
            (
                ColorInvariantPrintKind::Index,
                CS.idx(2, fundamental),
                "<msub><mi>T</mi><mi>R</mi></msub>",
                "TR",
            ),
        ] {
            assert_eq!(
                kind.print(atom.as_view(), &PrintOptions::mathml(), &PrintState::new())
                    .unwrap(),
                script
            );
            assert_eq!(
                kind.print(
                    atom.as_view(),
                    &options(SpensoPrintSettings::compact()),
                    &PrintState::new()
                )
                .unwrap(),
                format!("<mi mathvariant=\"normal\">{plain}</mi>")
            );
            let native = atom.printer(PrintOptions::mathml()).to_string();
            assert!(native.contains(script), "{native}");
            assert!(!native.contains("spenso::"));
            for exponent in [4, -4, -1] {
                let power = atom.pow(Atom::num(exponent));
                for (opt, alias) in [
                    (PrintOptions::mathml(), script),
                    (options(SpensoPrintSettings::compact()), plain),
                ] {
                    let native = power.printer(opt).to_string();
                    assert!(native.contains(alias), "{native}");
                    assert!(!native.contains("<mo>(</mo>"), "{native}");
                    assert!(
                        native.contains(if exponent == -1 { "<mfrac>" } else { "<msup>" }),
                        "{native}"
                    );
                }
            }
        }
        assert_eq!(
            ColorInvariantPrintKind::color_count(
                Atom::var(CS.nc).as_view(),
                &PrintOptions::mathml(),
                &PrintState::new()
            )
            .unwrap(),
            "<msub><mi>N</mi><mi>c</mi></msub>"
        );
        assert_eq!(
            ColorInvariantPrintKind::color_count(
                Atom::var(CS.nc).as_view(),
                &options(SpensoPrintSettings::compact()),
                &PrintState::new()
            )
            .unwrap(),
            "<mi mathvariant=\"normal\">Nc</mi>"
        );
    }

    #[test]
    fn mathml_explicit_invariants_keep_degree_dimensions_and_duality() {
        test_support::test_initialize();
        let fundamental = ColorFundamental {}.to_symbolic([Atom::num(3)]);
        let adjoint = ColorAdjoint {}.to_symbolic([Atom::num(8)]);
        let opt = options(SpensoPrintSettings {
            explicit_invariants: true,
            with_dim: true,
            ..SpensoPrintSettings::typst()
        });
        for (kind, atom, head, degree) in [
            (
                ColorInvariantPrintKind::Casimir,
                CS.cas(12, fundamental.clone()),
                "C",
                "12",
            ),
            (
                ColorInvariantPrintKind::Index,
                CS.idx(2, fundamental.clone()),
                "I",
                "2",
            ),
            (
                ColorInvariantPrintKind::Gram,
                CS.gram(4, fundamental.clone(), adjoint),
                "G",
                "4",
            ),
        ] {
            let source = kind
                .print(atom.as_view(), &opt, &PrintState::new())
                .unwrap();
            assert!(
                source.contains(&format!("<msub><mi>{head}</mi>")),
                "{source}"
            );
            assert!(source.contains(&format!("<mn>{degree}</mn>")), "{source}");
            assert!(source.contains("<msub><mi>F</mi>"), "{source}");
            assert!(source.contains("<mn>3</mn>"), "{source}");
            assert!(!source.contains("data-symbolica-path"));
            if matches!(kind, ColorInvariantPrintKind::Gram) {
                assert!(source.contains("<mo>,</mo><msub><mi>A</mi>"), "{source}");
                assert!(source.contains("<mn>8</mn>"));
            }
        }
        let dual = function!(AIND_SYMBOLS.dind, fundamental.clone());
        let atom = CS.cas(2, dual);
        let source = ColorInvariantPrintKind::Casimir
            .print(atom.as_view(), &opt, &PrintState::new())
            .unwrap();
        assert!(source.contains("<msub><mover accent=\"true\"><mi>F</mi><mo>¯</mo></mover>"));
        let unscripted = options(SpensoPrintSettings {
            explicit_invariants: true,
            with_dim: true,
            ..SpensoPrintSettings::compact()
        });
        let source = ColorInvariantPrintKind::Casimir
            .print(
                CS.cas(2, fundamental).as_view(),
                &unscripted,
                &PrintState::new(),
            )
            .unwrap();
        assert!(source.contains("<mi mathvariant=\"normal\">cas</mi>"));
        assert!(source.contains("<mo>,</mo><mrow><mi>F</mi><mo>[</mo>"));
        assert!(source.contains("<mn>3</mn>"));
    }

    #[test]
    fn mathml_representation_metadata_and_exponent_bases_are_grouped() {
        test_support::test_initialize();
        let label = IndexDisplay::math(
            "frac",
            vec![
                IndexDisplay::Symbol("alpha".into()),
                IndexDisplay::Number(2),
            ],
        )
        .unwrap();
        let symbol = SymbolBuilder::new(NamespacedSymbol::parse(
            "idenso_color_mathml_tests::fractional_representation",
        ))
        .with_tags(vec![T.representation.clone(), T.dualizable.clone()])
        .with_user_data(
            RepresentationMetadata {
                class: RepresentationClass::Dualizable,
                label: label.clone(),
                index_palette: IndexPalette::Numeric,
                index_row: IndexRow::Top,
            }
            .to_user_data(),
        )
        .build()
        .unwrap();
        let atom = CS.cas(4, function!(symbol, 7));
        let source = ColorInvariantPrintKind::Casimir
            .print(atom.as_view(), &PrintOptions::mathml(), &PrintState::new())
            .unwrap();
        assert!(source.contains(&label.to_mathml().unwrap()));
        let state = PrintState {
            in_exp_base: true,
            ..PrintState::new()
        };
        assert_eq!(
            ColorInvariantPrintKind::Casimir
                .print(atom.as_view(), &PrintOptions::mathml(), &state)
                .unwrap(),
            ColorInvariantPrintKind::mathml_group(&source)
        );
        let native = atom
            .pow(Atom::num(2))
            .printer(PrintOptions::mathml())
            .to_string();
        assert!(native.contains("<msup>"), "{native}");
        assert!(native.contains("<mo>(</mo>"));
        assert!(native.contains("<mfrac><mi>α</mi><mn>2</mn></mfrac>"));
    }

    #[test]
    fn mathml_dimension_callbacks_preserve_literal_path_text() {
        test_support::test_initialize();
        let dimension = SymbolBuilder::new(NamespacedSymbol::parse(
            "idenso_color_mathml_tests::literal_path_dimension",
        ))
        .with_print_function(|_, options, _| {
            options.mode.is_mathml().then(|| {
                r#"<mrow data-symbolica-path="dimension-owned"><mtext>literal data-symbolica-path="keep"</mtext></mrow>"#.to_owned()
            })
        })
        .build()
        .unwrap();
        let representation = ColorFundamental {}.to_symbolic([Atom::var(dimension)]);
        let atom = CS.cas(4, representation);
        let original = atom.as_view().get_data().to_vec();
        let opt = options(SpensoPrintSettings {
            with_dim: true,
            ..SpensoPrintSettings::typst()
        });
        let source = ColorInvariantPrintKind::Casimir
            .print(atom.as_view(), &opt, &PrintState::new())
            .unwrap();
        assert_eq!(source.matches("data-symbolica-path").count(), 1, "{source}");
        for source in [source, atom.printer(opt).to_string()] {
            assert!(
                source.contains(r#"<mtext>literal data-symbolica-path="keep"</mtext>"#),
                "{source}"
            );
            assert!(
                !source.contains(r#"data-symbolica-path="dimension-owned""#),
                "{source}"
            );
        }
        assert_eq!(atom.as_view().get_data(), original);
    }

    #[test]
    fn mathml_invariant_callbacks_reject_large_or_malformed_inputs() {
        test_support::test_initialize();
        let representation = ColorFundamental {}.to_symbolic([Atom::num(3)]);
        let degree = (0..70)
            .fold(
                FunctionBuilder::new(symbol!("idenso_color_mathml_tests::large_degree")),
                |function, value| function.add_arg(Atom::num(value)),
            )
            .finish();
        let atom = CS.cas(degree, representation);
        assert!(
            ColorInvariantPrintKind::Casimir
                .print(atom.as_view(), &PrintOptions::mathml(), &PrintState::new())
                .is_none()
        );
        let malformed = function!(CS.cas, 2);
        assert!(
            ColorInvariantPrintKind::Casimir
                .print(
                    malformed.as_view(),
                    &PrintOptions::mathml(),
                    &PrintState::new()
                )
                .is_none()
        );
    }
}
