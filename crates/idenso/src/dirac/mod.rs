use std::sync::LazyLock;

use spenso::{
    network::library::symbolic::{ETS, ExplicitKey},
    shadowing::symbolica_utils::SpensoPrintSettings,
    structure::{
        dimension::Dimension,
        representation::{LibraryRep, Minkowski, RepName},
        slot::{AbsInd, DummyAind, ParseableAind},
    },
    utils::to_superscript,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    function,
    id::{Context, Replacement},
    printer::PrintState,
    symbol,
    utils::Settable,
};

use crate::{IndexTooling, color::SelectiveExpand, metric::MetricSimplifier, rep_symbols::RS};
use eyre::Result;

use super::representations::Bispinor;

pub struct GammaLibrary {
    pub gamma: Symbol,
    pub gammaconj: Symbol,
    pub gammaadj: Symbol,
    pub gamma0: Symbol,
    pub projp: Symbol,
    pub projm: Symbol,
    pub gamma5: Symbol,
    pub sigma: Symbol,
}

impl GammaLibrary {
    pub fn replace_with(&self, rep: &Self) -> Vec<Replacement> {
        vec![
            Replacement::new(
                function!(self.gamma, RS.i__).to_pattern(),
                function!(rep.gamma, RS.i__),
            ),
            Replacement::new(
                function!(self.projp, RS.i__).to_pattern(),
                function!(rep.projp, RS.i__),
            ),
            Replacement::new(
                function!(self.projm, RS.i__).to_pattern(),
                function!(rep.projm, RS.i__),
            ),
            Replacement::new(
                function!(self.gamma5, RS.i__).to_pattern(),
                function!(rep.gamma5, RS.i__),
            ),
            Replacement::new(
                function!(self.sigma, RS.i__).to_pattern(),
                function!(rep.sigma, RS.i__),
            ),
        ]
    }
}

pub struct GammaSymbolsInternal {
    pub gamma_chain: Symbol,
    pub gamma_trace: Symbol,
}

impl GammaSymbolsInternal {
    pub fn chain_pattern<'a, It: Into<AtomOrView<'a>>>(
        &self,
        i: impl Into<AtomOrView<'a>>,
        j: impl Into<AtomOrView<'a>>,
        mus: impl IntoIterator<Item = It>,
    ) -> Atom {
        let mus: Vec<_> = mus
            .into_iter()
            .enumerate()
            .map(|(idx, a)| {
                if idx % 2 == 0 {
                    Minkowski {}.to_symbolic([a])
                } else {
                    Bispinor {}.to_symbolic([a])
                }
            })
            .collect();
        if mus.len() % 2 == 0 {
            panic!("Should intersperse dummies")
        }
        FunctionBuilder::new(self.gamma_chain)
            .add_args(&mus)
            .add_arg(Bispinor {}.to_symbolic([i]))
            .add_arg(Bispinor {}.to_symbolic([j]))
            .finish()
    }
}

pub struct PolSymbols {
    pub eps: Symbol,
    pub ebar: Symbol,
    pub u: Symbol,
    pub ubar: Symbol,
    pub v: Symbol,
    pub vbar: Symbol,
}

pub static PS: LazyLock<PolSymbols> = LazyLock::new(|| PolSymbols {
    eps: symbol!("spenso::ϵ"),
    ebar: symbol!("spenso::ϵbar"),
    u: symbol!("spenso::u"),
    ubar: symbol!("spenso::ubar"),
    v: symbol!("spenso::v"),
    vbar: symbol!("spenso::vbar"),
});

pub fn pol_conj_impl(expression: AtomView) -> Atom {
    let expr = expression.to_owned();

    expr.replace_multiple(&[
        Replacement::new(
            function!(PS.ebar, RS.i__).to_pattern(),
            function!(PS.eps, RS.i__),
        ),
        Replacement::new(
            function!(PS.eps, RS.i__).to_pattern(),
            function!(PS.ebar, RS.i__),
        ),
        Replacement::new(
            function!(PS.u, RS.i__).to_pattern(),
            function!(PS.ubar, RS.i__),
        ),
        Replacement::new(
            function!(PS.ubar, RS.i__).to_pattern(),
            function!(PS.u, RS.i__),
        ),
        Replacement::new(
            function!(PS.v, RS.i__).to_pattern(),
            function!(PS.vbar, RS.i__),
        ),
        Replacement::new(
            function!(PS.vbar, RS.i__).to_pattern(),
            function!(PS.v, RS.i__),
        ),
    ])
}

pub fn gamma_conj_impl(expression: AtomView) -> Atom {
    let expr = expression.to_owned();

    expr.replace(crate::gamma!(RS.a__, RS.i__, RS.j__).to_pattern())
        .with((-crate::gamma!(RS.a__, RS.j__, RS.i__)).to_pattern())
        .replace(crate::gamma5!(RS.i__, RS.j__).to_pattern())
        .with(crate::gamma5!(RS.j__, RS.i__).to_pattern())
}
pub static GS: LazyLock<GammaSymbolsInternal> = LazyLock::new(|| GammaSymbolsInternal {
    gamma_chain: symbol!("spenso::gamma_chain"),
    gamma_trace: symbol!("spenso::gamma_trace"),
});

pub static AGS: LazyLock<GammaLibrary> = LazyLock::new(|| GammaLibrary {
    gamma: symbol!(
        "spenso::gamma",
        print = |a, opt| {
            match opt.custom_print_mode {
                Some(("spenso", i)) => {
                    let SpensoPrintSettings {
                        parens,
                        symbol_scripts,
                        commas,
                        ..
                    } = SpensoPrintSettings::from(i);

                    let AtomView::Fun(f) = a else {
                        return None;
                    };
                    if f.get_nargs() != 3 {
                        return None;
                    }
                    let mut argitem = f.iter();
                    let i = argitem.next().unwrap();
                    let j = argitem.next().unwrap();
                    let mu = argitem.next().unwrap();

                    let mut out = "γ".to_string();
                    if symbol_scripts {
                        out.push('^');
                    }
                    if opt.color_builtin_symbols {
                        out = nu_ansi_term::Color::Magenta.paint(out).to_string();
                    }

                    if parens {
                        out.push('(');
                    }
                    mu.format(&mut out, opt, PrintState::new()).unwrap();

                    out.push(',');

                    i.format(&mut out, opt, PrintState::new()).unwrap();
                    if commas {
                        out.push(',');
                    } else {
                        out.push(' ');
                    }
                    j.format(&mut out, opt, PrintState::new()).unwrap();
                    if parens {
                        out.push(')');
                    }

                    Some(out)
                }
                _ => None,
            }
        }
    ),
    gammaadj: symbol!("spenso::gammaadj"),
    projp: symbol!(
        "spenso::projp",
        print = |a, opt| {
            match opt.custom_print_mode {
                Some(("spenso", i)) => {
                    let settings = SpensoPrintSettings::from(i);
                    let is_typst = settings.is_typst();
                    let SpensoPrintSettings {
                        parens,
                        symbol_scripts,
                        commas,
                        ..
                    } = settings;

                    let AtomView::Fun(f) = a else {
                        return None;
                    };
                    if f.get_nargs() != 2 {
                        return None;
                    }
                    let mut argitem = f.iter();
                    let i = argitem.next().unwrap();
                    let j = argitem.next().unwrap();

                    let mut out = if is_typst { "ℙ_p" } else { "ℙₚ" }.to_string();

                    if symbol_scripts {
                        out.push('^');
                    }
                    if opt.color_builtin_symbols {
                        out = nu_ansi_term::Color::Magenta.paint(out).to_string();
                    }

                    if parens {
                        out.push('(');
                    }
                    i.format(&mut out, opt, PrintState::new()).unwrap();
                    if commas {
                        out.push(',');
                    } else {
                        out.push(' ');
                    }
                    j.format(&mut out, opt, PrintState::new()).unwrap();
                    if parens {
                        out.push(')');
                    }

                    Some(out)
                }
                _ => None,
            }
        }
    ),
    projm: symbol!(
        "spenso::projm",
        print = |a, opt| {
            match opt.custom_print_mode {
                Some(("spenso", i)) => {
                    let settings = SpensoPrintSettings::from(i);
                    let is_typst = settings.is_typst();
                    let SpensoPrintSettings {
                        parens,
                        symbol_scripts,
                        commas,
                        ..
                    } = settings;

                    let AtomView::Fun(f) = a else {
                        return None;
                    };
                    if f.get_nargs() != 2 {
                        return None;
                    }
                    let mut argitem = f.iter();
                    let i = argitem.next().unwrap();
                    let j = argitem.next().unwrap();
                    let mut out = if is_typst { "ℙ_m" } else { "ℙₘ" }.to_string();

                    if symbol_scripts {
                        out.push('^');
                    }
                    if opt.color_builtin_symbols {
                        out = nu_ansi_term::Color::Magenta.paint(out).to_string();
                    }

                    if parens {
                        out.push('(');
                    }
                    i.format(&mut out, opt, PrintState::new()).unwrap();
                    if commas {
                        out.push(',');
                    } else {
                        out.push(' ');
                    }
                    j.format(&mut out, opt, PrintState::new()).unwrap();
                    if parens {
                        out.push(')');
                    }

                    Some(out)
                }
                _ => None,
            }
        }
    ),
    gamma5: symbol!(
        "spenso::gamma5",
        print = |a, opt| {
            match opt.custom_print_mode {
                Some(("spenso", i)) => {
                    let settings = SpensoPrintSettings::from(i);
                    let is_typst = settings.is_typst();
                    let SpensoPrintSettings {
                        parens,
                        symbol_scripts,
                        commas,
                        ..
                    } = settings;

                    let AtomView::Fun(f) = a else {
                        return None;
                    };
                    if f.get_nargs() != 2 {
                        return None;
                    }
                    let mut argitem = f.iter();
                    let i = argitem.next().unwrap();
                    let j = argitem.next().unwrap();

                    let mut out = "γ".to_string();
                    if is_typst {
                        out.push_str("_5");
                    } else {
                        out.push_str(&to_superscript(5));
                    }
                    if symbol_scripts {
                        out.push('^');
                    }
                    if opt.color_builtin_symbols {
                        out = nu_ansi_term::Color::Magenta.paint(out).to_string();
                    }

                    if parens {
                        out.push('(');
                    }
                    i.format(&mut out, opt, PrintState::new()).unwrap();
                    if commas {
                        out.push(',');
                    } else {
                        out.push(' ');
                    }
                    j.format(&mut out, opt, PrintState::new()).unwrap();
                    if parens {
                        out.push(')');
                    }

                    Some(out)
                }
                _ => None,
            }
        }
    ),
    sigma: symbol!("spenso::sigma"),
    gamma0: symbol!("spenso::gamma0";Real,Symmetric;print = |a, opt| {
        match opt.custom_print_mode {
            Some(("spenso", i)) => {
                let settings = SpensoPrintSettings::from(i);
                let is_typst = settings.is_typst();
                let SpensoPrintSettings {
                    parens,
                    symbol_scripts,
                    commas,
                    ..
                } = settings;

                let AtomView::Fun(f) = a else {
                    return None;
                };
                if f.get_nargs() != 2 {
                    return None;
                }
                let mut argitem = f.iter();
                let i = argitem.next().unwrap();
                let j = argitem.next().unwrap();

                let mut out = "γ".to_string();
                if is_typst {
                    out.push_str("_0");
                } else {
                    out.push_str(&to_superscript(0));
                }
                if symbol_scripts {
                    out.push('^');
                }
                if opt.color_builtin_symbols {
                    out = nu_ansi_term::Color::Magenta.paint(out).to_string();
                }

                if parens {
                    out.push('(');
                }
                i.format(&mut out, opt, PrintState::new()).unwrap();
                if commas {
                    out.push(',');
                } else {
                    out.push(' ');
                }
                j.format(&mut out, opt, PrintState::new()).unwrap();
                if parens {
                    out.push(')');
                }

                Some(out)
            }
            _ => None,
        }
    }),
    gammaconj: symbol!("spenso::gammaconj"),
});

fn spinor_matrix_structure<Aind: AbsInd>(
    symbol: Symbol,
    dim: impl Into<Dimension>,
) -> ExplicitKey<Aind> {
    let dim = dim.into();
    ExplicitKey::from_iter(
        [Bispinor {}.new_rep(dim), Bispinor {}.new_rep(dim)],
        symbol,
        None,
    )
    .structure
}

fn gamma_matrix_structure<Aind: AbsInd>(
    symbol: Symbol,
    dim: impl Into<Dimension>,
) -> ExplicitKey<Aind> {
    ExplicitKey::from_iter(
        [
            LibraryRep::from(Minkowski {}).new_rep(dim),
            Bispinor {}.new_rep(4).cast(),
            Bispinor {}.new_rep(4).cast(),
        ],
        symbol,
        None,
    )
    .structure
}

macro_rules! gamma_matrix_structure_methods {
    ($($structure:ident, $field:ident;)*) => {
        $(
            pub fn $structure<Aind: AbsInd>(&self, dim: impl Into<Dimension>) -> ExplicitKey<Aind> {
                gamma_matrix_structure(self.$field, dim)
            }
        )*
    };
}

impl GammaLibrary {
    gamma_matrix_structure_methods! {
        gamma_strct, gamma;
        gamma_conj_strct, gammaconj;
        gamma_adj_strct, gammaadj;
    }

    pub fn gamma0_strct<Aind: AbsInd>(&self, dim: impl Into<Dimension>) -> ExplicitKey<Aind> {
        spinor_matrix_structure(self.gamma0, dim)
    }

    pub fn projm_strct<Aind: AbsInd>(&self, dim: impl Into<Dimension>) -> ExplicitKey<Aind> {
        spinor_matrix_structure(self.projm, dim)
    }

    pub fn projp_strct<Aind: AbsInd>(&self, dim: impl Into<Dimension>) -> ExplicitKey<Aind> {
        spinor_matrix_structure(self.projp, dim)
    }

    pub fn gamma5_strct<Aind: AbsInd>(&self, dim: impl Into<Dimension>) -> ExplicitKey<Aind> {
        spinor_matrix_structure(self.gamma5, dim)
    }
}

fn collect_gammas(expr: &mut Atom) {
    let reps: Vec<_> = [
        (
            function!(
                AGS.projp,
                Bispinor {}.to_symbolic([RS.a__]),
                Bispinor {}.to_symbolic([RS.b__])
            ),
            (Bispinor {}.id_atom([RS.a__], [RS.b__]) - crate::gamma5!(RS.a__, RS.b__)) / 2,
        ),
        (
            function!(
                AGS.projm,
                Bispinor {}.to_symbolic([RS.a__]),
                Bispinor {}.to_symbolic([RS.b__])
            ),
            (Bispinor {}.id_atom([RS.a__], [RS.b__]) + crate::gamma5!(RS.a__, RS.b__)) / 2,
        ),
        (
            crate::gamma!(RS.a__, RS.b__, RS.c__) * crate::gamma!(RS.d__, RS.c__, RS.e__),
            GS.chain_pattern(RS.b__, RS.e__, [RS.a__, RS.c__, RS.d__]),
        ),
        (crate::gamma!(RS.a__, RS.b__, RS.b__), Atom::Zero),
        (
            function!(GS.gamma_chain, RS.a__, RS.a_, RS.b_)
                * function!(GS.gamma_chain, RS.b__, RS.b_, RS.c_),
            function!(GS.gamma_chain, RS.a__, RS.b_, RS.b__, RS.a_, RS.c_),
        ),
        (
            function!(
                GS.gamma_chain,
                RS.a___,
                Bispinor {}.to_symbolic([RS.a__]),
                Bispinor {}.to_symbolic([RS.b__])
            ) * crate::gamma!(RS.y__, RS.b__, RS.c__),
            function!(
                GS.gamma_chain,
                RS.a___,
                Bispinor {}.to_symbolic([RS.b__]),
                Minkowski {}.to_symbolic([RS.y__]),
                Bispinor {}.to_symbolic([RS.a__]),
                Bispinor {}.to_symbolic([RS.c__])
            ),
        ),
        (
            crate::gamma!(RS.y__, RS.a__, RS.b__)
                * function!(
                    GS.gamma_chain,
                    RS.a___,
                    Bispinor {}.to_symbolic([RS.b__]),
                    Bispinor {}.to_symbolic([RS.c__])
                ),
            function!(
                GS.gamma_chain,
                RS.a___,
                Bispinor {}.to_symbolic([RS.b__]),
                Minkowski {}.to_symbolic([RS.y__]),
                Bispinor {}.to_symbolic([RS.a__]),
                Bispinor {}.to_symbolic([RS.c__])
            ),
        ),
    ]
    .iter()
    .map(|(a, b)| {
        // println!("{:#}->{:#}", a, b);

        Replacement::new(a.to_pattern(), b.to_pattern())
    })
    .collect();

    let mut atom = Atom::new();

    while expr.replace_multiple_into(&reps, &mut atom) {
        // println!("collecting:{atom}");
        std::mem::swap(expr, &mut atom);
        *expr = expr.expand();
        // println!("expanding::{expr}");
        *expr = expr.simplify_metrics();
        // println!("simplifying::{expr}");
    }
}

fn normalise_gammas(expr: &mut Atom) {
    // Uses the anti commutation rule of the gamma chain to sort the minkowski indices
    fn gamma_chain_normalisation(arg: AtomView, _context: &Context, out: &mut Settable<'_, Atom>) {
        if let AtomView::Fun(f) = arg
            && f.get_symbol() == GS.gamma_chain
        {
            // println!("Gamma chain:{}", arg);
            let mut args = f.iter().collect::<Vec<_>>();
            if args.len() >= 5 {
                let more = args.len() > 5;
                for i in 0..args.len().saturating_sub(4) {
                    if i % 2 == 0 {
                        // println!("{}", args[i]);
                        // println!("{}?{}", args[i], args[i + 1]);
                        if args[i] > args[i + 2] {
                            // println!("{}>{}", args[i], args[i + 1]);
                            args.swap(i, i + 2);
                            let swapped = FunctionBuilder::new(GS.gamma_chain)
                                .add_args(&args)
                                .finish();

                            // println!("{}", swapped);
                            let mu = args.remove(i);
                            let _bis = args.remove(i);
                            let nu = args.remove(i);
                            if more {
                                if i == 0 {
                                    args.remove(i);
                                } else {
                                    args.remove(i - 1);
                                }
                            }

                            // println!("mu:{}bis:{}nu:{}", mu, bis, nu);
                            let metric = function!(ETS.metric, mu, nu)
                                * 2
                                * FunctionBuilder::new(GS.gamma_chain)
                                    .add_args(&args)
                                    .finish();
                            **out = metric - swapped;
                            return;
                            // println!("{}->{}", a, c);
                        }
                    }
                }
            }
        }
    }

    loop {
        let new = expr.replace_map(&gamma_chain_normalisation);
        // .replace_multiple(&reps);
        if new == *expr {
            break;
        } else {
            *expr = new;
        }
    }
}

#[cfg(test)]
fn undo_gamma_chain(expr: &mut Atom) {
    let reps: Vec<_> = [
        (
            GS.chain_pattern(RS.b__, RS.e__, [RS.a__, RS.c__, RS.d__]),
            crate::gamma!(RS.a__, RS.b__, RS.c__) * crate::gamma!(RS.d__, RS.c__, RS.e__),
        ),
        (
            function!(
                GS.gamma_chain,
                RS.a___,
                Bispinor {}.to_symbolic([RS.b__]),
                Minkowski {}.to_symbolic([RS.y__]),
                Bispinor {}.to_symbolic([RS.a__]),
                Bispinor {}.to_symbolic([RS.c__])
            ),
            function!(
                GS.gamma_chain,
                RS.a___,
                Bispinor {}.to_symbolic([RS.a__]),
                Bispinor {}.to_symbolic([RS.b__])
            ) * crate::gamma!(RS.y__, RS.b__, RS.c__),
        ),
        (
            function!(
                GS.gamma_chain,
                RS.a___,
                Bispinor {}.to_symbolic([RS.b__]),
                Minkowski {}.to_symbolic([RS.y__]),
                Bispinor {}.to_symbolic([RS.a__]),
                Bispinor {}.to_symbolic([RS.c__])
            ),
            crate::gamma!(RS.y__, RS.a__, RS.b__)
                * function!(
                    GS.gamma_chain,
                    RS.a___,
                    Bispinor {}.to_symbolic([RS.b__]),
                    Bispinor {}.to_symbolic([RS.c__])
                ),
        ),
    ]
    .iter()
    .map(|(a, b)| Replacement::new(a.to_pattern(), b.to_pattern()))
    .collect();

    loop {
        let new = expr.replace_multiple(&reps);
        if new == *expr {
            break;
        }
        *expr = new;
    }
}

pub fn gamma_simplify_impl(expr: AtomView) -> Atom {
    let mink = Minkowski {};

    let mut coef_list_spin = expr.expand_mink_bis();

    for (expr, _) in &mut coef_list_spin {
        // println!("gs:{expr}");
        *expr = expr.simplify_metrics();
        collect_gammas(expr);
        // println!("collected:{expr}");
        // expr = normalise_gammas(expr);

        let reps: Vec<_> = [
            (
                function!(
                    GS.gamma_chain,
                    RS.a___,
                    Bispinor {}.to_symbolic([RS.i__]),
                    Minkowski {}.to_symbolic([RS.d_, RS.a_]),
                    Bispinor {}.to_symbolic([RS.j__]),
                    Minkowski {}.to_symbolic([RS.d_, RS.a_]),
                    RS.b__
                ),
                function!(GS.gamma_chain, RS.a___, RS.b__) * RS.d_,
            ),
            (
                function!(
                    GS.gamma_chain,
                    Minkowski {}.to_symbolic([RS.d_, RS.a_]),
                    Bispinor {}.to_symbolic([RS.b__]),
                    Minkowski {}.to_symbolic([RS.d_, RS.a_]),
                    Bispinor {}.to_symbolic([RS.i__]),
                    Bispinor {}.to_symbolic([RS.j__])
                ),
                function!(
                    GS.gamma_chain,
                    Bispinor {}.to_symbolic([RS.i__]),
                    Bispinor {}.to_symbolic([RS.j__])
                ) * RS.d_,
            ),
            (
                function!(
                    GS.gamma_chain,
                    Minkowski {}.to_symbolic([RS.a__]),
                    Bispinor {}.to_symbolic([RS.i__]),
                    Bispinor {}.to_symbolic([RS.j__])
                ),
                function!(
                    AGS.gamma,
                    Bispinor {}.to_symbolic([RS.i__]),
                    Bispinor {}.to_symbolic([RS.j__]),
                    Minkowski {}.to_symbolic([RS.a__])
                ),
            ),
            (
                function!(
                    AGS.gamma,
                    Bispinor {}.to_symbolic([RS.i__]),
                    Bispinor {}.to_symbolic([RS.i__]),
                    Minkowski {}.to_symbolic([RS.a__])
                ),
                Atom::Zero,
            ),
            (
                function!(
                    GS.gamma_chain,
                    Bispinor {}.to_symbolic([RS.a__]),
                    Bispinor {}.to_symbolic([RS.b__])
                ),
                Bispinor {}.id_atom([RS.a__], [RS.b__]),
            ),
        ]
        .iter()
        .map(|(a, b)| {
            // println!("{a:#}->\n\t{b:#}");
            Replacement::new(a.to_pattern(), b.to_pattern())
        })
        .collect();

        normalise_gammas(expr);
        loop {
            let new = expr.replace_multiple(&reps).expand().simplify_metrics();
            if new == *expr {
                break;
            } else {
                *expr = new;
            }
        }

        *expr = expr
            .replace(function!(GS.gamma_chain, RS.a__, RS.x_, RS.x_).to_pattern())
            .repeat()
            .with(function!(GS.gamma_trace, RS.a__).to_pattern())
            .replace(function!(
                GS.gamma_trace,
                RS.a___,
                Bispinor {}.to_symbolic([RS.a__]),
                RS.b___
            ))
            .repeat()
            .with(function!(GS.gamma_trace, RS.a___, RS.b___));

        // undo_gamma_chain(expr);
        // println!("Before tracer:{expr}");

        // expr = expr.replace_map(|term, ctx, out| {

        //     if let AtomView::Fun(f)= term{
        //         if f.get_symbol()= GS.gamma_chain{

        //         }
        //     }

        //     false});

        // //Chisholm identity:
        // expr.replace_all_repeat_mut(
        //     &(function!(AGS.gamma, RS.a_, RS.x_, RS.y_) * function!(gamma_trace, RS.a_, RS.a__)).to_pattern(),
        //     (function!(gamma_chain, RS.a__)).to_pattern(),
        //     None,
        //     None,
        // );
        //
        fn gamma_tracer(arg: AtomView, _context: &Context, out: &mut Settable<'_, Atom>) {
            let gamma_trace = GS.gamma_trace;

            if let AtomView::Fun(f) = arg
                && f.get_symbol() == gamma_trace
            {
                // println!("{arg}");

                let mut sum = Atom::Zero;

                if f.get_nargs() == 1 {
                    **out = Atom::Zero;
                }
                let args = f.iter().collect::<Vec<_>>();

                for i in 1..args.len() {
                    let sign = if i % 2 == 0 { -1 } else { 1 };

                    let mut gcn = FunctionBuilder::new(gamma_trace);
                    #[allow(clippy::needless_range_loop)]
                    for j in 1..args.len() {
                        if i != j {
                            gcn = gcn.add_arg(args[j]);
                        }
                    }

                    let metric = if args[0] == args[i] {
                        if let AtomView::Fun(f) = args[0].as_atom_view() {
                            f.iter().next().unwrap().to_owned()
                        } else {
                            panic!("aaaa")
                        }
                        // Atom::num(4)
                    } else {
                        function!(ETS.metric, args[0], args[i])
                    };
                    if args.len() == 2 {
                        sum += metric * sign * Atom::num(4);
                    } else {
                        sum += metric * gcn.finish() * sign;
                    }
                }
                **out = sum;

                // println!("{}->{}", arg, out);
            }
        }

        loop {
            let new = expr.replace_map(&gamma_tracer);
            if new == *expr {
                break;
            } else {
                *expr = new;
            }
        }

        *expr = expr
            .replace(
                function!(AGS.gamma, RS.a__, mink.to_symbolic([RS.d_, RS.b_]))
                    .pow(Atom::num(2))
                    .to_pattern(),
            )
            .repeat()
            .with(Atom::var(RS.d_) * 4)
            .expand()
            .simplify_metrics();
    }

    coef_list_spin
        .iter()
        .fold(Atom::Zero, |a, (c, s)| a + c * s)
}

/// Trait for simplifying expressions involving Dirac gamma matrices using Clifford algebra.
///
/// Implementors provide a method to apply gamma matrix identities, such as
/// anticommutation relations and trace evaluations.
pub trait GammaSimplifier {
    /// Simplifies gamma matrix structures within the expression.
    ///
    /// Uses the Clifford algebra relation `{gamma^mu, gamma^nu} = 2 * g^{mu nu}`
    /// and evaluates traces of products of gamma matrices. It handles intermediate
    /// simplification steps involving metric tensors.
    ///
    /// # Returns
    /// An [`Atom`] representing the expression after gamma matrix simplification.
    fn simplify_gamma(&self) -> Atom;

    /// Simplifies gamma matrices with explicit chain-ordering and trace settings.
    fn simplify_gamma_with(&self, settings: GammaSimplifySettings) -> Atom;

    fn simplify_gamma0(&self) -> Atom;

    fn simplify_gamma_conj<Aind: DummyAind + ParseableAind>(&self) -> eyre::Result<Atom>;
}

impl GammaSimplifier for Atom {
    fn simplify_gamma(&self) -> Atom {
        self.as_view().simplify_gamma()
    }

    fn simplify_gamma_with(&self, settings: GammaSimplifySettings) -> Atom {
        self.as_view().simplify_gamma_with(settings)
    }

    fn simplify_gamma0(&self) -> Atom {
        self.as_view().simplify_gamma0()
    }

    fn simplify_gamma_conj<Aind: DummyAind + ParseableAind>(&self) -> eyre::Result<Atom> {
        self.as_view().simplify_gamma_conj::<Aind>()
    }
}

impl GammaSimplifier for AtomView<'_> {
    fn simplify_gamma(&self) -> Atom {
        self.simplify_gamma_with(GammaSimplifySettings::default())
    }

    fn simplify_gamma_with(&self, settings: GammaSimplifySettings) -> Atom {
        simplify::simplify_dirac_chains_impl(*self, settings)
    }

    fn simplify_gamma0(&self) -> Atom {
        let repeated_gamma0 = crate::gamma0!(RS.a__, RS.b__) * crate::gamma0!(RS.b__, RS.c__);

        let gamma0_ia = crate::gamma0!([RS.d_, RS.i_], [RS.d_, RS.a_]);
        let gamma_ab = crate::gamma!(RS.a__, [RS.d_, RS.a_], [RS.d_, RS.b_]);
        let gamma0_bj = crate::gamma0!([RS.d_, RS.b_], [RS.d_, RS.j_]);

        let gmg = (Atom::var(RS.f_) * gamma0_ia.clone() * gamma_ab.clone() * gamma0_bj.clone()
            + Atom::var(RS.e_) * Bispinor {}.metric_atom([RS.d_, RS.j_], [RS.d_, RS.i_]))
        .to_pattern();

        let gmgrhs = (gamma0_ia.clone()
            * (Atom::var(RS.f_) * gamma_ab.clone()
                + Atom::var(RS.e_) * Bispinor {}.metric_atom([RS.d_, RS.a_], [RS.d_, RS.b_]))
            * gamma0_bj.clone())
        .to_pattern();

        let gmgn = (Atom::var(RS.f_) * gamma0_ia.clone() * gamma_ab.clone() * gamma0_bj.clone()
            + Bispinor {}.metric_atom([RS.d_, RS.j_], [RS.d_, RS.i_]))
        .to_pattern();

        let gmgnrhs = (gamma0_ia
            * (Atom::var(RS.f_) * gamma_ab
                + Bispinor {}.metric_atom([RS.d_, RS.a_], [RS.d_, RS.b_]))
            * gamma0_bj)
            .to_pattern();

        self.replace(gmg)
            .with(gmgrhs)
            .replace(gmgn)
            .with(gmgnrhs)
            .replace(repeated_gamma0)
            .repeat()
            .with(Bispinor {}.metric_atom([RS.a__], [RS.c__]))
    }

    fn simplify_gamma_conj<Aind: DummyAind + ParseableAind>(&self) -> Result<Atom> {
        let dummy = symbol!("dummy");

        let dummypati = function!(dummy, RS.i_).to_pattern();
        let dummypatj = function!(dummy, RS.j_).to_pattern();

        let conj_gamma = crate::gamma!(RS.a__, [RS.d_, RS.i_], [RS.d_, RS.j_]).spenso_conj();

        let conj_gamma_rhs =
            (crate::gamma0!([RS.d_, RS.j_], [Atom::var(RS.d_), function!(dummy, RS.j_)])
                * crate::gamma!(
                    RS.a__,
                    [Atom::var(RS.d_), function!(dummy, RS.j_)],
                    [Atom::var(RS.d_), function!(dummy, RS.i_)]
                )
                * crate::gamma0!([Atom::var(RS.d_), function!(dummy, RS.i_)], [RS.d_, RS.i_]))
            .to_pattern();

        Ok(self.replace(conj_gamma).with_map(move |m| {
            let a = conj_gamma_rhs.replace_wildcards_with_matches(m);
            let i = dummypati.replace_wildcards_with_matches(m);
            let j = dummypatj.replace_wildcards_with_matches(m);
            a.replace(i)
                .with(Aind::new_dummy().to_atom())
                .replace(j)
                .with(Aind::new_dummy().to_atom())
        }))
    }
}

pub fn id_atom(i: impl Into<Atom>, j: impl Into<Atom>) -> Atom {
    function!(ETS.metric, i.into(), j.into())
}

mod macros;
mod simplify;
pub use simplify::{GammaChainOrdering, GammaSimplifySettings};
#[cfg(test)]
mod test;
