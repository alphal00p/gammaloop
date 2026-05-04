use std::sync::LazyLock;

use spenso::{
    network::{
        library::symbolic::{ETS, ExplicitKey},
        tags::SPENSO_TAG as T,
    },
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

    expr.replace(AGS.gamma_pattern(RS.a__, RS.i__, RS.j__).to_pattern())
        .with((-AGS.gamma_pattern(RS.a__, RS.j__, RS.i__)).to_pattern())
        .replace(AGS.gamma5_pattern(RS.i__, RS.j__).to_pattern())
        .with(AGS.gamma5_pattern(RS.j__, RS.i__).to_pattern())
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

impl GammaLibrary {
    pub fn chain_gamma<'a>(&self, mu: impl Into<AtomOrView<'a>>) -> Atom {
        FunctionBuilder::new(self.gamma)
            .add_arg(Atom::var(T.chain_in))
            .add_arg(Atom::var(T.chain_out))
            .add_arg(mu)
            .finish()
    }

    pub fn chain_gamma5(&self) -> Atom {
        FunctionBuilder::new(self.gamma5)
            .add_arg(Atom::var(T.chain_in))
            .add_arg(Atom::var(T.chain_out))
            .finish()
    }

    pub fn projp<'a, 'b>(
        &self,
        a: impl Into<AtomOrView<'a>>,
        b: impl Into<AtomOrView<'b>>,
    ) -> Atom {
        function!(self.projp, a.into().as_view(), b.into().as_view())
    }

    pub fn projp_pattern<'a>(
        &self,
        i: impl Into<AtomOrView<'a>>,
        j: impl Into<AtomOrView<'a>>,
    ) -> Atom {
        function!(
            self.projp,
            Bispinor {}.to_symbolic([i]),
            Bispinor {}.to_symbolic([j])
        )
    }

    pub fn projm<'a, 'b>(
        &self,
        a: impl Into<AtomOrView<'a>>,
        b: impl Into<AtomOrView<'b>>,
    ) -> Atom {
        function!(self.projm, a.into().as_view(), b.into().as_view())
    }
    pub fn projm_pattern<'a>(
        &self,
        i: impl Into<AtomOrView<'a>>,
        j: impl Into<AtomOrView<'a>>,
    ) -> Atom {
        function!(
            self.projm,
            Bispinor {}.to_symbolic([i]),
            Bispinor {}.to_symbolic([j])
        )
    }

    pub fn gamma_strct<Aind: AbsInd>(&self, dim: impl Into<Dimension>) -> ExplicitKey<Aind> {
        let gamma = ExplicitKey::from_iter(
            [
                LibraryRep::from(Minkowski {}).new_rep(dim),
                Bispinor {}.new_rep(4).cast(),
                Bispinor {}.new_rep(4).cast(),
            ],
            self.gamma,
            None,
        );
        gamma.structure
    }

    pub fn gamma_conj_strct<Aind: AbsInd>(&self, dim: impl Into<Dimension>) -> ExplicitKey<Aind> {
        let gamma = ExplicitKey::from_iter(
            [
                LibraryRep::from(Minkowski {}).new_rep(dim),
                Bispinor {}.new_rep(4).cast(),
                Bispinor {}.new_rep(4).cast(),
            ],
            self.gammaconj,
            None,
        );
        gamma.structure
    }

    pub fn gamma_adj_strct<Aind: AbsInd>(&self, dim: impl Into<Dimension>) -> ExplicitKey<Aind> {
        let gamma = ExplicitKey::from_iter(
            [
                LibraryRep::from(Minkowski {}).new_rep(dim),
                Bispinor {}.new_rep(4).cast(),
                Bispinor {}.new_rep(4).cast(),
            ],
            self.gammaadj,
            None,
        );
        gamma.structure
    }

    pub fn gamma_pattern<'a>(
        &self,
        mu: impl Into<AtomOrView<'a>>,
        i: impl Into<AtomOrView<'a>>,
        j: impl Into<AtomOrView<'a>>,
    ) -> Atom {
        function!(
            self.gamma,
            Bispinor {}.to_symbolic([i]),
            Bispinor {}.to_symbolic([j]),
            Minkowski {}.to_symbolic([mu])
        )
    }
    pub fn gamma5_pattern<'a>(
        &self,
        i: impl Into<AtomOrView<'a>>,
        j: impl Into<AtomOrView<'a>>,
    ) -> Atom {
        function!(
            self.gamma5,
            Bispinor {}.to_symbolic([i]),
            Bispinor {}.to_symbolic([j])
        )
    }

    pub fn gamma0_pattern<'a>(
        &self,
        i: impl Into<AtomOrView<'a>>,
        j: impl Into<AtomOrView<'a>>,
    ) -> Atom {
        function!(
            self.gamma0,
            Bispinor {}.to_symbolic([i]),
            Bispinor {}.to_symbolic([j])
        )
    }

    pub fn gamma5_strct<Aind: AbsInd>(&self, dim: impl Into<Dimension>) -> ExplicitKey<Aind> {
        let dim = dim.into();
        let gamma5 = ExplicitKey::from_iter(
            [Bispinor {}.new_rep(dim), Bispinor {}.new_rep(dim)],
            self.gamma5,
            None,
        );
        gamma5.structure
    }

    pub fn gamma0_strct<Aind: AbsInd>(&self, dim: impl Into<Dimension>) -> ExplicitKey<Aind> {
        let dim = dim.into();
        let gamma5 = ExplicitKey::from_iter(
            [Bispinor {}.new_rep(dim), Bispinor {}.new_rep(dim)],
            self.gamma0,
            None,
        );
        gamma5.structure
    }

    pub fn projm_strct<Aind: AbsInd>(&self, dim: impl Into<Dimension>) -> ExplicitKey<Aind> {
        let dim = dim.into();
        let projm = ExplicitKey::from_iter(
            [Bispinor {}.new_rep(dim), Bispinor {}.new_rep(dim)],
            self.projm,
            None,
        );
        projm.structure
    }

    pub fn projp_strct<Aind: AbsInd>(&self, dim: impl Into<Dimension>) -> ExplicitKey<Aind> {
        let dim = dim.into();
        let projp_strct = ExplicitKey::from_iter(
            [Bispinor {}.new_rep(dim), Bispinor {}.new_rep(dim)],
            self.projp,
            None,
        );

        projp_strct.structure
    }
}

fn collect_gammas(expr: &mut Atom) {
    let reps: Vec<_> = [
        (
            AGS.projp_pattern(RS.a__, RS.b__),
            (Bispinor {}.id_atom([RS.a__], [RS.b__]) - AGS.gamma5_pattern(RS.a__, RS.b__)) / 2,
        ),
        (
            AGS.projm_pattern(RS.a__, RS.b__),
            (Bispinor {}.id_atom([RS.a__], [RS.b__]) + AGS.gamma5_pattern(RS.a__, RS.b__)) / 2,
        ),
        (
            AGS.gamma_pattern(RS.a__, RS.b__, RS.c__) * AGS.gamma_pattern(RS.d__, RS.c__, RS.e__),
            GS.chain_pattern(RS.b__, RS.e__, [RS.a__, RS.c__, RS.d__]),
        ),
        (AGS.gamma_pattern(RS.a__, RS.b__, RS.b__), Atom::Zero),
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
            ) * AGS.gamma_pattern(RS.y__, RS.b__, RS.c__),
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
            AGS.gamma_pattern(RS.y__, RS.a__, RS.b__)
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

fn undo_gamma_chain(expr: &mut Atom) {
    let reps: Vec<_> = [
        (
            GS.chain_pattern(RS.b__, RS.e__, [RS.a__, RS.c__, RS.d__]),
            AGS.gamma_pattern(RS.a__, RS.b__, RS.c__) * AGS.gamma_pattern(RS.d__, RS.c__, RS.e__),
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
            ) * AGS.gamma_pattern(RS.y__, RS.b__, RS.c__),
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
            AGS.gamma_pattern(RS.y__, RS.a__, RS.b__)
                * function!(
                    GS.gamma_chain,
                    RS.a___,
                    Bispinor {}.to_symbolic([RS.b__]),
                    Bispinor {}.to_symbolic([RS.c__])
                ),
        ),
    ]
    .iter()
    .map(|(a, b)| {
        // println!("{}->{}", a, b);

        Replacement::new(a.to_pattern(), b.to_pattern())
    })
    .collect();

    loop {
        let new = expr.replace_multiple(&reps);
        if new == *expr {
            break;
        } else {
            *expr = new;
        }
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

        undo_gamma_chain(expr);
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

    fn simplify_gamma0(&self) -> Atom;

    fn simplify_gamma_conj<Aind: DummyAind + ParseableAind>(&self) -> eyre::Result<Atom>;
}

impl GammaSimplifier for Atom {
    fn simplify_gamma(&self) -> Atom {
        gamma_simplify_impl(self.as_atom_view())
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
        gamma_simplify_impl(self.as_atom_view())
    }

    fn simplify_gamma0(&self) -> Atom {
        let repeated_gamma0 =
            AGS.gamma0_pattern(RS.a__, RS.b__) * AGS.gamma0_pattern(RS.b__, RS.c__);

        let gmg = (Atom::var(RS.f_)
            * function!(
                AGS.gamma0,
                Bispinor {}.to_symbolic([RS.d_, RS.i_]),
                Bispinor {}.to_symbolic([RS.d_, RS.a_])
            )
            * function!(
                AGS.gamma,
                Bispinor {}.to_symbolic([RS.d_, RS.a_]),
                Bispinor {}.to_symbolic([RS.d_, RS.b_]),
                Minkowski {}.to_symbolic([RS.a__])
            )
            * function!(
                AGS.gamma0,
                Bispinor {}.to_symbolic([RS.d_, RS.b_]),
                Bispinor {}.to_symbolic([RS.d_, RS.j_])
            )
            + Atom::var(RS.e_) * Bispinor {}.metric_atom([RS.d_, RS.j_], [RS.d_, RS.i_]))
        .to_pattern();

        let gmgrhs = (function!(
            AGS.gamma0,
            Bispinor {}.to_symbolic([RS.d_, RS.i_]),
            Bispinor {}.to_symbolic([RS.d_, RS.a_])
        ) * (Atom::var(RS.f_)
            * function!(
                AGS.gamma,
                Bispinor {}.to_symbolic([RS.d_, RS.a_]),
                Bispinor {}.to_symbolic([RS.d_, RS.b_]),
                Minkowski {}.to_symbolic([RS.a__])
            )
            + Atom::var(RS.e_) * Bispinor {}.metric_atom([RS.d_, RS.a_], [RS.d_, RS.b_]))
            * function!(
                AGS.gamma0,
                Bispinor {}.to_symbolic([RS.d_, RS.b_]),
                Bispinor {}.to_symbolic([RS.d_, RS.j_])
            ))
        .to_pattern();

        let gmgn = (Atom::var(RS.f_)
            * function!(
                AGS.gamma0,
                Bispinor {}.to_symbolic([RS.d_, RS.i_]),
                Bispinor {}.to_symbolic([RS.d_, RS.a_])
            )
            * function!(
                AGS.gamma,
                Bispinor {}.to_symbolic([RS.d_, RS.a_]),
                Bispinor {}.to_symbolic([RS.d_, RS.b_]),
                Minkowski {}.to_symbolic([RS.a__])
            )
            * function!(
                AGS.gamma0,
                Bispinor {}.to_symbolic([RS.d_, RS.b_]),
                Bispinor {}.to_symbolic([RS.d_, RS.j_])
            )
            + Bispinor {}.metric_atom([RS.d_, RS.j_], [RS.d_, RS.i_]))
        .to_pattern();

        let gmgnrhs = (function!(
            AGS.gamma0,
            Bispinor {}.to_symbolic([RS.d_, RS.i_]),
            Bispinor {}.to_symbolic([RS.d_, RS.a_])
        ) * (Atom::var(RS.f_)
            * function!(
                AGS.gamma,
                Bispinor {}.to_symbolic([RS.d_, RS.a_]),
                Bispinor {}.to_symbolic([RS.d_, RS.b_]),
                Minkowski {}.to_symbolic([RS.a__])
            )
            + Bispinor {}.metric_atom([RS.d_, RS.a_], [RS.d_, RS.b_]))
            * function!(
                AGS.gamma0,
                Bispinor {}.to_symbolic([RS.d_, RS.b_]),
                Bispinor {}.to_symbolic([RS.d_, RS.j_])
            ))
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

        let conj_gamma = function!(
            AGS.gamma,
            Bispinor {}.to_symbolic([RS.d_, RS.i_]),
            Bispinor {}.to_symbolic([RS.d_, RS.j_]),
            Minkowski {}.to_symbolic([RS.a__])
        )
        .spenso_conj();

        let conj_gamma_rhs = (function!(
            AGS.gamma0,
            Bispinor {}.to_symbolic([RS.d_, RS.j_]),
            Bispinor {}.to_symbolic([Atom::var(RS.d_), function!(dummy, RS.j_)])
        ) * function!(
            AGS.gamma,
            Bispinor {}.to_symbolic([Atom::var(RS.d_), function!(dummy, RS.j_)]),
            Bispinor {}.to_symbolic([Atom::var(RS.d_), function!(dummy, RS.i_)]),
            Minkowski {}.to_symbolic([RS.a__])
        ) * function!(
            AGS.gamma0,
            Bispinor {}.to_symbolic([Atom::var(RS.d_), function!(dummy, RS.i_)]),
            Bispinor {}.to_symbolic([RS.d_, RS.i_])
        ))
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

/// Builds an identity atom between two indices.
///
/// The arguments are parsed as Symbolica literals, preserving the historical
/// `id!(i, j)` shorthand used in idenso tests and examples.
#[macro_export]
macro_rules! id {
    ($i: expr, $j: expr) => {{
        let i = symbolica::parse_lit!($i);
        let j = symbolica::parse_lit!($j);
        id_atom(i, j)
    }};
}

/// Builds a gamma matrix.
///
/// With one argument, this builds a chain factor using the placeholder indices
/// `in` and `out`; use this form only as a factor inside `chain!` or `trace!`.
/// With three arguments, this builds the ordinary gamma tensor with explicit
/// spinor endpoints and a Lorentz slot.
///
/// Arguments are converted through `spenso::symbolica_atom::IntoAtom`, so they
/// can be typed slots, atoms, or atom views.
///
/// # Examples
///
/// ```ignore
/// use idenso::gamma;
/// use spenso::{chain, slot};
///
/// let factor = gamma!(slot!(mink4, mu));
/// let chain_expr = chain!(slot!(bis4, a), slot!(bis4, b), factor);
/// let default_tensor = gamma!(mu, a, b);
/// let indexed_tensor = gamma!(1, 2, 3);
/// let mixed_tensor = gamma!(mu, slot!(bis_d, a), 1);
/// let explicit_tensor = gamma!(slot!(mink_d, mu), slot!(bis_d, a), slot!(bis_d, b));
/// ```
#[macro_export]
macro_rules! gamma {
    ($mu:expr) => {
        $crate::gamma::AGS.chain_gamma(spenso::symbolica_atom::IntoAtom::into_atom($mu))
    };
    ($mu:ident, $($rest:tt)*) => {{
        let mink = spenso::structure::representation::RepName::new_rep(
            &spenso::structure::representation::Minkowski {},
            4,
        );
        $crate::gamma!(@tensor spenso::slot!(mink, $mu); $($rest)*)
    }};
    ($mu:literal, $($rest:tt)*) => {{
        let mink = spenso::structure::representation::RepName::new_rep(
            &spenso::structure::representation::Minkowski {},
            4,
        );
        $crate::gamma!(@tensor spenso::slot!(mink, $mu); $($rest)*)
    }};
    ($mu:expr, $($rest:tt)*) => {
        $crate::gamma!(@tensor $mu; $($rest)*)
    };
    (@tensor $mu:expr; $i:ident, $($rest:tt)*) => {{
        let bis = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::Bispinor {},
            4,
        );
        $crate::gamma!(@tensor2 $mu, spenso::slot!(bis, $i); $($rest)*)
    }};
    (@tensor $mu:expr; $i:literal, $($rest:tt)*) => {{
        let bis = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::Bispinor {},
            4,
        );
        $crate::gamma!(@tensor2 $mu, spenso::slot!(bis, $i); $($rest)*)
    }};
    (@tensor $mu:expr; $i:expr, $($rest:tt)*) => {
        $crate::gamma!(@tensor2 $mu, $i; $($rest)*)
    };
    (@tensor2 $mu:expr, $i:expr; $j:ident) => {{
        let bis = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::Bispinor {},
            4,
        );
        $crate::gamma!(@tensor_done $mu, $i, spenso::slot!(bis, $j))
    }};
    (@tensor2 $mu:expr, $i:expr; $j:literal) => {{
        let bis = spenso::structure::representation::RepName::new_rep(
            &$crate::representations::Bispinor {},
            4,
        );
        $crate::gamma!(@tensor_done $mu, $i, spenso::slot!(bis, $j))
    }};
    (@tensor2 $mu:expr, $i:expr; $j:expr) => {
        $crate::gamma!(@tensor_done $mu, $i, $j)
    };
    (@tensor_done $mu:expr, $i:expr, $j:expr) => {
        symbolica::atom::FunctionBuilder::new($crate::gamma::AGS.gamma)
            .add_arg(spenso::symbolica_atom::IntoAtom::into_atom($i))
            .add_arg(spenso::symbolica_atom::IntoAtom::into_atom($j))
            .add_arg(spenso::symbolica_atom::IntoAtom::into_atom($mu))
            .finish()
    };
}

/// Builds a gamma-five factor for use inside `chain!` or `trace!`.
///
/// The factor uses the chain placeholder indices `in` and `out`; the surrounding
/// chain or trace owns the physical endpoints.
#[macro_export]
macro_rules! gamma5 {
    () => {
        $crate::gamma::AGS.chain_gamma5()
    };
}
#[cfg(test)]
mod test;
