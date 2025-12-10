use std::sync::LazyLock;

use crate::graph::parse::string_utils::ToOrderedSimple;
use crate::utils::{GS, W_, symbolica_ext::CallSymbol};
use idenso::color::CS;
use idenso::gamma::AGS;
use idenso::representations::{
    Bispinor, ColorAdjoint, ColorAntiFundamental, ColorAntiSextet, ColorFundamental, ColorSextet,
};
use linnet::half_edge::involution::{EdgeIndex, Flow};
use spenso::network::library::symbolic::ETS;
use spenso::structure::representation::{LibraryRep, Minkowski, RepName};
use spenso::structure::slot::{IsAbstractSlot, Slot};
use spenso::structure::{OrderedStructure, TensorStructure};
use symbolica::atom::FunctionBuilder;
use symbolica::atom::{AtomCore, AtomView};
use symbolica::domains::float::FloatLike;
use symbolica::domains::rational::Rational;
use symbolica::symbol;
use symbolica::{
    atom::{Atom, Symbol},
    function,
    id::Replacement,
};
use tracing::debug;

use super::aind::Aind;

#[allow(dead_code)]
pub struct UFOSymbols {
    pub identity: Symbol,
    pub identityl: Symbol,
    pub gamma: Symbol,
    pub gamma5: Symbol,
    pub projm: Symbol,
    pub projp: Symbol,
    pub sigma: Symbol,
    pub charge_conj: Symbol,
    pub metric: Symbol,
    pub momentum: Symbol,
    pub levicivita: Symbol,
    pub t: Symbol,
    pub f: Symbol,
    pub d: Symbol,
    pub antilevicivita: Symbol,
    pub t6: Symbol,
    pub k6: Symbol,
    pub idx: Symbol,
    pub dummy: Symbol,
    pub k6bar: Symbol,
    pub complexconjugate: Symbol,
    pub pslash: Symbol,
    pub complex: Symbol,
}

pub static UFO: LazyLock<UFOSymbols> = LazyLock::new(|| UFOSymbols {
    identity: symbol!("UFO::Identity"),
    identityl: symbol!("UFO::IdentityL"),
    gamma: symbol!("UFO::Gamma"),
    complexconjugate: symbol!(
        "UFO::complexconjugate",
        norm = |f, out| {
            if let AtomView::Fun(ff) = f
                && ff.get_nargs() == 1
            {
                **out = ff.iter().next().unwrap().conj();
            }
        }
    ),
    idx: symbol!("UFO::idx"),
    dummy: symbol!("UFO::dummy"),
    gamma5: symbol!("UFO::Gamma5"),
    projm: symbol!("UFO::ProjM"),
    projp: symbol!("UFO::ProjP"),
    sigma: symbol!("UFO::Sigma"),
    charge_conj: symbol!("UFO::C"),
    metric: symbol!("UFO::Metric"),
    momentum: symbol!("UFO::P"),
    levicivita: symbol!("UFO::Epsilon"),
    t: symbol!("UFO::T"),
    f: symbol!("UFO::f"),
    d: symbol!("UFO::d"),
    antilevicivita: symbol!("UFO::EpsilonBar"),
    t6: symbol!("UFO::T6"),
    k6: symbol!("UFO::K6"),
    k6bar: symbol!("UFO::K6Bar"),
    pslash: symbol!("UFO::PSlash"),
    complex: symbol!(
        "UFO::complex",
        norm = |f, out| {
            if let AtomView::Fun(ff) = f {
                let mut re = Rational::zero();
                let mut im = Rational::zero();
                let mut count = 0;
                for i in ff.iter() {
                    if let Ok(i) = Rational::try_from(i) {
                        if count == 0 {
                            re = i;
                        } else {
                            im = i;
                        }
                        count += 1;
                    } else if let Ok(i) = i64::try_from(i) {
                        if count == 0 {
                            re = re.from_i64(i);
                        } else {
                            im = im.from_i64(i);
                        }
                        count += 1;
                    }

                    if count > 1 {
                        break;
                    }
                }

                if count == 2 {
                    **out = Atom::num(symbolica::domains::float::Complex { re, im });
                }
            }
        }
    ),
});

impl UFOSymbols {
    pub fn idx(&self, id: usize, shift: usize) -> Atom {
        function!(self.idx, shift, id)
    }
    pub(crate) fn reindex_spin(
        &self,
        slots: &[&OrderedStructure<LibraryRep, Aind>],
        momenta: &[(Flow, EdgeIndex)],
        mut atom: Atom,
        dummy: impl Fn(usize) -> Aind,
    ) -> Atom {
        let mink: LibraryRep = Minkowski {}.into();
        let bis: LibraryRep = Bispinor {}.into();

        // for s in slots {
        //     println!("{}", s);
        // }
        // println!("In{}", atom.to_ordered_simple());

        //`P(i,j)` = $$p_j^{\mu_i}$$
        // and P(i)= $$p^{\mu_i}$$
        atom = atom
            .replace(self.momentum.f(&[W_.a_, W_.b___]))
            .with(self.momentum.f((&[W_.b___], &[mink.to_symbolic([W_.a_])])))
            .replace(self.momentum.f((&[self.idx.f((&[1], &[W_.a_]))], &[W_.b_])))
            .with(self.momentum.f(&[W_.a_, W_.b_]));

        // Fill in representations
        let reps: Vec<_> = [
            (
                self.gamma.f(&[W_.a_, W_.i_, W_.j_]),
                self.gamma.f(&[
                    mink.to_symbolic([W_.a_]),
                    bis.to_symbolic([W_.i_]),
                    bis.to_symbolic([W_.j_]),
                ]),
            ),
            (
                self.identityl.f(&[W_.a_, W_.b_]),
                self.identityl
                    .f(&[mink.to_symbolic([W_.a_]), mink.to_symbolic([W_.b_])]),
            ),
            (
                self.metric.f(&[W_.a_, W_.b_]),
                self.metric
                    .f(&[mink.to_symbolic([W_.a_]), mink.to_symbolic([W_.b_])]),
            ),
            (
                self.sigma.f(&[W_.a_, W_.b_, W_.i_, W_.j_]),
                self.sigma.f(&[
                    mink.to_symbolic([W_.a_]),
                    mink.to_symbolic([W_.b_]),
                    bis.to_symbolic([W_.i_]),
                    bis.to_symbolic([W_.j_]),
                ]),
            ),
            (
                self.pslash.f(&[W_.i_, W_.j_, W_.a___]),
                self.pslash.f((
                    &[bis.to_symbolic([W_.i_]), bis.to_symbolic([W_.j_])],
                    &[W_.a___],
                )),
            ),
            (
                self.identity.f(&[W_.i_, W_.j_]),
                self.identity
                    .f(&[bis.to_symbolic([W_.i_]), bis.to_symbolic([W_.j_])]),
            ),
            (
                self.gamma5.f(&[W_.i_, W_.j_]),
                self.gamma5
                    .f(&[bis.to_symbolic([W_.i_]), bis.to_symbolic([W_.j_])]),
            ),
            (
                self.projm.f(&[W_.i_, W_.j_]),
                self.projm
                    .f(&[bis.to_symbolic([W_.i_]), bis.to_symbolic([W_.j_])]),
            ),
            (
                self.projp.f(&[W_.i_, W_.j_]),
                self.projp
                    .f(&[bis.to_symbolic([W_.i_]), bis.to_symbolic([W_.j_])]),
            ),
            (
                self.charge_conj.f(&[W_.i_, W_.j_]),
                self.charge_conj
                    .f(&[bis.to_symbolic([W_.i_]), bis.to_symbolic([W_.j_])]),
            ),
        ]
        .into_iter()
        .map(|(pat, rep)| Replacement::new(pat.to_pattern(), rep))
        .collect();

        atom = atom
            .replace_multiple(&reps)
            .replace(self.levicivita.f((&[W_.a___], &[W_.i_], &[W_.b___])))
            .repeat()
            .with(
                self.levicivita
                    .f((&[W_.a___], &[mink.to_symbolic([W_.i_])], &[W_.b___])),
            );

        for (i, s) in slots.iter().enumerate() {
            let i = i + 1;
            for (shift, r) in s.external_structure_iter().enumerate() {
                let rep = r.rep_name();
                let wrappedi = self.idx(i, shift + 1);
                let wrappedi = wrappedi.as_view();

                // println!("{}", rep.to_symbolic([i]));
                atom = atom
                    .replace(rep.to_symbolic([wrappedi]))
                    .level_range((1, Some(1)))
                    .with(r.to_atom())
                    .replace(rep.to_symbolic([i]))
                    .level_range((1, Some(1)))
                    .with(r.to_atom())
            }
        }
        // debug!(in = atom.printer(LOGPRINTOPTS).to_string());

        // Handle dummies:

        let mut max_dummy = 0;

        // debug!("Before dummies : {}", atom.to_ordered_simple());

        for rep in [LibraryRep::from(mink), bis.into()] {
            let mut max_dummy = 0;

            atom = atom.replace_map(|term, _, out| {
                let AtomView::Fun(f) = term else {
                    return;
                };

                let name = f.get_symbol();
                if name != rep.symbol() {
                    return;
                }

                let mut fbuilder = FunctionBuilder::new(name);

                if f.get_nargs() == 1 {
                    fbuilder = fbuilder.add_arg(4);
                    let arg = f.iter().next().unwrap();
                    let a = if let Ok(a) = i64::try_from(arg) {
                        if a < 0 {
                            let a = (-a) as usize;

                            a
                        } else {
                            return;
                        }
                    } else if let AtomView::Fun(f) = arg
                        && f.get_symbol() == self.dummy
                        && f.get_nargs() == 1
                        && let Ok(a) = usize::try_from(f.iter().next().unwrap())
                    {
                        a
                    } else {
                        return;
                    };
                    if a > max_dummy {
                        max_dummy = a;
                    }
                    fbuilder = fbuilder.add_arg(Atom::from(dummy(a)));
                    **out = fbuilder.finish();
                }
            });
        }

        let mink = Minkowski {}.new_rep(4);
        let bis = Bispinor {}.new_rep(4);
        // debug!(after_dummies = atom.printer(LOGPRINTOPTS).to_string());
        // debug!("After dummies{}", atom);

        let reps: Vec<_> = [
            (
                self.identity.f(&[bis.pattern(W_.a_), bis.pattern(W_.b_)]),
                bis.g(W_.a_, W_.b_),
            ),
            (
                self.identityl
                    .f(&[mink.pattern(W_.a_), mink.pattern(W_.b_)]),
                mink.g(W_.a_, W_.b_),
            ),
            (
                self.metric.f(&[mink.pattern(W_.a_), mink.pattern(W_.b_)]),
                mink.g(W_.a_, W_.b_),
            ),
            (
                self.gamma
                    .f(&[mink.pattern(W_.i_), bis.pattern(W_.a_), bis.pattern(W_.b_)]),
                AGS.gamma
                    .f(&[bis.pattern(W_.a_), bis.pattern(W_.b_), mink.pattern(W_.i_)]),
            ),
            (
                self.gamma5.f(&[bis.pattern(W_.a_), bis.pattern(W_.b_)]),
                AGS.gamma5.f(&[bis.pattern(W_.a_), bis.pattern(W_.b_)]),
            ),
            (
                self.projm.f(&[bis.pattern(W_.a_), bis.pattern(W_.b_)]),
                AGS.projm.f(&[bis.pattern(W_.a_), bis.pattern(W_.b_)]),
            ),
            (
                self.projp.f(&[bis.pattern(W_.a_), bis.pattern(W_.b_)]),
                AGS.projp.f(&[bis.pattern(W_.a_), bis.pattern(W_.b_)]),
            ),
            (
                self.sigma.f(&[
                    bis.pattern(W_.a_),
                    bis.pattern(W_.b_),
                    mink.pattern(W_.i_),
                    mink.pattern(W_.j_),
                ]),
                AGS.sigma.f(&[
                    bis.pattern(W_.a_),
                    bis.pattern(W_.b_),
                    mink.pattern(W_.i_),
                    mink.pattern(W_.j_),
                ]),
            ),
        ]
        .into_iter()
        .map(|(pat, rep)| Replacement::new(pat.to_pattern(), rep))
        .collect();

        // println!("out:{atom}");
        atom = atom.replace_multiple(&reps);

        atom = atom.replace_map(|term, _, out| {
            if let AtomView::Fun(f) = term {
                if f.get_symbol() == self.pslash {
                    let mut gamma = FunctionBuilder::new(AGS.gamma);

                    let mut count = 0;

                    for a in f.iter() {
                        count += 1;
                        if count <= 2 {
                            gamma = gamma.add_arg(a);
                        } else {
                            if let Ok(i) = i64::try_from(a) {
                                max_dummy += 1;

                                let minki: Slot<Minkowski, Aind> = mink.slot(dummy(max_dummy));

                                gamma = gamma.add_arg(minki.to_atom());

                                **out = gamma.finish()
                                    * GS.emr_mom(momenta[i as usize].1, minki.to_atom());
                                return;
                            }
                        }
                    }

                    if count == 2 {
                        max_dummy += 1;

                        let minki: Slot<Minkowski, Aind> = mink.slot(dummy(max_dummy));

                        gamma = gamma.add_arg(minki.to_atom());

                        **out = gamma.finish() * GS.emr_mom(momenta[0].1, minki.to_atom());
                    }
                }
            }
        });

        for (i, (f, e)) in momenta.iter().enumerate() {
            atom = atom
                .replace(function!(UFO.momentum, (i + 1) as i64, W_.a_))
                .with(match f {
                    Flow::Sink => GS.emr_mom(*e, W_.a_),
                    Flow::Source => -GS.emr_mom(*e, W_.a_),
                })
        }

        atom = atom
            .replace(function!(UFO.momentum, W_.a_))
            .with(GS.emr_mom(momenta[0].1, W_.a_));

        // debug!(out = atom.printer(LOGPRINTOPTS).to_string());
        // println!("out:{atom:#}");
        atom
    }

    pub(crate) fn reindex_color(
        &self,
        slots: &[&OrderedStructure<LibraryRep, Aind>],
        mut atom: Atom,
        dummy: impl Fn(usize) -> Aind,
    ) -> Atom {
        let adj = ColorAdjoint {};
        let fund = ColorFundamental {};
        let antifund = ColorAntiFundamental {};
        let sext = ColorSextet {};
        let antisext = ColorAntiSextet {};
        // for s in slots {
        //     debug!("{}", s);
        // }

        // debug!("In:{}", atom);
        // Fill in representations
        let reps: Vec<_> = [
            (
                self.t.f(&[W_.a_, W_.i_, W_.j_]),
                self.t.f(&[
                    adj.to_symbolic([W_.a_]),
                    fund.to_symbolic([W_.i_]),
                    antifund.to_symbolic([W_.j_]),
                ]),
            ),
            (
                self.f.f(&[W_.a_, W_.b_, W_.c_]),
                self.f.f(&[
                    adj.to_symbolic([W_.a_]),
                    adj.to_symbolic([W_.b_]),
                    adj.to_symbolic([W_.c_]),
                ]),
            ),
            (
                self.d.f(&[W_.a_, W_.b_, W_.c_]),
                self.d.f(&[
                    adj.to_symbolic([W_.a_]),
                    adj.to_symbolic([W_.b_]),
                    adj.to_symbolic([W_.c_]),
                ]),
            ),
            (
                self.levicivita.f(&[W_.i_, W_.j_, W_.k_]),
                self.levicivita.f(&[
                    fund.to_symbolic([W_.i_]),
                    fund.to_symbolic([W_.j_]),
                    fund.to_symbolic([W_.k_]),
                ]),
            ),
            (
                self.antilevicivita.f(&[W_.i_, W_.j_, W_.k_]),
                self.levicivita.f(&[
                    antifund.to_symbolic([W_.i_]),
                    antifund.to_symbolic([W_.j_]),
                    antifund.to_symbolic([W_.k_]),
                ]),
            ),
            (
                self.t6.f(&[W_.a_, W_.i_, W_.j_]),
                self.t6.f(&[
                    adj.to_symbolic([W_.a_]),
                    sext.to_symbolic([W_.i_]),
                    antisext.to_symbolic([W_.j_]),
                ]),
            ),
            (
                self.k6.f(&[W_.a_, W_.i_, W_.j_]),
                self.k6.f(&[
                    sext.to_symbolic([W_.a_]),
                    antifund.to_symbolic([W_.i_]),
                    antifund.to_symbolic([W_.j_]),
                ]),
            ),
            (
                self.k6bar.f(&[W_.a_, W_.i_, W_.j_]),
                self.k6bar.f(&[
                    antisext.to_symbolic([W_.a_]),
                    fund.to_symbolic([W_.i_]),
                    fund.to_symbolic([W_.j_]),
                ]),
            ),
        ]
        .into_iter()
        .map(|(pat, rep)| Replacement::new(pat.to_pattern(), rep))
        .collect();

        // println!("preid:{atom}");
        atom = atom.replace_multiple(&reps);

        // Now identify the kroneker reps:

        atom = atom
            .replace(
                function!(self.identity, W_.a_, W_.b_)
                    * function!(W_.f_, W_.a___, function!(W_.g_, W_.a_), W_.b___),
            )
            .repeat()
            .with(
                self.identity.f((&[W_.g_.f(&[W_.a_])], &[W_.b_]))
                    * W_.f_.f((&[W_.a___], &[W_.g_.f(&[W_.a_])], &[W_.b___])),
            );

        // debug!("Postid:{}", atom);

        // println!("postid:{atom}");

        for (i, s) in slots.iter().enumerate() {
            let i = i + 1;
            for (shift, r) in s.external_structure_iter().enumerate() {
                let rep = r.rep_name();
                let wrappedi = self.idx(i, shift + 1);
                let wrappedi = wrappedi.as_view();

                // println!("{}", rep.to_symbolic([i]));
                atom = atom
                    .replace(rep.to_symbolic([wrappedi]))
                    .level_range((1, Some(1)))
                    .with(r.to_atom())
                    .replace(self.identity.f((&[wrappedi], &[W_.i_])))
                    .with(self.identity.f((&[r.to_atom()], &[W_.i_])))
                    .replace(self.identity.f((&[W_.i_], &[wrappedi])))
                    .with(self.identity.f((&[W_.i_], &[r.to_atom()])))
                    .replace(rep.to_symbolic([i]))
                    .level_range((1, Some(1)))
                    .with(r.to_atom())
                    .replace(self.identity.f((&[i], &[W_.i_])))
                    .with(self.identity.f((&[r.to_atom()], &[W_.i_])))
                    .replace(self.identity.f((&[W_.i_], &[i])))
                    .with(self.identity.f((&[W_.i_], &[r.to_atom()])));
            }
        }

        // println!("slot:{atom}");
        // Handle dummies:
        //
        for rep in [
            LibraryRep::from(adj).new_rep(8),
            fund.new_rep(3).to_lib(),
            antifund.new_rep(3).to_lib(),
            sext.new_rep(6).to_lib(),
            antisext.new_rep(6).to_lib(),
        ] {
            let mut max_dummy = 0;

            atom = atom.replace_map(|term, _, out| {
                let AtomView::Fun(f) = term else {
                    return;
                };

                let name = f.get_symbol();
                if name != rep.rep.symbol() {
                    return;
                }

                let mut fbuilder = FunctionBuilder::new(name);

                if f.get_nargs() == 1 {
                    fbuilder = fbuilder.add_arg(rep.dim.to_symbolic());
                    let arg = f.iter().next().unwrap();
                    let a = if let Ok(a) = i64::try_from(arg) {
                        if a < 0 {
                            let a = (-a) as usize;

                            a
                        } else {
                            return;
                        }
                    } else if let AtomView::Fun(f) = arg
                        && f.get_symbol() == self.dummy
                        && f.get_nargs() == 1
                        && let Ok(a) = usize::try_from(f.iter().next().unwrap())
                    {
                        a
                    } else {
                        return;
                    };
                    if a > max_dummy {
                        max_dummy = a;
                    }
                    fbuilder = fbuilder.add_arg(Atom::from(dummy(a)));
                    **out = fbuilder.finish();
                }
            });
        }

        let reps: Vec<_> = [
            (
                self.identity.f(&[W_.a_, W_.b_]),
                ETS.metric.f(&[W_.a_, W_.b_]),
            ),
            (
                self.t.f(&[
                    adj.to_symbolic([W_.a__]),
                    fund.to_symbolic([W_.i__]),
                    antifund.to_symbolic([W_.j__]),
                ]),
                CS.t.f(&[
                    adj.to_symbolic([W_.a__]),
                    fund.to_symbolic([W_.i__]),
                    antifund.to_symbolic([W_.j__]),
                ]),
            ),
            (
                self.f.f(&[
                    adj.to_symbolic([W_.a__]),
                    adj.to_symbolic([W_.b__]),
                    adj.to_symbolic([W_.c__]),
                ]),
                CS.f.f(&[
                    adj.to_symbolic([W_.a__]),
                    adj.to_symbolic([W_.b__]),
                    adj.to_symbolic([W_.c__]),
                ]),
            ),
        ]
        .into_iter()
        .map(|(pat, rep)| {
            let a = Replacement::new(pat.to_pattern(), rep);
            a
        })
        .collect();

        atom = atom.replace_multiple(&reps);

        // println!("Out:{:#}", atom);
        atom
    }
}

#[cfg(test)]
pub mod test {

    #[test]
    fn ufo_spin_processing() {}
}
