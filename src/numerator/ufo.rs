use std::sync::LazyLock;

use idenso::color::CS;
use idenso::gamma::AGS;
use idenso::representations::{
    Bispinor, ColorAdjoint, ColorAntiFundamental, ColorAntiSextet, ColorFundamental, ColorSextet,
};
use linnet::half_edge::involution::EdgeIndex;
use spenso::network::library::symbolic::ETS;
use spenso::structure::representation::{BaseRepName, LibraryRep, Minkowski, RepName};
use spenso::structure::slot::{IsAbstractSlot, Slot};
use spenso::structure::{OrderedStructure, TensorStructure};
use symbolica::atom::{AtomCore, AtomView};
use symbolica::atom::FunctionBuilder;
use symbolica::domains::float::NumericalFloatLike;
use symbolica::domains::rational::Rational;
use symbolica::symbol;
use symbolica::{
    atom::{Atom, Symbol},
    function,
    id::Replacement,
};

use crate::symbolica_ext::CallSymbol;
use crate::utils::{GS, W_};

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
    pub k6bar: Symbol,
    pub pslash: Symbol,
    pub complex: Symbol,
}

#[allow(dead_code)]
pub static UFO: LazyLock<UFOSymbols> = LazyLock::new(|| UFOSymbols {
    identity: symbol!("Identity"),
    identityl: symbol!("IdentityL"),
    gamma: symbol!("Gamma"),
    gamma5: symbol!("Gamma5"),
    projm: symbol!("ProjM"),
    projp: symbol!("ProjP"),
    sigma: symbol!("Sigma"),
    charge_conj: symbol!("C"),
    metric: symbol!("Metric";Symmetric),
    momentum: symbol!("P"),
    levicivita: symbol!("Epsilon"),
    t: symbol!("T"),
    f: symbol!("f";Antisymmetric),
    d: symbol!("d"),
    antilevicivita: symbol!("EpsilonBar"),
    t6: symbol!("T6"),
    k6: symbol!("K6"),
    k6bar: symbol!("K6Bar"),
    pslash: symbol!("PSlash"),
    complex: symbol!("complex"),
});

impl UFOSymbols {
    pub fn normalize_complex(&self, a: impl AtomCore) -> Atom {
        a.replace_map(|term, _, out| {
            if let AtomView::Fun(f) = term {
                if f.get_symbol() == self.complex {
                    let mut re = Rational::zero();
                    let mut im = Rational::zero();
                    let mut count = 0;
                    for i in f.iter() {
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

                    // println!("{re}+i{im}");

                    if count == 2 {
                        *out = Atom::num(symbolica::domains::float::Complex { re, im });
                        return true;
                    }
                }
            }
            false
        })
    }

    pub fn reindex_spin(
        &self,
        slots: &[&OrderedStructure<LibraryRep, Aind>],
        momenta: &[EdgeIndex],
        mut atom: Atom,
        dummy: impl Fn(usize) -> Aind,
    ) -> Atom {
        let mink: LibraryRep = Minkowski {}.into();
        let bis: LibraryRep = Bispinor {}.into();

        // for s in slots {
        //     println!("{}", s);
        // }
        // println!("in:{atom}");
        //
        atom = self.normalize_complex(atom);
        for (i, s) in slots.iter().enumerate() {
            let i = (i + 1) as i64;
            for r in s.external_structure_iter() {
                let rep = r.rep_name();
                match rep {
                    _ if rep == mink => {
                        atom = atom
                            // Fill in gamma
                            .replace(self.gamma.f((&[i], &[W_.x_, W_.y_])))
                            .with(self.gamma.f((&[r.to_atom()], &[W_.x_, W_.y_])))
                            // Fill in kroneker
                            .replace(self.identityl.f((&[i], &[W_.a_])))
                            .with(self.identityl.f((&[r.to_atom()], &[W_.a_])))
                            .replace(self.identityl.f((&[W_.a_], &[i])))
                            .with(self.identityl.f((&[W_.a_], &[r.to_atom()])))
                            // Fill in metric
                            .replace(self.metric.f((&[i], &[W_.a_])))
                            .with(self.metric.f((&[r.to_atom()], &[W_.a_])))
                            .replace(self.metric.f((&[W_.a_], &[i])))
                            .with(self.metric.f((&[W_.a_], &[r.to_atom()])))
                            // Fill in momentum
                            .replace(self.momentum.f((&[i], &[W_.i_])))
                            .with(self.momentum.f((&[r.to_atom()], &[W_.i_])))
                            // Fill in sigma
                            .replace(self.sigma.f((&[i], &[W_.j_, W_.a_, W_.b_])))
                            .with(self.sigma.f((&[r.to_atom()], &[W_.j_, W_.a_, W_.b_])))
                            .replace(self.sigma.f((&[W_.j_], &[i], &[W_.a_, W_.b_])))
                            .with(self.sigma.f((&[W_.j_], &[r.to_atom()], &[W_.a_, W_.b_])))
                            // Fill in levi-civita
                            .replace(self.levicivita.f((&[W_.a___], &[i], &[W_.b___])))
                            .with(self.levicivita.f((&[W_.a___], &[r.to_atom()], &[W_.b___])))
                    }
                    _ if rep == bis => {
                        let reps: Vec<_> = [
                            // Fill in pslash
                            (
                                self.pslash.f((&[i], &[W_.i_, W_.a___])),
                                self.pslash.f((&[r.to_atom()], &[W_.i_, W_.a___])),
                            ),
                            (
                                self.pslash.f((&[W_.i_], &[i], &[W_.a___])),
                                self.pslash.f((&[W_.i_], &[r.to_atom()], &[W_.i_, W_.a___])),
                            ),
                            (
                                self.identity.f((&[i], &[W_.a_])),
                                self.identity.f((&[r.to_atom()], &[W_.a_])),
                            ),
                            (
                                self.identity.f((&[W_.a_], &[i])),
                                self.identity.f((&[W_.a_], &[r.to_atom()])),
                            ),
                            (
                                self.gamma.f((&[W_.i_], &[i], &[W_.a_])),
                                self.gamma.f((&[W_.i_], &[r.to_atom()], &[W_.a_])),
                            ),
                            (
                                self.gamma.f((&[W_.i_, W_.a_], &[i])),
                                self.gamma.f((&[W_.i_, W_.a_], &[r.to_atom()])),
                            ),
                            (
                                self.gamma5.f((&[i], &[W_.a_])),
                                self.gamma5.f((&[r.to_atom()], &[W_.a_])),
                            ),
                            (
                                self.gamma5.f((&[W_.a_], &[i])),
                                self.gamma5.f((&[W_.a_], &[r.to_atom()])),
                            ),
                            (
                                self.projm.f((&[i], &[W_.a_])),
                                self.projm.f((&[r.to_atom()], &[W_.a_])),
                            ),
                            (
                                self.projm.f((&[W_.a_], &[i])),
                                self.projm.f((&[W_.a_], &[r.to_atom()])),
                            ),
                            (
                                self.projp.f((&[i], &[W_.a_])),
                                self.projp.f((&[r.to_atom()], &[W_.a_])),
                            ),
                            (
                                self.projp.f((&[W_.a_], &[i])),
                                self.projp.f((&[W_.a_], &[r.to_atom()])),
                            ),
                            (
                                self.charge_conj.f((&[i], &[W_.a_])),
                                self.charge_conj.f((&[r.to_atom()], &[W_.a_])),
                            ),
                            (
                                self.charge_conj.f((&[W_.a_], &[i])),
                                self.charge_conj.f((&[W_.a_], &[r.to_atom()])),
                            ),
                            (
                                self.sigma.f((&[W_.i_, W_.j_], &[i], &[W_.a_])),
                                self.sigma.f((&[W_.i_, W_.j_], &[r.to_atom()], &[W_.a_])),
                            ),
                            (
                                self.sigma.f((&[W_.i_, W_.j_, W_.a_], &[i])),
                                self.sigma.f((&[W_.i_, W_.j_, W_.a_], &[r.to_atom()])),
                            ),
                        ]
                        .into_iter()
                        .map(|(pat, rep)| Replacement::new(pat.to_pattern(), rep))
                        .collect();

                        // for r in &reps {
                        //     println!("{:#}", r);
                        // }

                        atom = atom.replace_multiple(&reps);
                    }
                    _ => {}
                }
            }
        }

        // Handle dummies:

        let mink = Minkowski {}.new_rep(4);
        let bis = Bispinor {}.new_rep(4);

        let mut max_dummy = 0;

        atom = atom.replace_map(|term, ctx, out| {
            if let AtomView::Fun(f) = term {
                let name = f.get_symbol();
                if name == self.gamma {
                    let mut fbuilder = FunctionBuilder::new(self.gamma);
                    let mut first = true;
                    for a in f.iter() {
                        if let Ok(a) = i64::try_from(a) {
                            if a < 0 {
                                let a = (-a) as usize;

                                if first {
                                    if a > max_dummy {
                                        max_dummy = a;
                                    }
                                    first = false;
                                    fbuilder =
                                        fbuilder.add_arg(mink.slot::<Aind, _>(dummy(a)).to_atom());
                                } else {
                                    fbuilder =
                                        fbuilder.add_arg(bis.slot::<Aind, _>(dummy(a)).to_atom());
                                }
                            } else {
                                fbuilder = fbuilder.add_arg(a);
                            }
                        } else {
                            fbuilder = fbuilder.add_arg(a);
                        }
                    }

                    *out = fbuilder.finish();
                    true
                } else if name == self.charge_conj
                    || name == self.gamma5
                    || name == self.identity
                    || name == self.pslash
                    || name == self.projm
                    || name == self.projp
                {
                    let mut fbuilder = FunctionBuilder::new(name);
                    let mut count = 0;
                    for a in f.iter() {
                        count += 1;
                        if count <= 2 {
                            if let Ok(a) = i64::try_from(a) {
                                if a < 0 {
                                    let a = (-a) as usize;
                                    fbuilder =
                                        fbuilder.add_arg(bis.slot::<Aind, _>(dummy(a)).to_atom());
                                } else {
                                    fbuilder = fbuilder.add_arg(a);
                                }
                            } else {
                                fbuilder = fbuilder.add_arg(a);
                            }
                        }
                    }
                    *out = fbuilder.finish();
                    true
                } else if name == self.identityl
                    || name == self.levicivita
                    || name == self.momentum
                    || name == self.metric
                {
                    let mut fbuilder = FunctionBuilder::new(name);
                    for a in f.iter() {
                        if let Ok(a) = i64::try_from(a) {
                            if a < 0 {
                                let a = (-a) as usize;
                                if a > max_dummy {
                                    max_dummy = a;
                                }
                                fbuilder =
                                    fbuilder.add_arg(mink.slot::<Aind, _>(dummy(a)).to_atom());
                            } else {
                                fbuilder = fbuilder.add_arg(a);
                            }
                        } else {
                            fbuilder = fbuilder.add_arg(a);
                        }
                    }
                    *out = fbuilder.finish();
                    true
                } else if name == self.sigma {
                    let mut fbuilder = FunctionBuilder::new(self.gamma);
                    let mut count = 0;
                    for a in f.iter() {
                        if let Ok(a) = i64::try_from(a) {
                            if a < 0 {
                                let a = (-a) as usize;

                                if count < 2 {
                                    count += 1;
                                    if a > max_dummy {
                                        max_dummy = a;
                                    }
                                    fbuilder =
                                        fbuilder.add_arg(mink.slot::<Aind, _>(dummy(a)).to_atom());
                                } else {
                                    fbuilder =
                                        fbuilder.add_arg(bis.slot::<Aind, _>(dummy(a)).to_atom());
                                }
                            } else {
                                fbuilder = fbuilder.add_arg(a);
                            }
                        } else {
                            fbuilder = fbuilder.add_arg(a);
                        }
                    }

                    *out = fbuilder.finish();
                    true
                } else {
                    false
                }
            } else {
                false
            }
        });

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

        atom = atom.replace_multiple(&reps);

        atom = atom.replace_map(|term, ctx, out| {
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

                                *out = gamma.finish()
                                    * GS.emr_mom(momenta[i as usize], minki.to_atom());

                                return true;
                            }
                        }
                    }

                    if count == 2 {
                        max_dummy += 1;

                        let minki: Slot<Minkowski, Aind> = mink.slot(dummy(max_dummy));

                        gamma = gamma.add_arg(minki.to_atom());

                        *out = gamma.finish() * GS.emr_mom(momenta[0], minki.to_atom());

                        return true;
                    }
                }
            }
            false
        });

        for (i, e) in momenta.iter().enumerate() {
            atom = atom
                .replace(function!(UFO.momentum, i as i64, W_.a_))
                .with(GS.emr_mom(*e, W_.a_));
        }

        // println!("out:{atom}");
        atom
    }

    pub fn reindex_color(
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

        atom = atom.replace_multiple(&reps);

        // Now identify the kroneker reps:

        // println!("preid:{atom}");
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

        // println!("postid:{atom}");
        atom = self.normalize_complex(atom);
        for (i, s) in slots.iter().enumerate() {
            let i = (i + 1) as i64;
            for r in s.external_structure_iter() {
                let rep = r.rep_name();

                // println!("{r}");
                atom = atom
                    .replace(rep.to_symbolic([i]))
                    .with(r.to_atom())
                    .replace(self.identity.f((&[i], &[W_.i_])))
                    .with(self.identity.f((&[r.to_atom()], &[W_.i_])))
                    .replace(self.identity.f((&[W_.i_], &[i])))
                    .with(self.identity.f((&[W_.i_], &[r.to_atom()])));
            }
        }

        // println!("slot:{atom}");
        // Handle dummies:

        for rep in [
            LibraryRep::from(adj),
            fund.into(),
            antifund.into(),
            sext.into(),
            antisext.into(),
        ] {
            let mut max_dummy = 0;

            if let AtomView::Fun(r) = rep.to_symbolic::<Atom>([]).as_view() {
                atom = atom.replace_map(|term, ctx, out| {
                    if let AtomView::Fun(f) = term {
                        let name = f.get_symbol();
                        if name == r.get_symbol() {
                            let mut fbuilder = FunctionBuilder::new(name);
                            if f.get_nargs() == 1 {
                                let arg = f.iter().next().unwrap();
                                if let Ok(a) = i64::try_from(arg) {
                                    if a < 0 {
                                        let a = (-a) as usize;

                                        if a > max_dummy {
                                            max_dummy = a;
                                        }
                                        fbuilder = fbuilder
                                            .add_arg(rep.to_symbolic([Atom::from(dummy(a))]));
                                        *out = fbuilder.finish();
                                        return true;
                                    }
                                }
                            }
                        }
                    }
                    false
                });
            }
        }

        let reps: Vec<_> = [
            (
                self.identity.f(&[W_.a_, W_.b_]),
                ETS.metric.f(&[W_.a_, W_.b_]),
            ),
            (
                self.t.f(&[
                    adj.to_symbolic([W_.a_]),
                    fund.to_symbolic([W_.i_]),
                    antifund.to_symbolic([W_.j_]),
                ]),
                CS.t.f(&[
                    adj.to_symbolic([CS.nc * CS.nc - 1, Atom::var(W_.a_)]),
                    fund.to_symbolic([CS.nc, W_.i_]),
                    antifund.to_symbolic([CS.nc, W_.j_]),
                ]),
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
                    adj.to_symbolic([W_.a_]),
                    adj.to_symbolic([W_.b_]),
                    adj.to_symbolic([W_.c_]),
                ]),
                CS.f.f(&[
                    adj.to_symbolic([CS.nc * CS.nc - 1, Atom::var(W_.a_)]),
                    adj.to_symbolic([CS.nc * CS.nc - 1, Atom::var(W_.b_)]),
                    adj.to_symbolic([CS.nc * CS.nc - 1, Atom::var(W_.c_)]),
                ]),
            ),
        ]
        .into_iter()
        .map(|(pat, rep)| Replacement::new(pat.to_pattern(), rep))
        .collect();

        atom = atom.replace_multiple(&reps);

        atom
    }
}

pub fn preprocess_ufo_color_wrapped(atom: Atom) -> Atom {
    let a_ = symbol!("a_");
    let b_ = symbol!("b_");
    let c_ = symbol!("c_");

    // let coad = ColorAdjoint::rep(8);
    // let cof = ColorFundamental::rep(3);
    // let coaf = ColorFundamental::rep(3).dual();
    // let cos = ColorSextet::rep(6);
    // let coas = ColorSextet::rep(6).dual();

    // let reps = vec![
    //     (
    //         function!(UFO.t, a_, b_, c_),
    //         function!(
    //             UFO.t,
    //             coad.to_pattern_wrapped(a_),
    //             cof.to_pattern_wrapped(b_),
    //             coaf.to_pattern_wrapped(c_)
    //         ),
    //     ),
    //     (
    //         function!(UFO.f, a_, b_, c_),
    //         function!(
    //             UFO.f,
    //             coad.to_pattern_wrapped(a_),
    //             coad.to_pattern_wrapped(b_),
    //             coad.to_pattern_wrapped(c_)
    //         ),
    //     ),
    //     (
    //         function!(UFO.d, a_, b_, c_),
    //         function!(
    //             UFO.d,
    //             coad.to_pattern_wrapped(a_),
    //             coad.to_pattern_wrapped(b_),
    //             coad.to_pattern_wrapped(c_)
    //         ),
    //     ),
    //     (
    //         function!(UFO.levicivita, a_, b_, c_),
    //         function!(
    //             UFO.levicivita,
    //             cof.to_pattern_wrapped(a_),
    //             cof.to_pattern_wrapped(b_),
    //             cof.to_pattern_wrapped(c_)
    //         ),
    //     ),
    //     (
    //         function!(UFO.antilevicivita, a_, b_, c_),
    //         function!(
    //             UFO.levicivita,
    //             coaf.to_pattern_wrapped(a_),
    //             coaf.to_pattern_wrapped(b_),
    //             coaf.to_pattern_wrapped(c_)
    //         ),
    //     ),
    //     (
    //         function!(UFO.t6, a_, b_, c_),
    //         function!(
    //             UFO.t6,
    //             coad.to_pattern_wrapped(a_),
    //             cos.to_pattern_wrapped(b_),
    //             coas.to_pattern_wrapped(c_)
    //         ),
    //     ),
    //     (
    //         function!(UFO.k6, a_, b_, c_),
    //         function!(
    //             UFO.k6,
    //             coaf.to_pattern_wrapped(a_),
    //             coaf.to_pattern_wrapped(b_),
    //             cos.to_pattern_wrapped(c_)
    //         ),
    //     ),
    //     (
    //         function!(UFO.k6bar, a_, b_, c_),
    //         function!(
    //             UFO.k6bar,
    //             coas.to_pattern_wrapped(a_),
    //             cof.to_pattern_wrapped(b_),
    //             cof.to_pattern_wrapped(c_)
    //         ),
    //     ),
    // ];

    // let reps: Vec<Replacement> = []
    //     .into_iter()
    //     .map(|(a, b)| Replacement::new(a.to_pattern(), b.to_pattern()))
    //     .collect();

    // atom.replace_multiple(&reps)
    todo!()
}
pub fn preprocess_ufo_spin_wrapped(atom: Atom) -> Atom {
    // let a_ = symbol!("a_");
    // let b_ = symbol!("b_");
    // let c_ = symbol!("c_");
    // let d_ = symbol!("d_");
    // let wa_ = function!(symbol!("indexid"), Atom::var(symbol!("a_")));
    // let wb_ = function!(symbol!("indexid"), Atom::var(symbol!("b_")));
    // let wc_ = function!(symbol!("indexid"), Atom::var(symbol!("c_")));
    // let wd_ = function!(symbol!("indexid"), Atom::var(symbol!("d_")));

    // let bis = Bispinor {}.new_rep(4);
    // let mink = Minkowski {}.new_rep(4);

    // let reps = vec![
    //     (
    //         function!(UFO.identity, a_, b_),
    //         function!(ETS.id, bis.pattern(&wa_), bis.pattern(&wb_)),
    //     ),
    //     (
    //         function!(UFO.identityl, a_, b_),
    //         function!(ETS.id, mink.pattern(&wa_), mink.pattern(&wb_)),
    //     ),
    //     (
    //         function!(UFO.gamma, a_, b_, c_),
    //         function!(
    //             ETS.gamma,
    //             mink.pattern(&wa_),
    //             bis.pattern(&wb_),
    //             bis.pattern(&wc_)
    //         ),
    //     ),
    //     (
    //         function!(UFO.gamma5, a_, b_),
    //         function!(ETS.gamma5, bis.pattern(&wa_), bis.pattern(&wb_)),
    //     ),
    //     (
    //         function!(UFO.projm, a_, b_),
    //         function!(ETS.proj_m, bis.pattern(&wa_), bis.pattern(&wb_)),
    //     ),
    //     (
    //         function!(UFO.projp, a_, b_),
    //         function!(ETS.proj_p, bis.pattern(&wa_), bis.pattern(&wb_)),
    //     ),
    //     (
    //         function!(UFO.sigma, a_, b_, c_, d_),
    //         function!(
    //             ETS.sigma,
    //             mink.pattern(&wa_),
    //             mink.pattern(&wb_),
    //             bis.pattern(&wc_),
    //             bis.pattern(&wd_)
    //         ),
    //     ),
    //     (
    //         function!(UFO.charge_conj, a_, b_),
    //         function!(UFO.charge_conj, bis.pattern(&wa_), bis.pattern(&wb_)),
    //     ),
    //     (
    //         function!(UFO.metric, a_, b_),
    //         function!(ETS.metric, mink.pattern(&wa_), mink.pattern(&wb_)),
    //     ),
    // ];

    // let settings = MatchSettings {
    //     rhs_cache_size: 0,
    //     ..Default::default()
    // };
    // let reps: Vec<_> = reps
    //     .into_iter()
    //     .map(|(a, b)| {
    //         Replacement::new(a.to_pattern(), b.to_pattern()).with_settings(settings.clone())
    //     })
    //     .collect();

    // atom.replace_multiple(&reps)
    //
    todo!()
}

#[cfg(test)]
pub mod test {

    #[test]
    fn ufo_spin_processing() {}
}
