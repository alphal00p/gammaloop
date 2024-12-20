use std::sync::LazyLock;

use spenso::shadowing::ETS;
use spenso::structure::representation::{
    BaseRepName, Bispinor, ColorAdjoint, ColorFundamental, ColorSextet, Minkowski,
};
use symbolica::id::MatchSettings;
use symbolica::state::FunctionAttribute;
use symbolica::{
    atom::{Atom, Symbol},
    fun,
    id::{Pattern, PatternOrMap, Replacement},
    state::State,
    symb,
};

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
}

#[allow(dead_code)]
pub static UFO: LazyLock<UFOSymbols> = LazyLock::new(|| UFOSymbols {
    identity: symb!("Identity"),
    identityl: symb!("IdentityL"),
    gamma: symb!("Gamma"),
    gamma5: symb!("Gamma5"),
    projm: symb!("ProjM"),
    projp: symb!("ProjP"),
    sigma: symb!("Sigma"),
    charge_conj: symb!("C"),
    metric: State::get_symbol_with_attributes("Metric", &[FunctionAttribute::Symmetric]).unwrap(),
    momentum: symb!("P"),
    levicivita: symb!("Epsilon"),
    t: symb!("T"),
    f: State::get_symbol_with_attributes("f", &[FunctionAttribute::Antisymmetric]).unwrap(),
    d: symb!("d"),
    antilevicivita: symb!("EpsilonBar"),
    t6: symb!("T6"),
    k6: symb!("K6"),
    k6bar: symb!("K6Bar"),
    pslash: symb!("PSlash"),
});

pub fn preprocess_ufo_color_wrapped(atom: Atom) -> Atom {
    let a_ = symb!("a_");
    let b_ = symb!("b_");
    let c_ = symb!("c_");

    let coad = ColorAdjoint::rep(8);
    let cof = ColorFundamental::rep(3);
    let coaf = ColorFundamental::rep(3).dual();
    let cos = ColorSextet::rep(6);
    let coas = ColorSextet::rep(6).dual();

    let reps = vec![
        (
            fun!(UFO.t, a_, b_, c_),
            fun!(
                UFO.t,
                coad.to_pattern_wrapped(a_),
                cof.to_pattern_wrapped(b_),
                coaf.to_pattern_wrapped(c_)
            ),
        ),
        (
            fun!(UFO.f, a_, b_, c_),
            fun!(
                UFO.f,
                coad.to_pattern_wrapped(a_),
                coad.to_pattern_wrapped(b_),
                coad.to_pattern_wrapped(c_)
            ),
        ),
        (
            fun!(UFO.d, a_, b_, c_),
            fun!(
                UFO.d,
                coad.to_pattern_wrapped(a_),
                coad.to_pattern_wrapped(b_),
                coad.to_pattern_wrapped(c_)
            ),
        ),
        (
            fun!(UFO.levicivita, a_, b_, c_),
            fun!(
                UFO.levicivita,
                cof.to_pattern_wrapped(a_),
                cof.to_pattern_wrapped(b_),
                cof.to_pattern_wrapped(c_)
            ),
        ),
        (
            fun!(UFO.antilevicivita, a_, b_, c_),
            fun!(
                UFO.levicivita,
                coaf.to_pattern_wrapped(a_),
                coaf.to_pattern_wrapped(b_),
                coaf.to_pattern_wrapped(c_)
            ),
        ),
        (
            fun!(UFO.t6, a_, b_, c_),
            fun!(
                UFO.t6,
                coad.to_pattern_wrapped(a_),
                cos.to_pattern_wrapped(b_),
                coas.to_pattern_wrapped(c_)
            ),
        ),
        (
            fun!(UFO.k6, a_, b_, c_),
            fun!(
                UFO.k6,
                coaf.to_pattern_wrapped(a_),
                coaf.to_pattern_wrapped(b_),
                cos.to_pattern_wrapped(c_)
            ),
        ),
        (
            fun!(UFO.k6bar, a_, b_, c_),
            fun!(
                UFO.k6bar,
                coas.to_pattern_wrapped(a_),
                cof.to_pattern_wrapped(b_),
                cof.to_pattern_wrapped(c_)
            ),
        ),
    ];

    let reps: Vec<(Pattern, PatternOrMap)> = reps
        .into_iter()
        .map(|(a, b)| (a.into_pattern(), b.into_pattern().into()))
        .collect();

    let reps: Vec<_> = reps.iter().map(|(a, b)| Replacement::new(a, b)).collect();

    atom.replace_all_multiple(&reps)
}
pub fn preprocess_ufo_spin_wrapped(atom: Atom) -> Atom {
    let a_ = symb!("a_");
    let b_ = symb!("b_");
    let c_ = symb!("c_");
    let d_ = symb!("d_");
    let wa_ = fun!(symb!("indexid"), Atom::new_var(symb!("a_")));
    let wb_ = fun!(symb!("indexid"), Atom::new_var(symb!("b_")));
    let wc_ = fun!(symb!("indexid"), Atom::new_var(symb!("c_")));
    let wd_ = fun!(symb!("indexid"), Atom::new_var(symb!("d_")));

    let bis = Bispinor::rep(4);
    let mink = Minkowski::rep(4);

    let reps = vec![
        (
            fun!(UFO.identity, a_, b_),
            fun!(ETS.id, bis.pattern(&wa_), bis.pattern(&wb_)),
        ),
        (
            fun!(UFO.identityl, a_, b_),
            fun!(ETS.id, mink.pattern(&wa_), mink.pattern(&wb_)),
        ),
        (
            fun!(UFO.gamma, a_, b_, c_),
            fun!(
                ETS.gamma,
                mink.pattern(&wa_),
                bis.pattern(&wb_),
                bis.pattern(&wc_)
            ),
        ),
        (
            fun!(UFO.gamma5, a_, b_),
            fun!(ETS.gamma5, bis.pattern(&wa_), bis.pattern(&wb_)),
        ),
        (
            fun!(UFO.projm, a_, b_),
            fun!(ETS.proj_m, bis.pattern(&wa_), bis.pattern(&wb_)),
        ),
        (
            fun!(UFO.projp, a_, b_),
            fun!(ETS.proj_p, bis.pattern(&wa_), bis.pattern(&wb_)),
        ),
        (
            fun!(UFO.sigma, a_, b_, c_, d_),
            fun!(
                ETS.sigma,
                mink.pattern(&wa_),
                mink.pattern(&wb_),
                bis.pattern(&wc_),
                bis.pattern(&wd_)
            ),
        ),
        (
            fun!(UFO.charge_conj, a_, b_),
            fun!(UFO.charge_conj, bis.pattern(&wa_), bis.pattern(&wb_)),
        ),
        (
            fun!(UFO.metric, a_, b_),
            fun!(ETS.metric, mink.pattern(&wa_), mink.pattern(&wb_)),
        ),
    ];

    let reps: Vec<(Pattern, PatternOrMap)> = reps
        .into_iter()
        .map(|(a, b)| (a.into_pattern(), b.into_pattern().into()))
        .collect();

    let settings = MatchSettings {
        rhs_cache_size: 0,
        ..Default::default()
    };
    let reps: Vec<_> = reps
        .iter()
        .map(|(a, b)| Replacement::new(a, b).with_settings(&settings))
        .collect();

    atom.replace_all_multiple(&reps)
}

#[cfg(test)]
pub mod test {

    #[test]
    fn ufo_spin_processing() {}
}
