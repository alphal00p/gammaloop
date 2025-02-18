use std::sync::LazyLock;

use spenso::shadowing::ETS;
use spenso::structure::representation::{
    BaseRepName, Bispinor, ColorAdjoint, ColorFundamental, ColorSextet, Minkowski,
};
use symbolica::atom::AtomCore;
use symbolica::atom::FunctionAttribute;
use symbolica::id::MatchSettings;
use symbolica::{
    atom::{Atom, Symbol},
    function,
    id::Replacement,
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
    metric: Symbol::new_with_attributes("Metric", &[FunctionAttribute::Symmetric]).unwrap(),
    momentum: symb!("P"),
    levicivita: symb!("Epsilon"),
    t: symb!("T"),
    f: Symbol::new_with_attributes("f", &[FunctionAttribute::Antisymmetric]).unwrap(),
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
            function!(UFO.t, a_, b_, c_),
            function!(
                UFO.t,
                coad.to_pattern_wrapped(a_),
                cof.to_pattern_wrapped(b_),
                coaf.to_pattern_wrapped(c_)
            ),
        ),
        (
            function!(UFO.f, a_, b_, c_),
            function!(
                UFO.f,
                coad.to_pattern_wrapped(a_),
                coad.to_pattern_wrapped(b_),
                coad.to_pattern_wrapped(c_)
            ),
        ),
        (
            function!(UFO.d, a_, b_, c_),
            function!(
                UFO.d,
                coad.to_pattern_wrapped(a_),
                coad.to_pattern_wrapped(b_),
                coad.to_pattern_wrapped(c_)
            ),
        ),
        (
            function!(UFO.levicivita, a_, b_, c_),
            function!(
                UFO.levicivita,
                cof.to_pattern_wrapped(a_),
                cof.to_pattern_wrapped(b_),
                cof.to_pattern_wrapped(c_)
            ),
        ),
        (
            function!(UFO.antilevicivita, a_, b_, c_),
            function!(
                UFO.levicivita,
                coaf.to_pattern_wrapped(a_),
                coaf.to_pattern_wrapped(b_),
                coaf.to_pattern_wrapped(c_)
            ),
        ),
        (
            function!(UFO.t6, a_, b_, c_),
            function!(
                UFO.t6,
                coad.to_pattern_wrapped(a_),
                cos.to_pattern_wrapped(b_),
                coas.to_pattern_wrapped(c_)
            ),
        ),
        (
            function!(UFO.k6, a_, b_, c_),
            function!(
                UFO.k6,
                coaf.to_pattern_wrapped(a_),
                coaf.to_pattern_wrapped(b_),
                cos.to_pattern_wrapped(c_)
            ),
        ),
        (
            function!(UFO.k6bar, a_, b_, c_),
            function!(
                UFO.k6bar,
                coas.to_pattern_wrapped(a_),
                cof.to_pattern_wrapped(b_),
                cof.to_pattern_wrapped(c_)
            ),
        ),
    ];

    let reps: Vec<Replacement> = reps
        .into_iter()
        .map(|(a, b)| Replacement::new(a.to_pattern(), b.to_pattern()))
        .collect();

    atom.replace_all_multiple(&reps)
}
pub fn preprocess_ufo_spin_wrapped(atom: Atom) -> Atom {
    let a_ = symb!("a_");
    let b_ = symb!("b_");
    let c_ = symb!("c_");
    let d_ = symb!("d_");
    let wa_ = function!(symb!("indexid"), Atom::new_var(symb!("a_")));
    let wb_ = function!(symb!("indexid"), Atom::new_var(symb!("b_")));
    let wc_ = function!(symb!("indexid"), Atom::new_var(symb!("c_")));
    let wd_ = function!(symb!("indexid"), Atom::new_var(symb!("d_")));

    let bis = Bispinor::rep(4);
    let mink = Minkowski::rep(4);

    let reps = vec![
        (
            function!(UFO.identity, a_, b_),
            function!(ETS.id, bis.pattern(&wa_), bis.pattern(&wb_)),
        ),
        (
            function!(UFO.identityl, a_, b_),
            function!(ETS.id, mink.pattern(&wa_), mink.pattern(&wb_)),
        ),
        (
            function!(UFO.gamma, a_, b_, c_),
            function!(
                ETS.gamma,
                mink.pattern(&wa_),
                bis.pattern(&wb_),
                bis.pattern(&wc_)
            ),
        ),
        (
            function!(UFO.gamma5, a_, b_),
            function!(ETS.gamma5, bis.pattern(&wa_), bis.pattern(&wb_)),
        ),
        (
            function!(UFO.projm, a_, b_),
            function!(ETS.proj_m, bis.pattern(&wa_), bis.pattern(&wb_)),
        ),
        (
            function!(UFO.projp, a_, b_),
            function!(ETS.proj_p, bis.pattern(&wa_), bis.pattern(&wb_)),
        ),
        (
            function!(UFO.sigma, a_, b_, c_, d_),
            function!(
                ETS.sigma,
                mink.pattern(&wa_),
                mink.pattern(&wb_),
                bis.pattern(&wc_),
                bis.pattern(&wd_)
            ),
        ),
        (
            function!(UFO.charge_conj, a_, b_),
            function!(UFO.charge_conj, bis.pattern(&wa_), bis.pattern(&wb_)),
        ),
        (
            function!(UFO.metric, a_, b_),
            function!(ETS.metric, mink.pattern(&wa_), mink.pattern(&wb_)),
        ),
    ];

    let settings = MatchSettings {
        rhs_cache_size: 0,
        ..Default::default()
    };
    let reps: Vec<_> = reps
        .into_iter()
        .map(|(a, b)| {
            Replacement::new(a.to_pattern(), b.to_pattern()).with_settings(settings.clone())
        })
        .collect();

    atom.replace_all_multiple(&reps)
}

#[cfg(test)]
pub mod test {

    #[test]
    fn ufo_spin_processing() {}
}
