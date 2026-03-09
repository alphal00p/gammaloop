use idenso::{
    gamma::AGS,
    metric::{MetricSimplifier, PermuteWithMetric},
    representations::{Bispinor, initialize},
};
use spenso::{
    network::{
        ExecutionResult, Sequential, SmallestDegree, TensorNetworkError, TensorOrScalarOrKey,
        library::{TensorLibraryData, symbolic::ExplicitKey},
        parsing::ParseSettings,
    },
    structure::{
        IndexlessNamedStructure, ScalarTensor,
        abstract_index::AbstractIndex,
        permuted::Perm,
        representation::{Minkowski, RepName},
    },
    tensors::{parametric::ParamTensor, symbolic::SymbolicTensor},
};
use spenso_hep_lib::{FUN_LIB, HEP_LIB, HepNet, HepTensor};
use symbolica::{
    atom::{Atom, AtomView, Symbol},
    function, symbol,
};

pub fn test_initialize() {
    let _a = HEP_LIB.get(&AGS.gamma_strct(4)).unwrap();

    initialize();
}

pub trait NetExt {
    fn execute_and_res(
        self,
    ) -> Result<HepTensor<AbstractIndex>, TensorNetworkError<ExplicitKey<AbstractIndex>, Symbol>>;
}

impl NetExt for HepNet<AbstractIndex> {
    fn execute_and_res(
        mut self,
    ) -> Result<HepTensor<AbstractIndex>, TensorNetworkError<ExplicitKey<AbstractIndex>, Symbol>>
    {
        self.execute::<Sequential, SmallestDegree, _, _, _>(&*HEP_LIB, &*FUN_LIB)?;
        match self.result()? {
            ExecutionResult::One => Ok(HepTensor::Param(ParamTensor::new_scalar(Atom::one()))),
            ExecutionResult::Zero => Ok(HepTensor::Param(ParamTensor::new_scalar(Atom::zero()))),
            ExecutionResult::Val(a) => match a {
                TensorOrScalarOrKey::Scalar(a) => {
                    Ok(HepTensor::Param(ParamTensor::new_scalar(a.clone())))
                }
                _ => match self.result_tensor(&*HEP_LIB)? {
                    ExecutionResult::One => {
                        Ok(HepTensor::Param(ParamTensor::new_scalar(Atom::one())))
                    }
                    ExecutionResult::Zero => {
                        Ok(HepTensor::Param(ParamTensor::new_scalar(Atom::zero())))
                    }
                    ExecutionResult::Val(t) => Ok(t.into_owned()),
                },
            },
        }
    }
}
pub trait HepAtomExt {
    #[allow(clippy::result_large_err)]
    fn parse_to_hep_net(
        &self,
        settings: &ParseSettings,
    ) -> Result<HepNet<AbstractIndex>, TensorNetworkError<ExplicitKey<AbstractIndex>, Symbol>>;
}

impl HepAtomExt for Atom {
    fn parse_to_hep_net(
        &self,
        settings: &ParseSettings,
    ) -> Result<HepNet<AbstractIndex>, TensorNetworkError<ExplicitKey<AbstractIndex>, Symbol>> {
        self.as_view().parse_to_hep_net(settings)
    }
}

impl HepAtomExt for AtomView<'_> {
    fn parse_to_hep_net(
        &self,
        settings: &ParseSettings,
    ) -> Result<HepNet<AbstractIndex>, TensorNetworkError<ExplicitKey<AbstractIndex>, Symbol>> {
        HepNet::try_from_view(*self, &*HEP_LIB, settings)
    }
}

#[allow(non_snake_case)]
#[allow(unused)]
pub fn A(
    i: impl Into<AbstractIndex>,
    j: impl Into<AbstractIndex>,
    k: impl Into<AbstractIndex>,
) -> Atom {
    let a_strct = IndexlessNamedStructure::<Symbol, ()>::from_iter(
        [
            Bispinor {}.new_rep(2).to_lib(),
            Bispinor {}.new_rep(2).cast(),
            Bispinor {}.new_rep(2).cast(),
        ],
        symbol!("A"),
        None,
    );
    a_strct
        .reindex([i.into(), j.into(), k.into()])
        .unwrap()
        .map_structure(|a| SymbolicTensor::from_named(&a).unwrap())
        .permute_inds()
        .expression
        .simplify_metrics()
}
#[allow(unused)]
#[allow(non_snake_case)]
pub fn B(
    i: impl Into<AbstractIndex>,
    j: impl Into<AbstractIndex>,
    k: impl Into<AbstractIndex>,
) -> Atom {
    let a_strct = IndexlessNamedStructure::<Symbol, ()>::from_iter(
        [
            Bispinor {}.new_rep(2).to_lib(),
            Bispinor {}.new_rep(2).cast(),
            Bispinor {}.new_rep(2).cast(),
        ],
        symbol!("B"),
        None,
    );
    a_strct
        .reindex([i.into(), j.into(), k.into()])
        .unwrap()
        .map_structure(|a| SymbolicTensor::from_named(&a).unwrap())
        .permute_inds()
        .expression
        .simplify_metrics()
}

pub fn gamma(
    i: impl Into<AbstractIndex>,
    j: impl Into<AbstractIndex>,
    mu: impl Into<AbstractIndex>,
) -> Atom {
    let gamma_strct = IndexlessNamedStructure::<Symbol, ()>::from_iter(
        [
            Bispinor {}.new_rep(4).to_lib(),
            Bispinor {}.new_rep(4).cast(),
            Minkowski {}.new_rep(4).cast(),
        ],
        AGS.gamma,
        None,
    );
    gamma_strct
        .reindex([i.into(), j.into(), mu.into()])
        .unwrap()
        .permute_with_metric()
}

pub fn gammaadj(
    i: impl Into<AbstractIndex>,
    j: impl Into<AbstractIndex>,
    mu: impl Into<AbstractIndex>,
) -> Atom {
    let gamma_strct = IndexlessNamedStructure::<Symbol, ()>::from_iter(
        [
            Bispinor {}.new_rep(4).to_lib(),
            Bispinor {}.new_rep(4).cast(),
            Minkowski {}.new_rep(4).cast(),
        ],
        AGS.gammaadj,
        None,
    );
    gamma_strct
        .reindex([i.into(), j.into(), mu.into()])
        .unwrap()
        .permute_with_metric()
}

pub fn gamma0(i: impl Into<AbstractIndex>, j: impl Into<AbstractIndex>) -> Atom {
    let gamma_strct = IndexlessNamedStructure::<Symbol, ()>::from_iter(
        [
            Bispinor {}.new_rep(4).to_lib(),
            Bispinor {}.new_rep(4).cast(),
        ],
        AGS.gamma0,
        None,
    );
    gamma_strct
        .reindex([i.into(), j.into()])
        .unwrap()
        .permute_with_metric()
}

pub fn gammaconj(
    i: impl Into<AbstractIndex>,
    j: impl Into<AbstractIndex>,
    mu: impl Into<AbstractIndex>,
) -> Atom {
    let gamma_strct = IndexlessNamedStructure::<Symbol, ()>::from_iter(
        [
            Bispinor {}.new_rep(4).to_lib(),
            Bispinor {}.new_rep(4).cast(),
            Minkowski {}.new_rep(4).cast(),
        ],
        AGS.gammaconj,
        None,
    );
    gamma_strct
        .reindex([i.into(), j.into(), mu.into()])
        .unwrap()
        .permute_with_metric()
}

pub fn p(m: impl Into<AbstractIndex>) -> Atom {
    let m_atom: AbstractIndex = m.into();
    let m_atom: Atom = m_atom.into();
    let mink = Minkowski {}.new_rep(4);
    function!(symbol!("spenso::p"), mink.to_symbolic([m_atom]))
}
pub fn q(m: impl Into<AbstractIndex>) -> Atom {
    let m_atom: AbstractIndex = m.into();
    let m_atom: Atom = m_atom.into();
    let mink = Minkowski {}.new_rep(4);
    function!(symbol!("spenso::q"), mink.to_symbolic([m_atom]))
}

pub fn u(i: usize, m: impl Into<AbstractIndex>) -> Atom {
    let m_atom: AbstractIndex = m.into();
    let m_atom: Atom = m_atom.into();
    let mink = Bispinor {}.new_rep(4);
    function!(symbol!("spenso::u"), i, mink.to_symbolic([m_atom]))
}

pub fn ub(i: usize, m: impl Into<AbstractIndex>) -> Atom {
    let m_atom: AbstractIndex = m.into();
    let m_atom: Atom = m_atom.into();
    let mink = Bispinor {}.new_rep(4);
    function!(symbol!("spenso::ub"), i, mink.to_symbolic([m_atom]))
}
