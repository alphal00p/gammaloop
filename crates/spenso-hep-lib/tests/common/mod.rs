use idenso::{
    dirac::AGS,
    representations::{Bispinor, initialize},
    shorthands::metric::PermuteWithMetric,
};
use spenso::{
    network::{
        ExecutionResult, Sequential, SmallestDegree, TensorNetworkError,
        library::symbolic::ExplicitKey,
        parsing::{ParseSettings, StrictTensorFilter},
    },
    structure::{
        IndexlessNamedStructure, ScalarTensor, TensorStructure,
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName},
    },
    tensors::parametric::ParamTensor,
};
use spenso_hep_lib::{FUN_LIB, HEP_LIB, HepNet, HepTensor};
use symbolica::{
    atom::{Atom, AtomView, Symbol},
    function, symbol,
};

pub fn test_initialize() {
    let gamma = AGS.gamma_strct(4);
    let _a = HEP_LIB.get(gamma.canonical()).unwrap();

    initialize();
}
#[allow(unused)]
#[allow(clippy::result_large_err)]
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
        match self.result_tensor(&*HEP_LIB)? {
            ExecutionResult::One => Ok(HepTensor::Param(ParamTensor::new_scalar(Atom::one()))),
            ExecutionResult::Zero => Ok(HepTensor::Param(ParamTensor::new_scalar(Atom::zero()))),
            ExecutionResult::Val(tensor) => Ok(tensor.into_owned()),
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
        let settings = settings
            .clone()
            .with_strict_tensor_filter(StrictTensorFilter::ContainsReps);

        HepNet::try_from_view(*self, &*HEP_LIB, &settings)
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
    let logical_indices: [AbstractIndex; 3] = [i.into(), j.into(), k.into()];
    let storage_indices = a_strct.layout().logical_to_canonical(&logical_indices);
    a_strct
        .into_canonical()
        .reindex_storage(&storage_indices)
        .unwrap()
        .permute_with_metric()
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
    let logical_indices: [AbstractIndex; 3] = [i.into(), j.into(), k.into()];
    let storage_indices = a_strct.layout().logical_to_canonical(&logical_indices);
    a_strct
        .into_canonical()
        .reindex_storage(&storage_indices)
        .unwrap()
        .permute_with_metric()
}

#[allow(unused)]
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
    let logical_indices: [AbstractIndex; 3] = [i.into(), j.into(), mu.into()];
    let storage_indices = gamma_strct.layout().logical_to_canonical(&logical_indices);
    gamma_strct
        .into_canonical()
        .reindex_storage(&storage_indices)
        .unwrap()
        .permute_with_metric()
}
#[allow(unused)]
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
    let logical_indices: [AbstractIndex; 3] = [i.into(), j.into(), mu.into()];
    let storage_indices = gamma_strct.layout().logical_to_canonical(&logical_indices);
    gamma_strct
        .into_canonical()
        .reindex_storage(&storage_indices)
        .unwrap()
        .permute_with_metric()
}

#[allow(unused)]
pub fn gamma0(i: impl Into<AbstractIndex>, j: impl Into<AbstractIndex>) -> Atom {
    let gamma_strct = IndexlessNamedStructure::<Symbol, ()>::from_iter(
        [
            Bispinor {}.new_rep(4).to_lib(),
            Bispinor {}.new_rep(4).cast(),
        ],
        AGS.gamma0,
        None,
    );
    let logical_indices: [AbstractIndex; 2] = [i.into(), j.into()];
    let storage_indices = gamma_strct.layout().logical_to_canonical(&logical_indices);
    gamma_strct
        .into_canonical()
        .reindex_storage(&storage_indices)
        .unwrap()
        .permute_with_metric()
}
#[allow(unused)]
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
    let logical_indices: [AbstractIndex; 3] = [i.into(), j.into(), mu.into()];
    let storage_indices = gamma_strct.layout().logical_to_canonical(&logical_indices);
    gamma_strct
        .into_canonical()
        .reindex_storage(&storage_indices)
        .unwrap()
        .permute_with_metric()
}
#[allow(unused)]
pub fn p(m: impl Into<AbstractIndex>) -> Atom {
    let m_atom: AbstractIndex = m.into();
    let m_atom: Atom = m_atom.into();
    let mink = Minkowski {}.new_rep(4);
    function!(symbol!("spenso::p"), mink.to_symbolic([m_atom]))
}
#[allow(unused)]
pub fn q(m: impl Into<AbstractIndex>) -> Atom {
    let m_atom: AbstractIndex = m.into();
    let m_atom: Atom = m_atom.into();
    let mink = Minkowski {}.new_rep(4);
    function!(symbol!("spenso::q"), mink.to_symbolic([m_atom]))
}
#[allow(unused)]
pub fn u(i: usize, m: impl Into<AbstractIndex>) -> Atom {
    let m_atom: AbstractIndex = m.into();
    let m_atom: Atom = m_atom.into();
    let mink = Bispinor {}.new_rep(4);
    function!(symbol!("spenso::u"), i, mink.to_symbolic([m_atom]))
}
#[allow(unused)]
pub fn ub(i: usize, m: impl Into<AbstractIndex>) -> Atom {
    let m_atom: AbstractIndex = m.into();
    let m_atom: Atom = m_atom.into();
    let mink = Bispinor {}.new_rep(4);
    function!(symbol!("spenso::ub"), i, mink.to_symbolic([m_atom]))
}
