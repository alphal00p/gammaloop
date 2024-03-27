use super::{
    Contract, FallibleAdd, HasName, IntoId, MixedTensor, Shadowable, Slot, StructureContract,
    TensorNetwork, TensorStructure, VecStructure,
};

use symbolica::representations::{AddView, Atom, AtomView, MulView, Symbol};

/// A fully symbolic tensor, with no concrete values.
///
/// This tensor is used to represent the structure of a tensor, and is used to perform symbolic contraction.
/// Currently contraction is just a multiplication of the atoms, but in the future this will ensure that internal indices are independent accross the contraction.
///
/// Additionally, this can also be used as a tensor structure, that tracks the history, much like [`HistoryStructure`].
#[derive(Debug, Clone)]
pub struct SymbolicTensor {
    structure: VecStructure,
    expression: symbolica::representations::Atom,
}

impl TensorStructure for SymbolicTensor {
    type Structure = VecStructure;

    fn structure(&self) -> &Self::Structure {
        self.structure.structure()
    }

    fn mut_structure(&mut self) -> &mut Self::Structure {
        self.structure.mut_structure()
    }
    fn external_structure(&self) -> &[Slot] {
        self.structure.external_structure()
    }
}

impl StructureContract for SymbolicTensor {
    fn merge_at(&self, other: &Self, positions: (usize, usize)) -> Self {
        let structure = self.structure.merge_at(&other.structure, positions);
        // let mut out: Atom<Linear> = Atom::new();
        // other.expression.mul(state, ws, &self.expression, &mut out);

        SymbolicTensor {
            structure,
            expression: &other.expression * &self.expression,
        }
    }

    fn merge(&mut self, other: &Self) -> Option<usize> {
        self.expression = &other.expression * &self.expression;
        self.structure.merge(&other.structure)
    }

    fn trace_out(&mut self) {
        self.structure.trace_out();
    }

    fn trace(&mut self, i: usize, j: usize) {
        self.structure.trace(i, j);
    }
}

#[allow(dead_code)]
impl SymbolicTensor {
    pub fn from_named<N>(structure: &N) -> Option<Self>
    where
        N: TensorStructure + HasName,
        N::Name: IntoId + Clone,
    {
        Some(SymbolicTensor {
            expression: structure.to_symbolic()?,
            structure: structure.external_structure().to_vec().into(),
        })
    }

    #[must_use]
    pub fn get_atom(&self) -> &Atom {
        &self.expression
    }

    pub fn to_mixed(self) -> MixedTensor<VecStructure> {
        self.smart_shadow().unwrap()
    }

    fn add_view_to_tensor(add: AddView) -> Result<MixedTensor<VecStructure>, &'static str> {
        let mut terms = vec![];
        for t in add.iter() {
            match t {
                AtomView::Mul(m) => {
                    let mut term_net = Self::mul_to_network(m)?;
                    term_net.contract();
                    terms.push(term_net.result());
                }
                AtomView::Add(a) => {
                    terms.push(Self::add_view_to_tensor(a)?);
                }
                _ => {
                    return Err("Not a valid expression");
                }
            }
        }
        let sum = terms
            .into_iter()
            .reduce(|a, b| a.add_fallible(&b).unwrap())
            .unwrap();

        Ok(sum)
    }

    fn add_view_to_tracking_tensor(
        add: AddView,
    ) -> Result<MixedTensor<SymbolicTensor>, &'static str> {
        let mut terms = vec![];
        for t in add.iter() {
            match t {
                AtomView::Mul(m) => {
                    let mut term_net = Self::mul_to_tracking_network(m)?;
                    term_net.contract();
                    terms.push(term_net.result());
                }
                AtomView::Add(a) => {
                    terms.push(Self::add_view_to_tracking_tensor(a)?);
                }
                AtomView::Fun(f) => {
                    let mut a = Atom::new();
                    a.set_from_view(&f.as_view());
                    let structure: SymbolicTensor = a.try_into()?;

                    terms.push(structure.to_explicit_rep(f.get_symbol()));
                }
                _ => {
                    let mut a = Atom::new();
                    a.set_from_view(&t);
                    println!("Var: {}", a);
                    return Err("Not a valid expression");
                }
            }
        }
        let sum = terms
            .into_iter()
            .reduce(|a, b| a.add_fallible(&b).unwrap())
            .unwrap();

        Ok(sum)
    }

    pub fn mul_to_tracking_network(
        mul: MulView,
    ) -> Result<TensorNetwork<MixedTensor<SymbolicTensor>>, &'static str> {
        let mut network: TensorNetwork<MixedTensor<SymbolicTensor>> = TensorNetwork::new();
        for atom in mul.iter() {
            match atom {
                AtomView::Fun(f) => {
                    let mut a = Atom::new();
                    a.set_from_view(&f.as_view());
                    let structure: SymbolicTensor = a.try_into()?;

                    network.push(structure.to_explicit_rep(f.get_symbol()));
                }
                AtomView::Var(v) => {
                    let mut a = Atom::new();
                    a.set_from_view(&v.as_view());
                    network.scalar_mul(&a);
                }
                AtomView::Num(n) => {
                    let mut a = Atom::new();
                    a.set_from_view(&n.as_view());
                    network.scalar_mul(&a);
                }
                AtomView::Add(a) => {
                    let sum = Self::add_view_to_tracking_tensor(a)?;

                    network.push(sum);
                }
                _ => return Err("Not a valid expression"),
            }
        }
        Ok(network)
    }

    pub fn mul_to_network(
        mul: MulView,
    ) -> Result<TensorNetwork<MixedTensor<VecStructure>>, &'static str> {
        let mut network: TensorNetwork<MixedTensor<VecStructure>> = TensorNetwork::new();
        for atom in mul.iter() {
            match atom {
                AtomView::Fun(f) => {
                    let mut structure: Vec<Slot> = vec![];
                    let f_id = f.get_symbol();

                    for arg in f.iter() {
                        structure.push(arg.try_into()?);
                    }
                    let s: VecStructure = structure.into();
                    network.push(s.to_explicit_rep(f_id));
                }
                AtomView::Var(v) => {
                    let mut a = Atom::new();
                    a.set_from_view(&v.as_view());
                    network.scalar_mul(&a);
                }
                AtomView::Num(n) => {
                    let mut a = Atom::new();
                    a.set_from_view(&n.as_view());
                    network.scalar_mul(&a);
                }
                AtomView::Add(a) => {
                    let mut terms = vec![];
                    for t in a.iter() {
                        if let AtomView::Mul(m) = t {
                            let mut term_net = Self::mul_to_network(m)?;
                            term_net.contract();
                            terms.push(term_net.result());
                        }
                    }
                    let sum = terms
                        .into_iter()
                        .reduce(|a, b| a.add_fallible(&b).unwrap())
                        .unwrap();

                    network.push(sum);
                }
                _ => return Err("Not a valid expression"),
            }
        }
        Ok(network)
    }

    pub fn to_network(self) -> Result<TensorNetwork<MixedTensor<VecStructure>>, &'static str> {
        if let AtomView::Mul(mul) = self.expression.as_view() {
            Self::mul_to_network(mul)
        } else {
            Err("Not a valid expression")
        }
    }
}

impl TryFrom<Atom> for SymbolicTensor {
    type Error = &'static str;
    fn try_from(value: Atom) -> Result<Self, Self::Error> {
        Ok(SymbolicTensor {
            structure: value.as_view().try_into()?,
            expression: value,
        })
    }
}

impl HasName for SymbolicTensor {
    type Name = Symbol;
    fn name(&self) -> Option<std::borrow::Cow<Self::Name>> {
        if let AtomView::Fun(f) = self.expression.as_view() {
            Some(std::borrow::Cow::Owned(f.get_symbol()))
        } else {
            None
        }
    }

    fn set_name(&mut self, _name: &Self::Name) {
        unimplemented!("Cannot set name of a symbolic tensor")
    }
}

/// Symbolic contraction of two symbolic tensors is just a multiplication of the atoms.
///
impl Contract<SymbolicTensor> for SymbolicTensor {
    type LCM = SymbolicTensor;
    fn contract(&self, other: &SymbolicTensor) -> Option<Self::LCM> {
        let mut new_structure = self.structure.clone();

        let expression = &other.expression * &self.expression;
        new_structure.merge(&other.structure);
        Some(SymbolicTensor {
            expression,
            structure: new_structure,
        })
    }
}

impl std::fmt::Display for SymbolicTensor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.expression)
    }
}
