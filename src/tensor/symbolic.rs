use super::{
    Contract, HasName, IntoId, MixedTensor, Shadowable, Slot, StructureContract, TensorNetwork,
    TensorStructure, VecStructure,
};

use symbolica::representations::{Atom, AtomView, Symbol};

/// A fully symbolic tensor, with no concrete values.
///
/// This tensor is used to represent the structure of a tensor, and is used to perform symbolic contraction.
/// Currently contraction is just a multiplication of the atoms, but in the future this will ensure that internal indices are independent accross the contraction.
///
/// Additionally, this can also be used as a tensor structure, that tracks the history, much like [`HistoryStructure`].
#[derive(Debug)]
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

    pub fn to_network(self) -> Result<TensorNetwork<MixedTensor<VecStructure>>, &'static str> {
        let mut network: TensorNetwork<MixedTensor<VecStructure>> = TensorNetwork::new();

        if let AtomView::Mul(m) = self.expression.as_view() {
            for atom in m.iter() {
                if let AtomView::Fun(f) = atom {
                    let mut structure: Vec<Slot> = vec![];
                    let f_id = f.get_symbol();

                    for arg in f.iter() {
                        structure.push(arg.try_into()?);
                    }
                    let s: VecStructure = structure.into();
                    network.push(s.to_explicit_rep(f_id));
                }
            }
        }

        Ok(network)
    }
}

impl TryFrom<AtomView<'_>> for VecStructure {
    type Error = &'static str;
    fn try_from(value: AtomView) -> Result<Self, Self::Error> {
        let mut structure: Vec<Slot> = vec![];
        if let AtomView::Fun(f) = value {
            for arg in f.iter() {
                structure.push(arg.try_into()?);
            }
        }
        Ok(structure.into())
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
        println!("expression: {}", expression);
        new_structure.merge(&other.structure);
        Some(SymbolicTensor {
            expression,
            structure: new_structure,
        })
    }
}
