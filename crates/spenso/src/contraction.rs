use bitvec::vec::BitVec;
use log::trace;

use crate::{
    algebra::{
        algebraic_traits::{IsZero, RefZero},
        upgrading_arithmetic::{
            FallibleAddAssign, FallibleMul, FallibleSubAssign, TrySmallestUpgrade,
        },
    },
    structure::{
        HasStructure, MergeInfo, StructureContract, StructureError, TensorStructure,
        dimension::DimensionError,
    },
    tensors::data::DataTensor,
};

use thiserror::Error;

#[derive(Error, Debug)]
pub enum ContractionError {
    #[error("Sparse tensor is empty")]
    EmptySparse,
    #[error("Structure Error:{0}")]
    StructureError(#[from] StructureError),
    #[error("Dimension Error:{0}")]
    DimensionError(#[from] DimensionError),
    #[error(transparent)]
    Other(#[from] eyre::Error),
}

pub trait Contract<T = Self> {
    type LCM;
    fn contract(&self, other: &T) -> Result<Self::LCM, ContractionError>;
}

pub trait Trace {
    #[must_use]
    fn internal_contract(&self) -> Self;
}

pub trait ExteriorProduct<T> {
    type LCM: HasStructure;
    fn exterior_product(
        &self,
        other: &T,
        result_structure: <Self::LCM as HasStructure>::Structure,
    ) -> Result<Self::LCM, ContractionError>;
}

pub trait ExteriorProductInterleaved<T> {
    type LCM: HasStructure;
    fn exterior_product_interleaved(
        &self,
        other: &T,
        result_structure: <Self::LCM as HasStructure>::Structure,
        partition: BitVec,
    ) -> Result<Self::LCM, ContractionError>;
}

pub trait SingleContract<T> {
    type LCM: HasStructure;
    fn single_contract(
        &self,
        other: &T,
        resulting_structure: <Self::LCM as HasStructure>::Structure,
        i: usize,
        j: usize,
    ) -> Result<Self::LCM, ContractionError>;
}

pub trait SingleContractInterleaved<T> {
    type LCM: HasStructure;
    fn single_contract_interleaved(
        &self,
        other: &T,
        resulting_structure: <Self::LCM as HasStructure>::Structure,
        resulting_partition: BitVec,
        i: usize,
        j: usize,
    ) -> Result<Self::LCM, ContractionError>;
}

pub trait MultiContract<T> {
    type LCM: HasStructure;
    fn multi_contract(
        &self,
        other: &T,
        resulting_structure: <Self::LCM as HasStructure>::Structure,
    ) -> Result<Self::LCM, ContractionError>;
}

pub trait MultiContractInterleaved<T> {
    type LCM: HasStructure;
    fn multi_contract_interleaved(
        &self,
        other: &T,
        pos_self: BitVec,
        pos_other: BitVec,
        resulting_structure: <Self::LCM as HasStructure>::Structure,
        resulting_partition: BitVec,
    ) -> Result<Self::LCM, ContractionError>;
}

pub trait ContractableWith<T>
where
    Self: FallibleMul<T, Output = Self::Out> + Sized + TrySmallestUpgrade<T, LCM = Self::Out>,
{
    type Out: FallibleAddAssign<Self::Out> + FallibleSubAssign<Self::Out> + Clone + RefZero;
}

impl<T, U, Out> ContractableWith<T> for U
where
    U: FallibleMul<T, Output = Out> + TrySmallestUpgrade<T, LCM = Out>,
    // T: FallibleMul<U, Output = Out>,
    Out: FallibleAddAssign<Out> + FallibleSubAssign<Out> + Clone + RefZero,
{
    type Out = Out;
}

pub mod exteriorproduct;
pub mod multicontract;
pub mod singlecontract;
pub mod trace;

impl<T, U, O> Contract<T> for U
where
    U: SingleContract<T, LCM = O>
        + SingleContractInterleaved<T, LCM = O>
        + MultiContract<T, LCM = O>
        + MultiContractInterleaved<T, LCM = O>
        + ExteriorProduct<T, LCM = O>
        + ExteriorProductInterleaved<T, LCM = O>
        + HasStructure,
    // U::Structure: TensorStructure+ScalarStructure,
    T: HasStructure<Structure = U::Structure>,
    T: SingleContract<U, LCM = O>
        + SingleContractInterleaved<U, LCM = O>
        + MultiContract<U, LCM = O>
        + MultiContractInterleaved<U, LCM = O>
        + ExteriorProduct<U, LCM = O>
        + ExteriorProductInterleaved<U, LCM = O>
        + HasStructure<Structure = U::Structure>,
    U::Structure: StructureContract,
    O: HasStructure<Structure = U::Structure>,
{
    type LCM = O;
    fn contract(&self, other: &T) -> Result<Self::LCM, ContractionError> {
        let (resulting_structure, pos_self, pos_other, mergeinfo) =
            self.structure().merge(other.structure())?;

        // Count bits set to true in the BitVec
        let common_count_self = pos_self.count_ones();

        match common_count_self {
            0 => {
                trace!("exterior");
                match mergeinfo {
                    MergeInfo::Interleaved(partition) => {
                        self.exterior_product_interleaved(other, resulting_structure, partition)
                    }
                    MergeInfo::FirstBeforeSecond => {
                        self.exterior_product(other, resulting_structure)
                    }
                    MergeInfo::SecondBeforeFirst => {
                        other.exterior_product(self, resulting_structure)
                    }
                }
            }
            1 => {
                // Find indices where bits are set to true
                let i = pos_self.first_one().expect("Expected a bit to be set");
                let j = pos_other.first_one().expect("Expected a bit to be set");
                trace!("single");
                match mergeinfo {
                    MergeInfo::Interleaved(partition) => self.single_contract_interleaved(
                        other,
                        resulting_structure,
                        partition,
                        i,
                        j,
                    ),
                    MergeInfo::FirstBeforeSecond => {
                        self.single_contract(other, resulting_structure, i, j)
                    }
                    MergeInfo::SecondBeforeFirst => {
                        other.single_contract(self, resulting_structure, j, i)
                    }
                }
            }
            _ => {
                trace!("multi");
                match mergeinfo {
                    MergeInfo::Interleaved(partition) => self.multi_contract_interleaved(
                        other,
                        pos_self,
                        pos_other,
                        resulting_structure,
                        partition,
                    ),
                    MergeInfo::FirstBeforeSecond => self.multi_contract(other, resulting_structure),
                    MergeInfo::SecondBeforeFirst => other.multi_contract(self, resulting_structure),
                }
            }
        }
    }
}

impl<T, U, I, O> Contract<DataTensor<T, I>> for DataTensor<U, I>
where
    U: ContractableWith<T, Out = O>,
    T: ContractableWith<U, Out = O>,
    O: FallibleAddAssign<O> + FallibleSubAssign<O> + Clone + RefZero + IsZero,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DataTensor<U::Out, I>;
    fn contract(&self, other: &DataTensor<T, I>) -> Result<Self::LCM, ContractionError> {
        match (self, other) {
            (DataTensor::Dense(s), DataTensor::Dense(o)) => Ok(DataTensor::Dense(s.contract(o)?)),
            (DataTensor::Dense(s), DataTensor::Sparse(o)) => Ok(DataTensor::Dense(s.contract(o)?)),
            (DataTensor::Sparse(s), DataTensor::Dense(o)) => Ok(DataTensor::Dense(s.contract(o)?)),
            (DataTensor::Sparse(s), DataTensor::Sparse(o)) => {
                Ok(DataTensor::Sparse(s.contract(o)?))
            }
        }
    }
}
