use ahash::AHashMap;
#[cfg(feature = "shadowing")]
use ahash::{AHashSet, HashMap};
#[cfg(feature = "shadowing")]
use anyhow::anyhow;
use bincode::{Decode, Encode};
use bitvec::vec::BitVec;
use graph::{NAdd, NMul, NetworkGraph, NetworkLeaf};
use linnet::half_edge::involution::EdgeData;
use linnet::half_edge::subgraph::{BaseSubgraph, ModifySubgraph, SubGraph};
use linnet::half_edge::tree::SimpleTraversalTree;
use linnet::half_edge::{
    builder::HedgeGraphBuilder,
    involution::{Flow, Hedge, Orientation},
};
use linnet::half_edge::{HedgeGraph, HedgeGraphError, NodeIndex};
use linnet::tree::child_pointer::ParentChildStore;
use rayon::iter::IntoParallelIterator;
use ref_ops::{RefAdd, RefMul};
use std::borrow::Cow;
use std::fmt::format;
use std::ops::{Add, Mul, Neg};
#[cfg(feature = "shadowing")]
use std::sync::Arc;
use store::{NetworkStore, TensorScalarStoreMapping};

#[cfg(feature = "shadowing")]
use symbolica::atom::PowView;

use crate::algebraic_traits::{One, Zero};
use crate::structure::representation::LibraryRep;
use crate::structure::slot::{IsAbstractSlot, Slot};
// use log::trace;
use serde::{Deserialize, Serialize};
use slotmap::{new_key_type, DenseSlotMap, Key, SecondaryMap};
#[cfg(feature = "shadowing")]
use symbolica::{
    atom::{representation::FunView, AddView, Atom, AtomView, MulView, Symbol},
    coefficient::ConvertToRing,
    domains::{
        factorized_rational_polynomial::{
            FactorizedRationalPolynomial, FromNumeratorAndFactorizedDenominator,
        },
        float::{Complex as SymComplex, NumericalFloatLike, Real, SingleFloat},
        rational::Rational,
        rational_polynomial::{FromNumeratorAndDenominator, RationalPolynomial},
        EuclideanDomain,
    },
    evaluate::{
        CompileOptions, CompiledCode, CompiledEvaluator, EvalTree, EvaluationFn, ExportedCode,
        ExpressionEvaluator, FunctionMap, InlineASM,
    },
    id::Pattern,
    poly::{factor::Factorize, gcd::PolynomialGCD, polynomial::MultivariatePolynomial, Variable},
    utils::BorrowedOrOwned,
};

#[cfg(feature = "shadowing")]
use symbolica::{
    atom::{AtomCore, KeyLookup},
    id::BorrowReplacement,
    poly::PositiveExponent,
    symbol,
};

#[cfg(feature = "shadowing")]
use crate::{
    complex::{Complex, RealOrComplexTensor},
    contraction::RefZero,
    data::{DataIterator, DenseTensor, SetTensorData, SparseTensor},
    iterators::IteratableTensor,
    parametric::atomcore::ReplaceBuilderGeneric,
    parametric::{
        atomcore::PatternReplacement, AtomViewOrConcrete, CompiledEvalTensor, EvalTensor,
        EvalTreeTensor, MixedTensor, ParamTensor, SerializableCompiledCode,
        SerializableCompiledEvaluator, SerializableExportedCode,
    },
    shadowing::{ShadowMapping, Shadowable},
    structure::{StructureContract, ToSymbolic},
    symbolica_utils::{IntoArgs, IntoSymbol, SerializableAtom},
    tensor_library::{LibraryTensor, TensorLibrary},
    upgrading_arithmetic::{FallibleAdd, TrySmallestUpgrade},
};

use crate::{
    arithmetic::ScalarMul,
    contraction::{Contract, ContractionError, Trace},
    data::{DataTensor, GetTensorData, HasTensorData},
    structure::{
        representation::{LibrarySlot, RepName},
        slot::DualSlotTo,
        CastStructure, HasName, HasStructure, ScalarTensor, TensorStructure, TracksCount,
    },
    upgrading_arithmetic::FallibleMul,
};

#[cfg(feature = "shadowing")]
use crate::{
    data::StorageTensor,
    parametric::atomcore::{TensorAtomMaps, TensorAtomOps},
};

use super::store::TensorScalarStore;
use super::Network;
