use serde::{Deserialize, Serialize};

use super::{
    graph::ScalarRef,
    profile::{self, Counter, Timer},
};

pub trait TensorScalarStore: Default + TensorScalarStoreMapping {
    fn add_tensor(&mut self, tensor: Self::Tensor) -> usize;
    fn add_scalar(&mut self, scalar: Self::Scalar) -> usize;

    fn get_scalar(&self, index: usize) -> &Self::Scalar;
    fn get_scalar_ref(&self, scalar: ScalarRef) -> &Self::Scalar {
        self.get_scalar(scalar.index())
    }
    fn get_tensor(&self, index: usize) -> &Self::Tensor;

    fn n_tensors(&self) -> usize;
    fn n_scalars(&self) -> usize;
    fn extend(&mut self, other: Self);
}

#[doc(hidden)]
pub trait NetworkStoreAccess {
    type Tensor;
    type Scalar;

    fn tensor(&self, index: usize) -> &Self::Tensor;
    fn scalar(&self, index: usize) -> &Self::Scalar;
    fn scalar_ref(&self, scalar: ScalarRef) -> &Self::Scalar {
        self.scalar(scalar.index())
    }
    fn push_tensor(&mut self, tensor: Self::Tensor) -> usize;
    fn push_scalar(&mut self, scalar: Self::Scalar) -> usize;
}

pub trait TensorScalarStoreMapping: Sized {
    type Tensor;

    type Scalar;
    type Store<T, S>: TensorScalarStoreMapping<Tensor = T, Scalar = S>;
    fn iter_scalars(&self) -> impl Iterator<Item = &Self::Scalar>;

    fn iter_tensors(&self) -> impl Iterator<Item = &Self::Tensor>;
    fn iter_scalars_mut(&mut self) -> impl Iterator<Item = &mut Self::Scalar>;
    fn iter_tensors_mut(&mut self) -> impl Iterator<Item = &mut Self::Tensor>;
    fn map<U, V>(
        self,
        scalar_map: impl FnMut(Self::Scalar) -> U,
        tensor_map: impl FnMut(Self::Tensor) -> V,
    ) -> Self::Store<V, U>;

    // fn map_self(
    //     self,
    //     scalar_map: impl FnMut(Self::Scalar) -> Self::Scalar,
    //     tensor_map: impl FnMut(Self::Tensor) -> Self::Tensor,
    // ) -> Self;

    fn map_result<U, V, Er>(
        self,
        scalar_map: impl FnMut(Self::Scalar) -> Result<U, Er>,
        tensor_map: impl FnMut(Self::Tensor) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er>;

    // fn map_result_self<Er>(
    //     self,
    //     scalar_map: impl FnMut(Self::Scalar) -> Result<Self::Scalar, Er>,
    //     tensor_map: impl FnMut(Self::Tensor) -> Result<Self::Tensor, Er>,
    // ) -> Result<Self, Er>;

    fn map_ref<'a, U, V>(
        &'a self,
        scalar_map: impl FnMut(&'a Self::Scalar) -> U,
        tensor_map: impl FnMut(&'a Self::Tensor) -> V,
    ) -> Self::Store<V, U>;

    // fn map_ref_self(
    //     &self,
    //     scalar_map: impl FnMut(&Self::Scalar) -> Self::Scalar,
    //     tensor_map: impl FnMut(&Self::Tensor) -> Self::Tensor,
    // ) -> Self;

    fn map_ref_enumerate<U, V>(
        &self,
        scalar_map: impl FnMut((usize, &Self::Scalar)) -> U,
        tensor_map: impl FnMut((usize, &Self::Tensor)) -> V,
    ) -> Self::Store<V, U>;

    // fn map_ref_enumerate_self(
    //     &self,
    //     scalar_map: impl FnMut((usize, &Self::Scalar)) -> Self::Scalar,
    //     tensor_map: impl FnMut((usize, &Self::Tensor)) -> Self::Tensor,
    // ) -> Self;

    fn map_ref_result<U, V, Er>(
        &self,
        scalar_map: impl FnMut(&Self::Scalar) -> Result<U, Er>,
        tensor_map: impl FnMut(&Self::Tensor) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er>;

    // fn map_ref_result_self<Er>(
    //     &self,
    //     scalar_map: impl FnMut(&Self::Scalar) -> Result<Self::Scalar, Er>,
    //     tensor_map: impl FnMut(&Self::Tensor) -> Result<Self::Tensor, Er>,
    // ) -> Result<Self::Store<V, U>, Er>;

    fn map_ref_result_enumerate<U, V, Er>(
        &self,
        scalar_map: impl FnMut((usize, &Self::Scalar)) -> Result<U, Er>,
        tensor_map: impl FnMut((usize, &Self::Tensor)) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er>;

    fn map_ref_mut<U, V>(
        &mut self,
        scalar_map: impl FnMut(&mut Self::Scalar) -> U,
        tensor_map: impl FnMut(&mut Self::Tensor) -> V,
    ) -> Self::Store<V, U>;

    fn map_ref_mut_enumerate<U, V>(
        &mut self,
        scalar_map: impl FnMut((usize, &mut Self::Scalar)) -> U,
        tensor_map: impl FnMut((usize, &mut Self::Tensor)) -> V,
    ) -> Self::Store<V, U>;
    fn map_ref_mut_result<U, V, Er>(
        &mut self,
        scalar_map: impl FnMut(&mut Self::Scalar) -> Result<U, Er>,
        tensor_map: impl FnMut(&mut Self::Tensor) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er>;

    fn map_ref_mut_result_enumerate<U, V, Er>(
        &mut self,
        scalar_map: impl FnMut((usize, &mut Self::Scalar)) -> Result<U, Er>,
        tensor_map: impl FnMut((usize, &mut Self::Tensor)) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er>;
}

#[derive(
    Debug,
    Clone,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
    PartialEq,
    Eq,
)]
pub struct NetworkStore<T, S> {
    pub tensors: Vec<T>,
    // pub params: AHashSet<Atom>,
    pub scalar: Vec<S>,
    pub scalar_aliases: Vec<Option<S>>,
}

#[doc(hidden)]
pub struct NetworkStoreOverlay<'a, T, S> {
    base_tensors: &'a [T],
    base_scalars: &'a [S],
    base_scalar_aliases: &'a [Option<S>],
    pub tensors: Vec<T>,
    pub scalar: Vec<S>,
    pub scalar_aliases: Vec<Option<S>>,
}

impl<'a, T, S> NetworkStoreOverlay<'a, T, S> {
    pub fn new(base: &'a NetworkStore<T, S>) -> Self {
        Self {
            base_tensors: &base.tensors,
            base_scalars: &base.scalar,
            base_scalar_aliases: &base.scalar_aliases,
            tensors: Vec::new(),
            scalar: Vec::new(),
            scalar_aliases: Vec::new(),
        }
    }

    pub fn into_additions(self) -> (Vec<T>, Vec<S>, Vec<Option<S>>) {
        (self.tensors, self.scalar, self.scalar_aliases)
    }
}

impl<T, S> Default for NetworkStore<T, S> {
    fn default() -> Self {
        NetworkStore {
            tensors: vec![],
            scalar: vec![],
            scalar_aliases: vec![],
        }
    }
}

impl<T, S> NetworkStoreAccess for NetworkStore<T, S> {
    type Tensor = T;
    type Scalar = S;

    fn tensor(&self, index: usize) -> &Self::Tensor {
        &self.tensors[index]
    }

    fn scalar(&self, index: usize) -> &Self::Scalar {
        &self.scalar[index]
    }

    fn scalar_ref(&self, scalar: ScalarRef) -> &Self::Scalar {
        match scalar {
            ScalarRef::Store(index) => &self.scalar[index],
            ScalarRef::Alias(index) => self
                .scalar_aliases
                .get(index)
                .and_then(Option::as_ref)
                .unwrap_or(&self.scalar[index]),
        }
    }

    fn push_tensor(&mut self, tensor: Self::Tensor) -> usize {
        let index = self.tensors.len();
        self.tensors.push(tensor);
        index
    }

    fn push_scalar(&mut self, scalar: Self::Scalar) -> usize {
        let index = self.scalar.len();
        self.scalar.push(scalar);
        self.scalar_aliases.push(None);
        index
    }
}

impl<T, S> NetworkStoreAccess for NetworkStoreOverlay<'_, T, S> {
    type Tensor = T;
    type Scalar = S;

    fn tensor(&self, index: usize) -> &Self::Tensor {
        if index < self.base_tensors.len() {
            &self.base_tensors[index]
        } else {
            &self.tensors[index - self.base_tensors.len()]
        }
    }

    fn scalar(&self, index: usize) -> &Self::Scalar {
        if index < self.base_scalars.len() {
            &self.base_scalars[index]
        } else {
            &self.scalar[index - self.base_scalars.len()]
        }
    }

    fn scalar_ref(&self, scalar: ScalarRef) -> &Self::Scalar {
        let index = scalar.index();
        match scalar {
            ScalarRef::Store(_) => self.scalar(index),
            ScalarRef::Alias(_) => {
                if index < self.base_scalars.len() {
                    self.base_scalar_aliases
                        .get(index)
                        .and_then(Option::as_ref)
                        .unwrap_or(&self.base_scalars[index])
                } else {
                    let local = index - self.base_scalars.len();
                    self.scalar_aliases
                        .get(local)
                        .and_then(Option::as_ref)
                        .unwrap_or(&self.scalar[local])
                }
            }
        }
    }

    fn push_tensor(&mut self, tensor: Self::Tensor) -> usize {
        let index = self.base_tensors.len() + self.tensors.len();
        self.tensors.push(tensor);
        index
    }

    fn push_scalar(&mut self, scalar: Self::Scalar) -> usize {
        let index = self.base_scalars.len() + self.scalar.len();
        self.scalar.push(scalar);
        self.scalar_aliases.push(None);
        index
    }
}

impl<T, S> TensorScalarStore for NetworkStore<T, S> {
    fn n_tensors(&self) -> usize {
        self.tensors.len()
    }

    fn n_scalars(&self) -> usize {
        self.scalar.len()
    }

    fn extend(&mut self, other: Self) {
        let _span = profile::span(Timer::StoreExtend);
        profile::bump(Counter::StoreExtend, 1);
        self.tensors.extend(other.tensors);
        self.scalar.extend(other.scalar);
        self.scalar_aliases.extend(other.scalar_aliases);
    }

    fn add_scalar(&mut self, scalar: Self::Scalar) -> usize {
        let id = self.scalar.len();
        self.scalar.push(scalar);
        self.scalar_aliases.push(None);
        id
    }

    fn add_tensor(&mut self, tensor: Self::Tensor) -> usize {
        let id = self.tensors.len();
        self.tensors.push(tensor);
        id
    }

    fn get_scalar(&self, index: usize) -> &Self::Scalar {
        &self.scalar[index]
    }

    fn get_scalar_ref(&self, scalar: ScalarRef) -> &Self::Scalar {
        match scalar {
            ScalarRef::Store(index) => &self.scalar[index],
            ScalarRef::Alias(index) => self
                .scalar_aliases
                .get(index)
                .and_then(Option::as_ref)
                .unwrap_or(&self.scalar[index]),
        }
    }

    fn get_tensor(&self, index: usize) -> &Self::Tensor {
        &self.tensors[index]
    }
}

impl<T, S> TensorScalarStoreMapping for NetworkStore<T, S> {
    type Store<U, V> = NetworkStore<U, V>;
    type Scalar = S;
    type Tensor = T;

    fn iter_scalars(&self) -> impl Iterator<Item = &Self::Scalar> {
        self.scalar.iter()
    }

    fn iter_tensors(&self) -> impl Iterator<Item = &Self::Tensor> {
        self.tensors.iter()
    }

    fn iter_scalars_mut(&mut self) -> impl Iterator<Item = &mut Self::Scalar> {
        self.scalar.iter_mut()
    }

    fn iter_tensors_mut(&mut self) -> impl Iterator<Item = &mut Self::Tensor> {
        self.tensors.iter_mut()
    }

    fn map<U, V>(
        self,
        mut scalar_map: impl FnMut(S) -> U,
        tensor_map: impl FnMut(T) -> V,
    ) -> NetworkStore<V, U> {
        NetworkStore {
            tensors: self.tensors.into_iter().map(tensor_map).collect(),
            scalar: self.scalar.into_iter().map(&mut scalar_map).collect(),
            scalar_aliases: self
                .scalar_aliases
                .into_iter()
                .map(|alias| alias.map(&mut scalar_map))
                .collect(),
        }
    }

    fn map_result<U, V, Er>(
        self,
        mut scalar_map: impl FnMut(S) -> Result<U, Er>,
        tensor_map: impl FnMut(T) -> Result<V, Er>,
    ) -> Result<NetworkStore<V, U>, Er> {
        Ok(NetworkStore {
            tensors: self
                .tensors
                .into_iter()
                .map(tensor_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar: self
                .scalar
                .into_iter()
                .map(&mut scalar_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar_aliases: self
                .scalar_aliases
                .into_iter()
                .map(|alias| alias.map(&mut scalar_map).transpose())
                .collect::<Result<Vec<_>, Er>>()?,
        })
    }

    fn map_ref<'a, U, V>(
        &'a self,
        mut scalar_map: impl FnMut(&'a S) -> U,
        tensor_map: impl FnMut(&'a T) -> V,
    ) -> NetworkStore<V, U> {
        NetworkStore {
            tensors: self.tensors.iter().map(tensor_map).collect(),
            scalar: self.scalar.iter().map(&mut scalar_map).collect(),
            scalar_aliases: self
                .scalar_aliases
                .iter()
                .map(|alias| alias.as_ref().map(&mut scalar_map))
                .collect(),
        }
    }

    fn map_ref_enumerate<U, V>(
        &self,
        mut scalar_map: impl FnMut((usize, &S)) -> U,
        tensor_map: impl FnMut((usize, &T)) -> V,
    ) -> NetworkStore<V, U> {
        NetworkStore {
            tensors: self.tensors.iter().enumerate().map(tensor_map).collect(),
            scalar: self
                .scalar
                .iter()
                .enumerate()
                .map(&mut scalar_map)
                .collect(),
            scalar_aliases: self
                .scalar_aliases
                .iter()
                .enumerate()
                .map(|(index, alias)| alias.as_ref().map(|scalar| scalar_map((index, scalar))))
                .collect(),
        }
    }

    fn map_ref_result<U, V, Er>(
        &self,
        mut scalar_map: impl FnMut(&S) -> Result<U, Er>,
        tensor_map: impl FnMut(&T) -> Result<V, Er>,
    ) -> Result<NetworkStore<V, U>, Er> {
        Ok(NetworkStore {
            tensors: self
                .tensors
                .iter()
                .map(tensor_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar: self
                .scalar
                .iter()
                .map(&mut scalar_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar_aliases: self
                .scalar_aliases
                .iter()
                .map(|alias| alias.as_ref().map(&mut scalar_map).transpose())
                .collect::<Result<Vec<_>, Er>>()?,
        })
    }

    fn map_ref_result_enumerate<U, V, Er>(
        &self,
        mut scalar_map: impl FnMut((usize, &S)) -> Result<U, Er>,
        tensor_map: impl FnMut((usize, &T)) -> Result<V, Er>,
    ) -> Result<NetworkStore<V, U>, Er> {
        Ok(NetworkStore {
            tensors: self
                .tensors
                .iter()
                .enumerate()
                .map(tensor_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar: self
                .scalar
                .iter()
                .enumerate()
                .map(&mut scalar_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar_aliases: self
                .scalar_aliases
                .iter()
                .enumerate()
                .map(|(index, alias)| {
                    alias
                        .as_ref()
                        .map(|scalar| scalar_map((index, scalar)))
                        .transpose()
                })
                .collect::<Result<Vec<_>, Er>>()?,
        })
    }

    fn map_ref_mut<U, V>(
        &mut self,
        mut scalar_map: impl FnMut(&mut S) -> U,
        tensor_map: impl FnMut(&mut T) -> V,
    ) -> NetworkStore<V, U> {
        NetworkStore {
            tensors: self.tensors.iter_mut().map(tensor_map).collect(),
            scalar: self.scalar.iter_mut().map(&mut scalar_map).collect(),
            scalar_aliases: self
                .scalar_aliases
                .iter_mut()
                .map(|alias| alias.as_mut().map(&mut scalar_map))
                .collect(),
        }
    }

    fn map_ref_mut_enumerate<U, V>(
        &mut self,
        mut scalar_map: impl FnMut((usize, &mut S)) -> U,
        tensor_map: impl FnMut((usize, &mut T)) -> V,
    ) -> NetworkStore<V, U> {
        NetworkStore {
            tensors: self
                .tensors
                .iter_mut()
                .enumerate()
                .map(tensor_map)
                .collect(),
            scalar: self
                .scalar
                .iter_mut()
                .enumerate()
                .map(&mut scalar_map)
                .collect(),
            scalar_aliases: self
                .scalar_aliases
                .iter_mut()
                .enumerate()
                .map(|(index, alias)| alias.as_mut().map(|scalar| scalar_map((index, scalar))))
                .collect(),
        }
    }

    fn map_ref_mut_result<U, V, Er>(
        &mut self,
        mut scalar_map: impl FnMut(&mut S) -> Result<U, Er>,
        tensor_map: impl FnMut(&mut T) -> Result<V, Er>,
    ) -> Result<NetworkStore<V, U>, Er> {
        Ok(NetworkStore {
            tensors: self
                .tensors
                .iter_mut()
                .map(tensor_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar: self
                .scalar
                .iter_mut()
                .map(&mut scalar_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar_aliases: self
                .scalar_aliases
                .iter_mut()
                .map(|alias| alias.as_mut().map(&mut scalar_map).transpose())
                .collect::<Result<Vec<_>, Er>>()?,
        })
    }

    fn map_ref_mut_result_enumerate<U, V, Er>(
        &mut self,
        mut scalar_map: impl FnMut((usize, &mut S)) -> Result<U, Er>,
        tensor_map: impl FnMut((usize, &mut T)) -> Result<V, Er>,
    ) -> Result<NetworkStore<V, U>, Er> {
        Ok(NetworkStore {
            tensors: self
                .tensors
                .iter_mut()
                .enumerate()
                .map(tensor_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar: self
                .scalar
                .iter_mut()
                .enumerate()
                .map(&mut scalar_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar_aliases: self
                .scalar_aliases
                .iter_mut()
                .enumerate()
                .map(|(index, alias)| {
                    alias
                        .as_mut()
                        .map(|scalar| scalar_map((index, scalar)))
                        .transpose()
                })
                .collect::<Result<Vec<_>, Er>>()?,
        })
    }
}
