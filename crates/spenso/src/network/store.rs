use serde::{Deserialize, Serialize};

pub trait TensorScalarStore: Default + TensorScalarStoreMapping {
    fn add_tensor(&mut self, tensor: Self::Tensor) -> usize;
    fn add_scalar(&mut self, scalar: Self::Scalar) -> usize;

    fn get_scalar(&self, index: usize) -> &Self::Scalar;
    fn get_tensor(&self, index: usize) -> &Self::Tensor;

    fn n_tensors(&self) -> usize;
    fn n_scalars(&self) -> usize;
    fn extend(&mut self, other: Self);
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
}

impl<T, S> Default for NetworkStore<T, S> {
    fn default() -> Self {
        NetworkStore {
            tensors: vec![],
            scalar: vec![],
        }
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
        self.tensors.extend(other.tensors);
        self.scalar.extend(other.scalar);
    }

    fn add_scalar(&mut self, scalar: Self::Scalar) -> usize {
        let id = self.scalar.len();
        self.scalar.push(scalar);
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
        scalar_map: impl FnMut(S) -> U,
        tensor_map: impl FnMut(T) -> V,
    ) -> NetworkStore<V, U> {
        NetworkStore {
            tensors: self.tensors.into_iter().map(tensor_map).collect(),
            scalar: self.scalar.into_iter().map(scalar_map).collect(),
        }
    }

    fn map_result<U, V, Er>(
        self,
        scalar_map: impl FnMut(S) -> Result<U, Er>,
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
                .map(scalar_map)
                .collect::<Result<Vec<_>, Er>>()?,
        })
    }

    fn map_ref<'a, U, V>(
        &'a self,
        scalar_map: impl FnMut(&'a S) -> U,
        tensor_map: impl FnMut(&'a T) -> V,
    ) -> NetworkStore<V, U> {
        NetworkStore {
            tensors: self.tensors.iter().map(tensor_map).collect(),
            scalar: self.scalar.iter().map(scalar_map).collect(),
        }
    }

    fn map_ref_enumerate<U, V>(
        &self,
        scalar_map: impl FnMut((usize, &S)) -> U,
        tensor_map: impl FnMut((usize, &T)) -> V,
    ) -> NetworkStore<V, U> {
        NetworkStore {
            tensors: self.tensors.iter().enumerate().map(tensor_map).collect(),
            scalar: self.scalar.iter().enumerate().map(scalar_map).collect(),
        }
    }

    fn map_ref_result<U, V, Er>(
        &self,
        scalar_map: impl FnMut(&S) -> Result<U, Er>,
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
                .map(scalar_map)
                .collect::<Result<Vec<_>, Er>>()?,
        })
    }

    fn map_ref_result_enumerate<U, V, Er>(
        &self,
        scalar_map: impl FnMut((usize, &S)) -> Result<U, Er>,
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
                .map(scalar_map)
                .collect::<Result<Vec<_>, Er>>()?,
        })
    }

    fn map_ref_mut<U, V>(
        &mut self,
        scalar_map: impl FnMut(&mut S) -> U,
        tensor_map: impl FnMut(&mut T) -> V,
    ) -> NetworkStore<V, U> {
        NetworkStore {
            tensors: self.tensors.iter_mut().map(tensor_map).collect(),
            scalar: self.scalar.iter_mut().map(scalar_map).collect(),
        }
    }

    fn map_ref_mut_enumerate<U, V>(
        &mut self,
        scalar_map: impl FnMut((usize, &mut S)) -> U,
        tensor_map: impl FnMut((usize, &mut T)) -> V,
    ) -> NetworkStore<V, U> {
        NetworkStore {
            tensors: self
                .tensors
                .iter_mut()
                .enumerate()
                .map(tensor_map)
                .collect(),
            scalar: self.scalar.iter_mut().enumerate().map(scalar_map).collect(),
        }
    }

    fn map_ref_mut_result<U, V, Er>(
        &mut self,
        scalar_map: impl FnMut(&mut S) -> Result<U, Er>,
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
                .map(scalar_map)
                .collect::<Result<Vec<_>, Er>>()?,
        })
    }

    fn map_ref_mut_result_enumerate<U, V, Er>(
        &mut self,
        scalar_map: impl FnMut((usize, &mut S)) -> Result<U, Er>,
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
                .map(scalar_map)
                .collect::<Result<Vec<_>, Er>>()?,
        })
    }
}
