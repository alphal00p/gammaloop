#[derive(Debug, Clone)]
struct A {}
#[derive(Debug, Clone)]
struct B {}
pub trait Upgradable {}

impl Upgradable for A {}
impl Upgradable for B {}

trait SmallestUpgrade<T> {
    type Output;
    fn upgrade(self) -> Self::Output;
}

// impl SmallestUpgrade<A> for B {
//     type Output = B;
// }

impl<T: Upgradable> SmallestUpgrade<T> for T {
    type Output = T;
    fn upgrade(self) -> Self::Output {
        println!("Upgrading to self");
        self
    }
}

impl SmallestUpgrade<B> for A {
    type Output = B;
    fn upgrade(self) -> Self::Output {
        println!("Upgrading A and B to B");
        B {}
    }
}

impl SmallestUpgrade<A> for B {
    type Output = B;
    fn upgrade(self) -> Self::Output {
        println!("Upgrading B and A to B");

        self
    }
}

#[derive(Debug, Clone)]
struct Dense<T> {
    data: T,
}
#[derive(Debug, Clone)]
struct Sparse<T> {
    data: T,
}

// #[derive(Debug,Clone)]
// enum Upgrading {
//     A(A),
//     B(B)
// }

// impl SmallestUpgrade<Upgrading> for Upgrading {
//     type Output = Upgrading;
//     fn upgrade(self) -> Self::Output {
//         match self {
//             Upgrading::A(A)=> Upgrading::
//         }
//     }
// }

#[derive(Debug, Clone)]
enum NumTensor<T> {
    Dense(Dense<T>),
    Sparse(Sparse<T>),
}

impl<T> NumTensor<T> {
    pub fn new_dense(data: T) -> Self {
        NumTensor::Dense(Dense { data })
    }

    pub fn new_sparse(data: T) -> Self {
        NumTensor::Sparse(Sparse { data })
    }
}

impl<T, U> Contract<NumTensor<T>> for NumTensor<U>
where
    T: Clone,
    U: SmallestUpgrade<T> + Clone,
{
    type Output = NumTensor<U::Output>;
    fn contract(&self, other: &NumTensor<T>) -> Option<NumTensor<U::Output>> {
        match (self, other) {
            (NumTensor::Dense(s), NumTensor::Dense(o)) => Some(NumTensor::Dense(s.contract(o)?)),
            (NumTensor::Dense(s), NumTensor::Sparse(o)) => Some(NumTensor::Dense(s.contract(o)?)),
            (NumTensor::Sparse(s), NumTensor::Dense(o)) => Some(NumTensor::Dense(s.contract(o)?)),
            (NumTensor::Sparse(s), NumTensor::Sparse(o)) => Some(NumTensor::Sparse(s.contract(o)?)),
        }
    }
}

trait Contract<T> {
    type Output;
    fn contract(&self, rhs: &T) -> Option<Self::Output>;
}

impl<T: Clone, U: SmallestUpgrade<T> + Clone> Contract<Dense<T>> for Dense<U> {
    type Output = Dense<U::Output>;
    fn contract(&self, _rhs: &Dense<T>) -> Option<Self::Output> {
        println!("Contracting Dense Dense");
        Some(Dense {
            data: self.data.clone().upgrade(),
        })
    }
}

impl<T: Clone, U: SmallestUpgrade<T> + Clone> Contract<Dense<T>> for Sparse<U> {
    type Output = Dense<U::Output>;
    fn contract(&self, _rhs: &Dense<T>) -> Option<Self::Output> {
        println!("Contracting Sparse Dense");
        Some(Dense {
            data: self.data.clone().upgrade(),
        })
    }
}

impl<T: Clone, U: SmallestUpgrade<T> + Clone> Contract<Sparse<T>> for Sparse<U> {
    type Output = Sparse<U::Output>;
    fn contract(&self, _rhs: &Sparse<T>) -> Option<Self::Output> {
        println!("Contracting Sparse Sparse");
        Some(Sparse {
            data: self.data.clone().upgrade(),
        })
    }
}

impl<T: Clone, U: SmallestUpgrade<T> + Clone> Contract<Sparse<T>> for Dense<U> {
    type Output = Dense<U::Output>;
    fn contract(&self, _rhs: &Sparse<T>) -> Option<Self::Output> {
        println!("Contracting Dense Sparse");
        Some(Dense {
            data: self.data.clone().upgrade(),
        })
    }
}

// #[derive(Debug, Clone)]
// enum NumTensors {
//     DenseA(Dense<A>),
//     DenseB(Dense<B>),
//     SparseA(Sparse<A>),
//     SparseB(Sparse<B>),
// }

// impl Contract<NumTensors> for NumTensors {
//     type Output = NumTensors;
//     fn contract(&self, rhs: &NumTensors) -> Option<Self::Output> {
//         match (self, rhs) {
//             (NumTensors::DenseA(s), NumTensors::DenseA(o)) => {
//                 Some(NumTensors::DenseA(s.contract(o)?))
//             }
//             (NumTensors::DenseA(s), NumTensors::SparseA(o)) => {
//                 Some(NumTensors::DenseA(s.contract(o)?))
//             }
//             (NumTensors::DenseA(s), NumTensors::DenseB(o)) => {
//                 Some(NumTensors::DenseB(s.contract(o)?))
//             }
//             (NumTensors::DenseA(s), NumTensors::SparseB(o)) => {
//                 Some(NumTensors::DenseB(s.contract(o)?))
//             }
//             (NumTensors::DenseB(s), NumTensors::DenseA(o)) => {
//                 Some(NumTensors::DenseB(s.contract(o)?))
//             }
//             (NumTensors::DenseB(s), NumTensors::SparseA(o)) => {
//                 Some(NumTensors::DenseB(s.contract(o)?))
//             }
//             (NumTensors::DenseB(s), NumTensors::DenseB(o)) => {
//                 Some(NumTensors::DenseB(s.contract(o)?))
//             }
//             (NumTensors::DenseB(s), NumTensors::SparseB(o)) => {
//                 Some(NumTensors::DenseB(s.contract(o)?))
//             }
//             (NumTensors::SparseA(s), NumTensors::DenseA(o)) => {
//                 Some(NumTensors::DenseA(s.contract(o)?))
//             }
//             (NumTensors::SparseA(s), NumTensors::SparseA(o)) => {
//                 Some(NumTensors::SparseA(s.contract(o)?))
//             }
//             (NumTensors::SparseA(s), NumTensors::DenseB(o)) => {
//                 Some(NumTensors::DenseB(s.contract(o)?))
//             }
//             (NumTensors::SparseA(s), NumTensors::SparseB(o)) => {
//                 Some(NumTensors::SparseB(s.contract(o)?))
//             }
//             (NumTensors::SparseB(s), NumTensors::DenseA(o)) => {
//                 Some(NumTensors::DenseB(s.contract(o)?))
//             }
//             (NumTensors::SparseB(s), NumTensors::SparseA(o)) => {
//                 Some(NumTensors::SparseB(s.contract(o)?))
//             }
//             (NumTensors::SparseB(s), NumTensors::DenseB(o)) => {
//                 Some(NumTensors::DenseB(s.contract(o)?))
//             }
//             (NumTensors::SparseB(s), NumTensors::SparseB(o)) => {
//                 Some(NumTensors::SparseB(s.contract(o)?))
//             }
//         }
//     }
// }

#[derive(Debug, Clone)]
enum NumTensors {
    A(NumTensor<A>),
    B(NumTensor<B>),
}

trait NumTensorable<T> {
    fn new_dense(data: T) -> NumTensors;
    fn new_sparse(data: T) -> NumTensors;
}

impl NumTensorable<A> for NumTensors {
    fn new_dense(data: A) -> NumTensors {
        NumTensors::A(NumTensor::new_dense(data))
    }
    fn new_sparse(data: A) -> NumTensors {
        NumTensors::A(NumTensor::new_sparse(data))
    }
}

impl NumTensorable<B> for NumTensors {
    fn new_dense(data: B) -> NumTensors {
        NumTensors::B(NumTensor::new_dense(data))
    }
    fn new_sparse(data: B) -> NumTensors {
        NumTensors::B(NumTensor::new_sparse(data))
    }
}

impl Contract<NumTensors> for NumTensors {
    type Output = NumTensors;
    fn contract(&self, rhs: &NumTensors) -> Option<Self::Output> {
        match (self, rhs) {
            (NumTensors::A(s), NumTensors::A(o)) => Some(NumTensors::A(s.contract(o)?)),
            (NumTensors::A(s), NumTensors::B(o)) => Some(NumTensors::B(s.contract(o)?)),
            (NumTensors::B(s), NumTensors::A(o)) => Some(NumTensors::B(s.contract(o)?)),
            (NumTensors::B(s), NumTensors::B(o)) => Some(NumTensors::B(s.contract(o)?)),
        }
    }
}

fn main() {
    let dense_a = Dense { data: A {} };
    let dense_b = Dense { data: B {} };
    let sparse_a = Sparse { data: A {} };
    let sparse_b = Sparse { data: B {} };
    // let dense_b = NumTensors::DenseB(Dense { data: B {} });

    // let dense_a = NumTensors::DenseA(Dense { data: A {} });
    // let sparse_a = NumTensors::SparseA(Sparse { data: A {} });
    // let sparse_b = NumTensors::SparseB(Sparse { data: B {} });

    // // let a = dense_b.contract(&sparse_a);
    // // let b = dense_b.contract(&dense_b);
    // // let c = sparse_a.contract(&sparse_b);
    // // let d = sparse_a.contract(&dense_b);

    // let tensor_net = vec![dense_b.clone(), sparse_b, dense_a, sparse_a];

    // println!("{:?}", tensor_net);

    // tensor_net.iter().try_fold(dense_b, |x, y| x.contract(y));

    let dense_a = NumTensors::new_dense(A {});
    let dense_b = NumTensors::new_dense(B {});
    let sparse_a = NumTensors::new_sparse(A {});
    let sparse_b = NumTensors::new_sparse(B {});

    let tensor_net2 = vec![dense_a.clone(), sparse_b, dense_b, sparse_a];

    println!("{:?}", tensor_net2);
    tensor_net2.iter().try_fold(dense_a, |x, y| x.contract(y));
}
