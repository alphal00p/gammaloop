use std::hash::Hash;
pub trait IndexLike:
    Copy + PartialEq + Eq + PartialOrd + Ord + Hash + From<usize> + Into<usize>
{
}

impl<ID: Copy + PartialEq + Eq + PartialOrd + Ord + Hash + From<usize> + Into<usize>> IndexLike
    for ID
{
}

#[macro_export]
macro_rules! define_indexed_vec {
    (
        $(#[$idx_meta:meta])*
        $idx_vis:vis struct $Idx:ident ;

        $(#[$vec_meta:meta])*
        $vec_vis:vis struct $Vec:ident ;
    ) => {
        /* ——————————————————— index new‑type ——————————————————— */

        $(#[$idx_meta])*
        #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
        #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
        #[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
        $idx_vis struct $Idx(pub usize);

        impl ::std::convert::From<usize> for $Idx {
            fn from(value: usize) -> Self {
                $Idx(value)
            }
        }



        impl ::std::ops::Add<$Idx> for $Idx {
            type Output = $Idx;

            fn add(self, other: $Idx) -> Self::Output {
                $Idx(self.0 + other.0)
            }
        }

        impl ::std::ops::Sub<$Idx> for $Idx {
            type Output = $Idx;

            fn sub(self, other: $Idx) -> Self::Output {
                $Idx(self.0 - other.0)
            }
        }

        impl ::std::ops::AddAssign<$Idx> for $Idx {
            fn add_assign(&mut self, other: $Idx) {
                self.0 += other.0;
            }
        }

        impl ::std::ops::SubAssign<$Idx> for $Idx {
            fn sub_assign(&mut self, other: $Idx) {
                self.0 -= other.0;
            }
        }

        impl ::std::convert::From<$Idx> for usize {
            fn from(value: $Idx) -> Self {
                value.0
            }
        }





        /* ——————————————————— vector new‑type ——————————————————— */

        $(#[$vec_meta])*
        #[derive(Clone, Debug, Default, Hash, PartialEq, Eq,PartialOrd,Ord)]
        #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
        #[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
        $vec_vis struct $Vec<T>(::std::vec::Vec<T>);

        /* --- Restricted indexing -------------------------------------------------- */

        impl<T> ::std::ops::Index<$Idx> for $Vec<T> {
            type Output = T;
            #[inline] fn index(&self, i: $Idx) -> &Self::Output { &self.0[i.0] }
        }
        impl<T> ::std::ops::IndexMut<$Idx> for $Vec<T> {
            #[inline] fn index_mut(&mut self, i: $Idx) -> &mut Self::Output { &mut self.0[i.0] }
        }

        impl<T>  ::std::convert::AsMut<[T]> for $Vec<T>{
            fn as_mut(&mut self)->&mut [T]{
                self.0.as_mut()
            }
        }

        impl<T> $crate::half_edge::swap::Swap<$Idx> for $Vec<T>{
            fn swap(&mut self, a: $Idx, b: $Idx) {
                self.0.swap(a.0, b.0);
            }

            // type Item = T;

            // fn filter(&self, id: &$Idx, filter: &impl Fn(&$Idx, &Self::Item) -> bool) -> bool {
            //     filter(id,&self[*id])
            // }

            fn len(&self) -> $Idx {
                $Idx(self.0.len())
            }
            #[inline] fn is_zero_length(&self) -> bool { self.0.is_empty() }

        }

        /* --- Delegated Vec<T> API ------------------------------------------------- */

        ::paste::paste! {
                    #[derive(Clone)]
                    $vec_vis struct [<$Vec Iter>]<'a, T>(std::iter::Map<
                        std::iter::Enumerate<std::slice::Iter<'a, T>>,
                        fn((usize, &T)) -> ($Idx, &T),
                    >);
                    impl<'a, T> ::std::iter::Iterator for [<$Vec Iter>]<'a, T> {
                        type Item = ($Idx, &'a T);
                        #[inline] fn next(&mut self) -> Option<Self::Item> { self.0.next() }
                        #[inline] fn size_hint(&self) -> (usize, Option<usize>) { self.0.size_hint() }
                    }
                    impl<'a, T> ::std::iter::DoubleEndedIterator for [<$Vec Iter>]<'a, T> {
                        #[inline] fn next_back(&mut self) -> Option<Self::Item> { self.0.next_back() }
                    }
                    impl<'a, T> ::std::iter::ExactSizeIterator for [<$Vec Iter>]<'a, T> {}
                    impl<'a, T> ::std::iter::FusedIterator     for [<$Vec Iter>]<'a, T> {}


                    impl<'a, T> ::std::iter::IntoIterator for &'a $Vec<T> {
                        type Item = ($Idx, &'a T);
                        type IntoIter = [<$Vec Iter>]<'a, T>;
                        fn into_iter(self) -> Self::IntoIter {
                            [<$Vec Iter>](self.0.iter().enumerate().map(|(u, t)| ($Idx(u), t)))
                        }
                    }

                    $vec_vis struct [<$Vec IterMut>]<'a, T>(::std::iter::Map<
                        std::iter::Enumerate<std::slice::IterMut<'a, T>>,
                        fn((usize, &mut T)) -> ($Idx, &mut T),
                    >);
                    impl<'a, T> ::std::iter::Iterator for [<$Vec IterMut>]<'a, T> {
                        type Item = ($Idx, &'a mut T);
                        #[inline] fn next(&mut self) -> Option<Self::Item> { self.0.next() }
                        #[inline] fn size_hint(&self) -> (usize, Option<usize>) { self.0.size_hint() }
                    }
                    impl<'a, T> ::std::iter::DoubleEndedIterator for [<$Vec IterMut>]<'a, T> {
                        #[inline] fn next_back(&mut self) -> Option<Self::Item> { self.0.next_back() }
                    }
                    impl<'a, T> ::std::iter::ExactSizeIterator for [<$Vec IterMut>]<'a, T> {}
                    impl<'a, T> ::std::iter::FusedIterator     for [<$Vec IterMut>]<'a, T> {}

                    impl<'a, T> ::std::iter::IntoIterator for &'a mut $Vec<T> {
                        type Item = ($Idx, &'a mut T);
                        type IntoIter = [<$Vec IterMut>]<'a, T>;
                        fn into_iter(self) -> Self::IntoIter {
                            [<$Vec IterMut>](self.0.iter_mut().enumerate().map(|(u, t)| ($Idx(u), t)))
                        }
                    }



                    $vec_vis struct [<$Vec IntoIter>]<T>(::std::iter::Map<std::iter::Enumerate<std::vec::IntoIter<T>>, fn((usize, T)) -> ($Idx, T)>);
                    impl<T> ::std::iter::Iterator for [<$Vec IntoIter>]<T> {
                        type Item = ($Idx,T);
                        #[inline] fn next(&mut self) -> Option<Self::Item> { self.0.next() }
                        #[inline] fn size_hint(&self) -> (usize, Option<usize>) { self.0.size_hint() }
                    }
                    impl<T> ::std::iter::DoubleEndedIterator for [<$Vec IntoIter>]<T> {
                        #[inline] fn next_back(&mut self) -> Option<Self::Item> { self.0.next_back() }
                    }
                    impl<T> ::std::iter::ExactSizeIterator for [<$Vec IntoIter>]<T> {}
                    impl<T> ::std::iter::FusedIterator     for [<$Vec IntoIter>]<T> {}

                    impl<T> ::std::iter::IntoIterator for $Vec<T> {
                        type Item = ($Idx,T);
                        type IntoIter =  [<$Vec IntoIter>]<T>;
                        #[inline] fn into_iter(self) -> Self::IntoIter { [<$Vec IntoIter>](self.0.into_iter().enumerate().map(|(u, t)| ($Idx(u), t))) }
                    }


                    impl<T> $Vec<T>{
                        #[inline] pub fn iter<'a>(&'a self) -> [<$Vec Iter>]<'a, T> { self.into_iter() }
                        #[inline] pub fn iter_mut<'a>(&'a mut self) -> [<$Vec IterMut>]<'a, T> { self.into_iter() }

                    }

        }

        impl<T> $Vec<T> {

            pub fn write_display<W: ::std::fmt::Write>(&self, writer: &mut W,formater:impl Fn(&T)->String) -> ::std::fmt::Result {
                writer.write_str("[")?;

                for (i,item) in self {
                    if i.0 != 0{
                        write!(writer, ", ")?;
                    }
                    write!(writer, "{}", formater(item))?;
                }
                writer.write_str("]")?;
                Ok(())
            }

            pub fn display_string(&self,formatter:impl Fn(&T)->String) -> String {
                let mut result = String::new();
                self.write_display(&mut result, formatter).unwrap();
                result
            }



            /* construction */
            #[inline] pub fn new() -> Self { Self(::std::vec::Vec::new()) }
            #[inline] pub fn with_capacity(c: usize) -> Self { Self(::std::vec::Vec::with_capacity(c)) }

            /* capacity */
            // #[inline] pub fn len(&self) -> usize { self.0.len() }

            #[inline] pub fn capacity(&self) -> usize { self.0.capacity() }
            #[inline] pub fn reserve(&mut self, n: usize) { self.0.reserve(n) }
            #[inline] pub fn reserve_exact(&mut self, n: usize) { self.0.reserve_exact(n) }
            #[inline] pub fn shrink_to_fit(&mut self) { self.0.shrink_to_fit() }

            /* push / pop */
            #[inline] pub fn push(&mut self, value: T) { self.0.push(value) }
            #[inline] pub fn pop(&mut self) -> Option<T> { self.0.pop() }


            #[inline] pub fn swap(&mut self, a: $Idx, b: $Idx) {
                self.0.swap(a.0, b.0);
            }

            #[inline] pub fn split_off(&mut self, at: $Idx) -> Self {
                Self(self.0.split_off(at.0))
            }

            /* insertion / removal with the index new‑type */
            #[inline] pub fn insert(&mut self, idx: $Idx, v: T) { self.0.insert(idx.0, v) }
            #[inline] pub fn remove(&mut self, idx: $Idx) -> T { self.0.remove(idx.0) }
            #[inline] pub fn swap_remove(&mut self, idx: $Idx) -> T { self.0.swap_remove(idx.0) }

            /* get APIs using the index new‑type */
            #[inline] pub fn get(&self, idx: $Idx) -> Option<&T> { self.0.get(idx.0) }
            #[inline] pub fn get_mut(&mut self, idx: $Idx) -> Option<&mut T> { self.0.get_mut(idx.0) }

            /* iteration */

            /* miscellaneous */
            #[inline] pub fn clear(&mut self) { self.0.clear() }
            #[inline] pub fn truncate(&mut self, len: usize) { self.0.truncate(len) }

            /* fall‑back escape hatch – intentionally *not* public: */
            #[inline] pub fn raw(&self) -> &::std::vec::Vec<T> { &self.0 }
        }

        /* --- standard trait impls ------------------------------------------------- */

        impl<T> ::std::iter::FromIterator<($Idx,T)> for $Vec<T> {
            #[inline] fn from_iter<I: ::std::iter::IntoIterator<Item = ($Idx,T)>>(it: I) -> Self {
                Self(::std::vec::Vec::from_iter(it.into_iter().map(|(_, val)|  val)))
            }
        }impl<T> ::std::iter::FromIterator<T> for $Vec<T> {
            #[inline] fn from_iter<I: ::std::iter::IntoIterator<Item = T>>(it: I) -> Self {
                Self(::std::vec::Vec::from_iter(it))
            }
        }

        impl<T> ::std::convert::AsRef<[T]> for $Vec<T> {
            #[inline] fn as_ref(&self) -> &[T] { &self.0 }
        }



        impl<T> ::std::iter::Extend<($Idx,T)> for $Vec<T> {
            #[inline] fn extend<I: ::std::iter::IntoIterator<Item = ($Idx,T)>>(&mut self, it: I) {
                self.0.extend(it.into_iter().map(|(_, val)|  val));
            }
        }


        impl<T> ::std::convert::From<::std::vec::Vec<T>> for $Vec<T> {
            #[inline] fn from(v: ::std::vec::Vec<T>) -> Self { Self(v) }
        }




        /// Permutation constructors

        impl $Vec<Option<$Idx>>{
            pub fn fill_in(&mut self,contained:impl Fn(&$Idx)->bool){
                let mut new_shifted = $Idx(0);

                for (_, new_e) in self {
                    if new_e.is_none() {
                        while contained(&new_shifted) {
                            new_shifted.0 += 1;
                        }
                        *new_e = Some(new_shifted);
                        new_shifted.0 += 1;
                    }
                }
            }
        }

        impl ::std::convert::TryFrom<$Vec<Option<$Idx>>> for $Vec<$Idx>{
            type Error = ();
            fn try_from(vec: $Vec<Option<$Idx>>) -> Result<Self, Self::Error> {
                vec.into_iter().map(|(i,e)| e.ok_or(()).map(|e|(i,e))).collect()
            }
        }

        impl ::std::convert::TryFrom<$Vec<Option<$Idx>>> for $crate::permutation::Permutation{
            type Error = ();
            fn try_from(vec: $Vec<Option<$Idx>>) -> Result<Self, Self::Error> {
                let new_vec:Vec<usize> = vec.into_iter().map(|(i,e)| e.ok_or(()).map(|e|(i,e))).collect::<Result<$Vec<$Idx>, ()>>()?.into_iter().map(|(_,x)| usize::from(x)).collect();

                Ok($crate::permutation::Permutation::from_map(new_vec))
            }
        }
    };
}
