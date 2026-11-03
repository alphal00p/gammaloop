#[macro_export]
macro_rules! define_index {
    (
$(#[$idx_meta:meta])*
$idx_vis:vis struct $Idx:ident ;) => {
    $(#[$idx_meta])*
    #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash,serde::Serialize, serde::Deserialize,bincode::Encode, bincode::Decode)]
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
};
}
