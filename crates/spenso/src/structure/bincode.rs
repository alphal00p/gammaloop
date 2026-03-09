use crate::structure::IndexLess;
use crate::structure::RepName;

impl<T: RepName, __Context> ::bincode::Decode<__Context> for IndexLess<T>
where
    T: ::bincode::Decode<__Context>,
{
    fn decode<__D: ::bincode::de::Decoder<Context = __Context>>(
        decoder: &mut __D,
    ) -> core::result::Result<Self, ::bincode::error::DecodeError> {
        core::result::Result::Ok(Self {
            structure: ::bincode::Decode::decode(decoder)?,
        })
    }
}
impl<'__de, T: RepName, __Context> ::bincode::BorrowDecode<'__de, __Context> for IndexLess<T>
where
    T: ::bincode::de::BorrowDecode<'__de, __Context>,
{
    fn borrow_decode<__D: ::bincode::de::BorrowDecoder<'__de, Context = __Context>>(
        decoder: &mut __D,
    ) -> core::result::Result<Self, ::bincode::error::DecodeError> {
        core::result::Result::Ok(Self {
            structure: ::bincode::BorrowDecode::<'_, __Context>::borrow_decode(decoder)?,
        })
    }
}
