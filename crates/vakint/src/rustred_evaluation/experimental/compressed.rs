//! Bounded transport for trusted generated native assets; not a CAS codec.

use std::io::{Read, Write};

use flate2::{Compression, GzBuilder, bufread::GzDecoder};

pub(super) fn encode(bytes: &[u8]) -> Result<Vec<u8>, String> {
    let mut encoder = GzBuilder::new()
        .mtime(0)
        .write(Vec::new(), Compression::best());
    encoder
        .write_all(bytes)
        .map_err(|error| error.to_string())?;
    encoder.finish().map_err(|error| error.to_string())
}

pub(crate) fn decode(input: &[u8], max_bytes: usize) -> Result<Vec<u8>, String> {
    if input.len() > max_bytes {
        return Err("compressed native asset exceeds byte limit".into());
    }
    let bound = u64::try_from(max_bytes)
        .ok()
        .and_then(|value| value.checked_add(1))
        .ok_or("native decompression limit overflow")?;
    let mut decoder = GzDecoder::new(input);
    let mut bytes = Vec::new();
    decoder
        .by_ref()
        .take(bound)
        .read_to_end(&mut bytes)
        .map_err(|error| format!("decompress native asset: {error}"))?;
    if bytes.len() > max_bytes {
        return Err("decompressed native asset exceeds byte limit".into());
    }
    if !decoder.into_inner().is_empty() {
        return Err("trailing bytes after compressed native asset".into());
    }
    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_single_member_transport_preserves_bytes() {
        let input = vec![42; 1024];
        let encoded = encode(&input).unwrap();
        assert_eq!(encoded, encode(&input).unwrap());
        assert_eq!(decode(&encoded, input.len()).unwrap(), input);
        assert!(
            decode(&encoded, input.len() - 1)
                .unwrap_err()
                .contains("byte limit")
        );
        assert!(
            decode(&encoded, encoded.len() - 1)
                .unwrap_err()
                .contains("byte limit")
        );
    }

    #[test]
    fn damaged_truncated_or_trailing_members_fail_closed() {
        let encoded = encode(&vec![7; 1024]).unwrap();
        for end in 0..encoded.len() {
            assert!(decode(&encoded[..end], 2048).is_err(), "truncated at {end}");
        }
        let mut corrupt = encoded.clone();
        let checksum = corrupt.len() - 8;
        corrupt[checksum] ^= 1;
        assert!(decode(&corrupt, 2048).is_err());
        let mut trailing = encoded.clone();
        trailing.push(0);
        assert!(decode(&trailing, 2048).unwrap_err().contains("trailing"));
        let mut concatenated = encoded.clone();
        concatenated.extend_from_slice(&encoded);
        assert!(
            decode(&concatenated, 2048)
                .unwrap_err()
                .contains("trailing")
        );
    }
}
