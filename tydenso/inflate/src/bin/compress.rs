use std::env;
use std::fs;
use std::io;
use std::path::PathBuf;

use miniz_oxide::deflate::compress_to_vec_zlib;
use miniz_oxide::inflate::decompress_to_vec_zlib;

fn required_path(args: &mut impl Iterator<Item = String>, label: &str) -> io::Result<PathBuf> {
    args.next().map(PathBuf::from).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("missing {label}; expected INPUT OUTPUT"),
        )
    })
}

fn main() -> io::Result<()> {
    let mut args = env::args();
    let _program = args.next();
    let input_path = required_path(&mut args, "input path")?;
    let output_path = required_path(&mut args, "output path")?;
    if args.next().is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "too many arguments; expected INPUT OUTPUT",
        ));
    }

    let input = fs::read(&input_path)?;
    let compressed = compress_to_vec_zlib(&input, 10);
    let round_trip = decompress_to_vec_zlib(&compressed).map_err(|error| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("could not verify compressed output: {error}"),
        )
    })?;
    if round_trip != input {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "compressed output did not round-trip to the input",
        ));
    }
    fs::write(&output_path, &compressed)?;
    eprintln!(
        "compressed {} bytes to {} bytes ({:.1}%)",
        input.len(),
        compressed.len(),
        compressed.len() as f64 * 100.0 / input.len() as f64,
    );
    Ok(())
}
