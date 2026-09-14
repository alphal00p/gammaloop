//! Independent small import check: an unchanged outer function still needs
//! its argument symbols remapped when their registration order changes.
use std::{
    error::Error,
    fs::File,
    io::{BufReader, BufWriter, Write},
};
use symbolica::{atom::Atom, function, symbol};

fn main() -> Result<(), Box<dyn Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    assert_eq!(args.len(), 2, "usage: import_remap write|read FILE");
    if args[0] == "write" {
        let x = Atom::var(symbol!("import_mre::x"));
        let y = Atom::var(symbol!("import_mre::y"));
        let expression = function!(symbol!("import_mre::f"), &x, &y);
        let mut writer = BufWriter::new(File::create_new(&args[1])?);
        expression.as_view().export(&mut writer)?;
        writer.flush()?;
        println!("exported {expression}");
    } else {
        assert_eq!(args[0], "read");
        let y = Atom::var(symbol!("import_mre::y"));
        let x = Atom::var(symbol!("import_mre::x"));
        let expected = function!(symbol!("import_mre::f"), &x, &y);
        let imported = Atom::import(&mut BufReader::new(File::open(&args[1])?), None)?;
        println!("imported {imported}; expected {expected}");
        assert_eq!(
            imported, expected,
            "Argument identities changed during import"
        );
    }
    Ok(())
}
