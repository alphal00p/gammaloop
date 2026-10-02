use symbolica::atom::Atom;

use crate::network::StructureLessDisplay;

impl StructureLessDisplay for Vec<Atom> {
    fn display(&self) -> String {
        self.iter()
            .map(|t| t.to_string())
            .collect::<Vec<String>>()
            .join(", ")
    }
}
