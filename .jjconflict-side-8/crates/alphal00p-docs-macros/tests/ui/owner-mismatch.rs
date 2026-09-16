use alphal00p_docs::func;

struct Actual;
struct Declared;

impl Actual {
    #[func(owner = "Declared")]
    fn create() -> Self {
        Self
    }
}

fn main() {}
