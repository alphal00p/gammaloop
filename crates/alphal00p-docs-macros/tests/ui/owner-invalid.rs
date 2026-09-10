use alphal00p_docs::func;

struct Actual;

impl Actual {
    #[func(owner = "Vec<")]
    fn create() -> Self {
        Self
    }
}

fn main() {}
