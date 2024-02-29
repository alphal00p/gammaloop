use std::{fs};

use symbolica::{
    poly::evaluate::ExpressionEvaluator,
    representations::{default::Linear, Atom, Identifier},
    state::{ResettableBuffer, State, Workspace},
};

fn main() {
    let mut state = State::new();
    let ws: Workspace<Linear> = Workspace::new();
    // let from_filemap: Vec<Vec<HashMap<String, String>>> =
    //     serde_yaml::from_reader(std::fs::File::open("outmap.yaml").unwrap()).unwrap();

    let from_file: Vec<Vec<(String, Vec<String>)>> =
        serde_yaml::from_reader(std::fs::File::open("out.yaml").unwrap()).unwrap();

    let levels: Vec<Vec<(Identifier, Vec<Atom>)>> = from_file
        .iter()
        .map(|x| {
            x.iter()
                .map(|(s, v)| {
                    (
                        state.get_or_insert_fn(s, None).unwrap(),
                        v.iter()
                            .map(|x| {
                                let a = Atom::parse(x, &mut state, &ws).unwrap();
                                let mut a_exp = Atom::new();
                                a.as_view().expand(&ws, &state, &mut a_exp);
                                a_exp
                            })
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let e = ExpressionEvaluator::new(levels, &state, 10);

    fs::write("eval.cpp", format!("{e}"));
}
