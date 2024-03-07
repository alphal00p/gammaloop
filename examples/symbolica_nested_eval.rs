use std::fs;

use symbolica::{
    poly::evaluate::ExpressionEvaluator,
    representations::{Atom, Symbol},
    state::{State, Workspace},
};

fn main() {
    let mut state = State::get_global_state().write().unwrap();
    let _ws = Workspace::new();
    // let from_filemap: Vec<Vec<HashMap<String, String>>> =
    //     serde_yaml::from_reader(std::fs::File::open("outmap.yaml").unwrap()).unwrap();

    let from_file: Vec<Vec<(String, Vec<String>)>> =
        serde_yaml::from_reader(std::fs::File::open("out.yaml").unwrap()).unwrap();

    let levels: Vec<Vec<(Symbol, Vec<Atom>)>> = from_file
        .iter()
        .map(|x| {
            x.iter()
                .map(|(s, v)| {
                    (
                        state.get_or_insert_fn(s, None).unwrap(),
                        v.iter()
                            .map(|x| {
                                let a = Atom::parse(x, &mut state).unwrap();
                                let _a_exp = Atom::new();
                                a.as_view().expand();
                                a
                            })
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let e = ExpressionEvaluator::new(levels, 10);

    fs::write("eval.cpp", format!("{e}")).unwrap();
}
