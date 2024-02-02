use std::{collections::HashMap, time::Instant};

use symbolica::{
    representations::{default::Linear, Atom, Identifier},
    state::{State, Workspace},
};

fn dump_c_with_func(_levels: Vec<Vec<(Identifier, Vec<Atom>)>>, _params: Vec<Atom>) {}
fn dump_c(_levels: Vec<Vec<HashMap<Atom, Atom>>>, _params: Vec<Atom>) {}

fn main() {
    let mut state = State::new();
    let ws: Workspace<Linear> = Workspace::new();
    let from_filemap: Vec<Vec<HashMap<String, String>>> =
        serde_yaml::from_reader(std::fs::File::open("outmap.yaml").unwrap()).unwrap();

    let paramstr: Vec<String> =
        serde_yaml::from_reader(std::fs::File::open("params.yaml").unwrap()).unwrap();

    let params: Vec<Atom> = paramstr
        .iter()
        .map(|x| Atom::parse(x, &mut state, &ws).unwrap())
        .collect();

    let startmap = Instant::now();

    let levelsmap = from_filemap
        .iter()
        .map(|x| {
            x.iter()
                .map(|x| {
                    let mut a = HashMap::new();
                    for (k, v) in x.iter() {
                        a.insert(
                            Atom::parse(k, &mut state, &ws).unwrap(),
                            Atom::parse(v, &mut state, &ws).unwrap(),
                        );
                    }
                    a
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let durationmap = startmap.elapsed();
    let from_file: Vec<Vec<(String, Vec<String>)>> =
        serde_yaml::from_reader(std::fs::File::open("out.yaml").unwrap()).unwrap();

    let start = Instant::now();
    let levels: Vec<Vec<(Identifier, Vec<Atom>)>> = from_file
        .iter()
        .map(|x| {
            x.iter()
                .map(|(s, v)| {
                    (
                        state.get_or_insert_fn(s, None).unwrap(),
                        v.iter()
                            .map(|x| Atom::parse(x, &mut state, &ws).unwrap())
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let duration = start.elapsed();

    println!("Time to parse map: {:?}", durationmap);
    println!("Time to parse vec: {:?}", duration);

    dump_c(levelsmap, params.clone());
    dump_c_with_func(levels, params);
}
