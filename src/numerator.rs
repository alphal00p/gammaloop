use crate::{
    graph::{EdgeType, Graph, LoopMomentumBasis},
    model::Model,
    momentum::FourMomentum,
    utils::{FloatLike, F},
};
use ahash::AHashMap;
use itertools::Itertools;
use log::debug;
use serde::{ser::SerializeStruct, Deserialize, Serialize};
use spenso::{
    AbstractIndex, NamedStructure, Representation, Shadowable, SymbolicTensor, TensorNetwork,
    TensorStructure,
};
use spenso::{Complex, ParamTensor};
use symbolica::{atom::AtomView, domains::float::Complex as SymComplex};
use symbolica::{
    atom::{Atom, FunctionBuilder, Symbol},
    state::State,
    printer::{AtomPrinter, PrintOptions}
};

pub fn apply_replacements(
    graph: &Graph,
    model: &Model,
    lmb: &LoopMomentumBasis,
    mut atom: Atom,
) -> Atom {
    atom = model.substitute_model_params(&atom);

    for edge in &graph.edges {
        atom = edge.substitute_lmb(atom, graph, lmb);
    }
    atom
}

#[derive(Debug, Clone)]
#[allow(clippy::type_complexity)]
pub struct Numerator {
    pub expression: Atom,
    pub network: Option<TensorNetwork<ParamTensor<NamedStructure<Symbol, Vec<Atom>>>, Atom>>,
    pub const_map: AHashMap<Atom, Complex<F<f64>>>,
}

impl Serialize for Numerator {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let expression =
            AtomPrinter::new_with_options(self.expression.as_view(), PrintOptions::file())
                .to_string();

        let const_map: AHashMap<String, Complex<F<f64>>> = self
            .const_map
            .iter()
            .map(|(k, &v)| {
                (
                    AtomPrinter::new_with_options(k.as_view(), PrintOptions::file()).to_string(),
                    v,
                )
            })
            .collect();

        let mut state = serializer.serialize_struct("Numerator", 3)?;
        state.serialize_field("expression", &expression)?;
        state.serialize_field("const_map", &const_map)?;
        state.end()
    }
}

impl<'de> Deserialize<'de> for Numerator {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct NumeratorData {
            expression: String,
            const_map: AHashMap<String, Complex<F<f64>>>,
        }

        let data = NumeratorData::deserialize(deserializer)?;

        let expression = Atom::parse(&data.expression).map_err(serde::de::Error::custom)?;

        let const_map: AHashMap<Atom, Complex<F<f64>>> = data
            .const_map
            .into_iter()
            .map(|(k, v)| (Atom::parse(&k).unwrap(), v))
            .collect();

        let sym_tensor: SymbolicTensor = expression
            .clone()
            .try_into()
            .map_err(serde::de::Error::custom)?;
        let network: TensorNetwork<
            spenso::ParamTensor<NamedStructure<symbolica::atom::Symbol, Vec<Atom>>>,
            Atom,
        > = sym_tensor
            .to_network()
            .map_err(serde::de::Error::custom)?
            .to_fully_parametric();

        Ok(Numerator {
            expression,
            network: Some(network),
            const_map,
        })
    }
}

impl Numerator {
    pub fn substitute_model_params(&mut self, model: &Model) {
        self.expression = model.substitute_model_params(&self.expression);
    }

    pub fn evaluate<T: FloatLike>(
        &self,
        emr: &[FourMomentum<F<T>>],
        graph: &Graph,
    ) -> Complex<F<T>> {
        let emr_params = graph
            .edges
            .iter()
            .enumerate()
            .map(|(i, _)| {
                let named_structur: NamedStructure<&str> = NamedStructure::from_iter(
                    [(AbstractIndex(i), Representation::Lorentz(4.into()))],
                    "Q",
                    Some(i),
                );
                debug!("Q{}", i);
                named_structur.to_shell().shadow().unwrap()
            })
            .collect_vec();

        let mut constmap: AHashMap<AtomView<'_>, SymComplex<F<T>>> = AHashMap::new();

        let mut pols = Vec::new();
        for (i, ext) in graph
            .edges
            .iter()
            .enumerate()
            .filter(|(_, e)| !matches!(e.edge_type, EdgeType::Virtual))
        {
            match ext.edge_type {
                EdgeType::Incoming => {
                    for (a, p) in ext
                        .particle
                        .incoming_polarization_match(i, &emr[i])
                        .into_iter()
                    {
                        pols.push((a, p));
                    }
                }

                EdgeType::Outgoing => {
                    for (a, p) in ext
                        .particle
                        .outgoing_polarization_match(i, &emr[i])
                        .into_iter()
                    {
                        pols.push((a, p));
                    }
                }
                _ => {}
            }
        }

        for (a, p) in &pols {
            constmap.insert(a.as_view(), p.clone().into());
        }

        for (key, value) in self.const_map.iter() {
            let v = Complex::new(F::<T>::from_ff64(value.re), F::<T>::from_ff64(value.im));
            constmap.insert(key.as_view(), v.into());
        }

        for (i, pe) in emr.iter().enumerate() {
            constmap.insert(
                emr_params[i][1.into()].as_view(),
                Complex::new_re(pe.spatial.px.clone()).into(),
            );

            constmap.insert(
                emr_params[i][2.into()].as_view(),
                Complex::new_re(pe.spatial.py.clone()).into(),
            );

            constmap.insert(
                emr_params[i][3.into()].as_view(),
                Complex::new_re(pe.spatial.pz.clone()).into(),
            );

            constmap.insert(
                emr_params[i][0.into()].as_view(),
                Complex::new_re(pe.temporal.value.clone()).into(),
            );
        }

        // let mut file = File::create("serializable_map.json").unwrap();

        // debug!("Serializable map: {}", serialized);
        // file.write_all(serialized.as_bytes()).unwrap();

        // println!(
        //     "Numerator evaluated: {}",
        //     numerator_network.scalar.as_ref().unwrap()
        // );

        let serializable_map: AHashMap<String, Complex<F<T>>> = constmap
            .iter()
            .map(|(k, v)| (k.to_string(), v.clone().into()))
            .collect();

        let serialized = serde_json::to_string(&serializable_map).unwrap();
        debug!("Serializable map: {:?}", serialized);

        let mut evaluated = self
            .network
            .as_ref()
            .unwrap()
            .evaluate::<SymComplex<F<T>>, _>(|c| c.into(), &constmap);

        evaluated.contract();

        evaluated.result().unwrap().into()
    }

    pub fn process(&mut self, graph: &Graph, model: &Model) {
        debug!("processing numerator for graph: {}", graph.name);

        let mut const_map = AHashMap::new();
        model.append_parameter_map(&mut const_map);

        let i_float = Complex::new_i();
        let i = Atom::new_var(State::I);
        const_map.insert(i, i_float);
        self.const_map.extend(const_map);

        self.substitute_model_params(model);

        let sym_tensor: SymbolicTensor = self.expression.clone().try_into().unwrap();

        let network = sym_tensor.to_network().unwrap().to_fully_parametric();
        self.network = Some(network);
    }

    pub fn generate(graph: &Graph) -> Self {
        debug!("generating numerator for graph: {}", graph.name);

        let vatoms: Vec<Atom> = graph
            .vertices
            .iter()
            .flat_map(|v| v.contracted_vertex_rule(graph))
            .collect();
        // graph
        //     .edges
        //     .iter()
        //     .filter(|e| e.edge_type == EdgeType::Virtual)
        //     .map(|e| e.particle);
        let eatoms: Vec<Atom> = graph.edges.iter().map(|e| e.numerator(graph)).collect();
        let mut builder = Atom::new_num(1);

        for v in &vatoms {
            builder = builder * v;
        }

        for e in &eatoms {
            builder = builder * e;
        }

        let i = Atom::new_var(State::I);
        let a = Atom::new_var(State::get_symbol("a_"));
        let b = Atom::new_var(State::get_symbol("b_"));

        let complex = FunctionBuilder::new(State::get_symbol("complex"))
            .add_arg(&a)
            .add_arg(&b)
            .finish();

        builder = complex.into_pattern().replace_all(
            builder.as_view(),
            &(&a + &b * &i).into_pattern(),
            None,
            None,
        );

        debug!("numerator: {}", builder);

        let expression = builder.clone();

        Numerator {
            expression,
            network: None,
            const_map: AHashMap::new(),
        }
    }
}

#[cfg(test)]
mod tests;
