#[cfg(feature = "shadowing")]
#[allow(dead_code)]
#[derive(Debug, Clone, Encode, Decode)]
#[bincode(decode_context = "symbolica::state::StateMap")]
pub struct Levels<
    T: HasStructure<Structure = S>,
    S: TensorStructure + HasName + Clone,
    Sc: AtomCore,
> {
    pub levels: Vec<Network<ParamTensor<S>, Sc>>,
    pub initial: Network<T, Sc>,
    // fn_map: FunctionMap<'static, Complex<T>>,
    params: Vec<Sc>,
}

#[cfg(feature = "shadowing")]
impl<T, S, Sc: AtomCore> From<Network<T, Sc>> for Levels<T, S, Sc>
where
    T: HasStructure<Structure = S>,
    S: TensorStructure + HasName + Clone,
{
    fn from(t: Network<T, Sc>) -> Self {
        Levels {
            initial: t,
            levels: vec![],
            // fn_map: FunctionMap::new(),
            params: vec![],
        }
    }
}

#[cfg(feature = "shadowing")]
impl<
        T: Clone + RefZero + Display,
        S: TensorStructure + Clone + TracksCount + Display,
        Sc: AtomCore + Clone,
    > Levels<MixedTensor<T, S>, S, Sc>
where
    MixedTensor<T, S>: Contract<LCM = MixedTensor<T, S>> + Trace,

    S: HasName<Name: IntoSymbol, Args: IntoArgs> + ToSymbolic + StructureContract,
{
    fn contract_levels(
        &mut self,
        depth: usize,
        // fn_map: &mut FunctionMap<'a, Complex<T>>,
    ) {
        let mut not_done = true;
        let level = self.levels.len();

        if let Some(current_level) = self.levels.last_mut() {
            current_level.namesym(&format!("L{level}"))
        } else {
            not_done = false;
        }

        let nl = if let Some(current_level) = self.levels.last() {
            let mut new_level = current_level.shadow();
            new_level
                .contract_algo(|tn| tn.edge_to_min_degree_node_with_depth(depth))
                .unwrap();
            if new_level.graph.n_nodes() == 1 {
                not_done = false;
            }
            Some(new_level)
        } else {
            None
        };

        if let Some(nl) = nl {
            self.levels.push(nl);
        }

        if not_done {
            self.contract_levels(depth)
        }
    }

    pub fn contract<R>(&mut self, depth: usize, fn_map: &mut FunctionMap<R>) -> ParamTensor<S>
    where
        R: From<T>,
    {
        self.initial
            .contract_algo(|tn| tn.edge_to_min_degree_node_with_depth(depth))
            .unwrap();

        self.initial.namesym("L0");
        if self.initial.graph.n_nodes() > 1 {
            let mut new_level = self.initial.shadow();
            new_level
                .contract_algo(|tn| tn.edge_to_min_degree_node_with_depth(depth))
                .unwrap();
            self.levels.push(new_level);

            self.contract_levels(depth);
            // println!("levels {}", self.levels.len());
            self.generate_fn_map(fn_map);
            self.levels.last().unwrap().result().unwrap().0
        } else {
            self.initial
                .result_tensor_ref()
                .unwrap()
                .expanded_shadow_with_map(fn_map)
                .unwrap()
        }
    }

    fn generate_fn_map<R>(&self, fn_map: &mut FunctionMap<R>)
    where
        R: From<T>,
    {
        self.initial.append_map(fn_map);
        for l in &self.levels {
            l.append_map(fn_map);
        }
    }
}

#[cfg(feature = "shadowing")]
impl<S: TensorStructure + Clone + TracksCount + Display, Sc: AtomCore + Clone>
    Levels<ParamTensor<S>, S, Sc>
where
    ParamTensor<S>: Contract<LCM = ParamTensor<S>>,

    S: HasName<Name = Symbol, Args: IntoArgs> + ToSymbolic + StructureContract,
{
    fn contract_levels(
        &mut self,
        depth: usize,
        // fn_map: &mut FunctionMap<'a, Complex<T>>,
    ) {
        let mut not_done = true;
        let level = self.levels.len();

        if let Some(current_level) = self.levels.last_mut() {
            current_level.namesym(&format!("L{level}"))
        } else {
            not_done = false;
        }

        let nl = if let Some(current_level) = self.levels.last() {
            let mut new_level = current_level.shadow();
            new_level
                .contract_algo(|tn| tn.edge_to_min_degree_node_with_depth(depth))
                .unwrap();
            if new_level.graph.n_nodes() == 1 {
                not_done = false;
            }
            Some(new_level)
        } else {
            None
        };

        if let Some(nl) = nl {
            self.levels.push(nl);
        }

        if not_done {
            self.contract_levels(depth)
        }
    }

    pub fn contract<R>(&mut self, depth: usize, fn_map: &mut FunctionMap<R>) -> ParamTensor<S> {
        self.initial
            .contract_algo(|tn| tn.edge_to_min_degree_node_with_depth(depth))
            .unwrap();

        self.initial.namesym("L0");
        if self.initial.graph.n_nodes() > 1 {
            let mut new_level = self.initial.shadow();
            new_level
                .contract_algo(|tn| tn.edge_to_min_degree_node_with_depth(depth))
                .unwrap();
            self.levels.push(new_level);

            self.contract_levels(depth);
            // println!("levels {}", self.levels.len());
            self.generate_fn_map(fn_map);
            self.levels.last().unwrap().result().unwrap().0
        } else {
            self.initial
                .result_tensor_ref()
                .unwrap()
                .expanded_shadow_with_map(fn_map)
                .unwrap()
        }
    }

    fn generate_fn_map<R>(&self, fn_map: &mut FunctionMap<R>) {
        self.initial.append_map(fn_map);
        for l in &self.levels {
            l.append_map(fn_map);
        }
    }
}
