//! Canonization keys: element payloads, or what the key callbacks make of them, as ids that
//! identify equal keys and as ranks in an order that only depends on the key values.
use std::cmp::Ordering;
use std::collections::hash_map::{Entry, HashMap};

use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyFloat, PyInt, PyString, PyTuple};

use super::PyHedgeGraph;

/// Symbolica expressions are keys of their own kind: compared, ordered and hashed structurally
/// as atoms, so never equal to the numbers they evaluate to.
#[cfg(feature = "community")]
mod expression {
    use pyo3::prelude::*;
    use symbolica::api::python::PythonExpression;
    pub(super) use symbolica::atom::Atom;

    use super::{PayloadKeys, Term};

    impl PayloadKeys<'_, '_> {
        /// The id of an expression key, or `None` for any other value.
        pub(super) fn expression(&mut self, value: &Bound<'_, PyAny>) -> Option<u32> {
            let expression = value.cast::<PythonExpression>().ok()?.borrow();
            let atom = &expression.expr;
            if let Some(&id) = self.atoms.get(atom) {
                return Some(id);
            }
            let id = self.push(Term::Expression(atom.clone()));
            self.atoms.insert(atom.clone(), id);
            Some(id)
        }
    }
}

/// Key kinds in their sort order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Kind {
    None,
    Number,
    Str,
    Bytes,
    Tuple,
    #[cfg(feature = "community")]
    Expression,
    Other,
}

const KINDS: usize = Kind::Other as usize + 1;

/// A distinct key value. Python values of the first four kinds are ordered by Python's `<`,
/// other Python values by qualified type name and repr, tuples element-wise and expressions as
/// atoms.
enum Term {
    /// A Python value, interned in the dictionary of its kind.
    Object(Kind),
    Tuple(Box<[u32]>),
    #[cfg(feature = "community")]
    Expression(expression::Atom),
}

impl Term {
    fn kind(&self) -> Kind {
        match self {
            Self::Object(kind) => *kind,
            Self::Tuple(_) => Kind::Tuple,
            #[cfg(feature = "community")]
            Self::Expression(_) => Kind::Expression,
        }
    }
}

/// Node, edge and half-edge keys, indexed like the graph's elements.
pub(crate) struct ElementKeys<T> {
    pub(crate) nodes: Vec<T>,
    pub(crate) edges: Vec<T>,
    pub(crate) hedges: Vec<T>,
}

/// Projects payloads, or the results of the key callbacks, to key ids shared by all graphs read
/// through it. Python values are interned by `__eq__` and `__hash__` within their kind, tuples
/// element-wise and expressions structurally, so that every distinct key is ordered only once.
pub(crate) struct PayloadKeys<'a, 'py> {
    py: Python<'py>,
    node_key: Option<&'a Bound<'py, PyAny>>,
    edge_key: Option<&'a Bound<'py, PyAny>>,
    half_edge_key: Option<&'a Bound<'py, PyAny>>,
    /// Python values by kind, mapped to their ids.
    objects: [Option<Bound<'py, PyDict>>; KINDS],
    /// Every occurrence of a value of the other kind, whose type name and repr must agree for
    /// equal values and differ for distinct ones to order them canonically.
    others: Vec<(&'static str, u32, Bound<'py, PyAny>)>,
    #[cfg(feature = "community")]
    atoms: HashMap<expression::Atom, u32>,
    tuples: HashMap<Box<[u32]>, u32>,
    terms: Vec<Term>,
    real: Option<Bound<'py, PyAny>>,
}

impl<'a, 'py> PayloadKeys<'a, 'py> {
    pub(crate) fn new(
        py: Python<'py>,
        node_key: Option<&'a Bound<'py, PyAny>>,
        edge_key: Option<&'a Bound<'py, PyAny>>,
        half_edge_key: Option<&'a Bound<'py, PyAny>>,
    ) -> Self {
        Self {
            py,
            node_key,
            edge_key,
            half_edge_key,
            objects: Default::default(),
            others: vec![],
            #[cfg(feature = "community")]
            atoms: HashMap::new(),
            tuples: HashMap::new(),
            terms: vec![],
            real: None,
        }
    }

    /// Key ids of every element: equal ids for equal keys, in no particular order.
    pub(crate) fn read(&mut self, graph: &PyHedgeGraph) -> PyResult<ElementKeys<u32>> {
        let nodes = graph
            .iter_nodes()
            .map(|(_, _, n)| self.key("node", self.node_key, &n.data))
            .collect::<PyResult<_>>()?;
        Ok(ElementKeys {
            nodes,
            ..self.read_edges(graph)?
        })
    }

    /// Edge and half-edge key ids only; node payloads are not read.
    pub(crate) fn read_edges(&mut self, graph: &PyHedgeGraph) -> PyResult<ElementKeys<u32>> {
        Ok(ElementKeys {
            nodes: vec![],
            edges: graph
                .iter_edges()
                .map(|(_, _, e)| self.key("edge", self.edge_key, &e.data.data))
                .collect::<PyResult<_>>()?,
            hedges: graph
                .iter_hedges()
                .map(|(_, h)| self.key("half_edge", self.half_edge_key, &h.data))
                .collect::<PyResult<_>>()?,
        })
    }

    fn key(
        &mut self,
        element: &'static str,
        callback: Option<&Bound<'py, PyAny>>,
        data: &Py<PyAny>,
    ) -> PyResult<u32> {
        let data = data.bind(self.py);
        match callback {
            Some(f) => self.term(element, &f.call1((data,))?),
            None => self.term(element, data),
        }
    }

    fn push(&mut self, term: Term) -> u32 {
        self.terms.push(term);
        self.terms.len() as u32 - 1
    }

    fn term(&mut self, element: &'static str, value: &Bound<'py, PyAny>) -> PyResult<u32> {
        #[cfg(feature = "community")]
        if let Some(id) = self.expression(value) {
            return Ok(id);
        }
        if let Ok(items) = value.cast::<PyTuple>() {
            let items = items
                .iter()
                .map(|item| self.term(element, &item))
                .collect::<PyResult<Box<[u32]>>>()?;
            if let Some(&id) = self.tuples.get(&items) {
                return Ok(id);
            }
            let id = self.push(Term::Tuple(items.clone()));
            self.tuples.insert(items, id);
            return Ok(id);
        }
        let kind = self.kind(value)?;
        let py = self.py;
        let objects = self.objects[kind as usize]
            .get_or_insert_with(|| PyDict::new(py))
            .clone();
        let found = objects
            .get_item(value)
            .map_err(|e| unhashable(element, value, e))?;
        let id = match found {
            Some(id) => id.extract()?,
            None => {
                if kind == Kind::Number && value.ne(value)? {
                    return Err(PyValueError::new_err(format!(
                        "NaN is not a valid {element} key; {}",
                        hint(element, "comparable")
                    )));
                }
                let id = self.push(Term::Object(kind));
                objects.set_item(value, id)?;
                id
            }
        };
        if kind == Kind::Other {
            self.others.push((element, id, value.clone()));
        }
        Ok(id)
    }

    /// The kind of a Python value that is neither a tuple nor an expression.
    fn kind(&mut self, value: &Bound<'py, PyAny>) -> PyResult<Kind> {
        Ok(if value.is_none() {
            Kind::None
        } else if value.is_instance_of::<PyInt>() || value.is_instance_of::<PyFloat>() {
            Kind::Number
        } else if value.is_instance_of::<PyString>() {
            Kind::Str
        } else if value.is_instance_of::<PyBytes>() {
            Kind::Bytes
        } else {
            let real = match &self.real {
                Some(real) => real,
                None => self
                    .real
                    .insert(self.py.import("numbers")?.getattr("Real")?),
            };
            if value.is_instance(real)? {
                Kind::Number
            } else {
                Kind::Other
            }
        })
    }

    /// Ids of the other values in the order of their qualified type names and reprs. Fails
    /// unless these identify the values' equality classes, so that the order is canonical.
    fn order_others(&self) -> PyResult<Vec<u32>> {
        let default_repr = self
            .py
            .import("builtins")?
            .getattr("object")?
            .getattr("__repr__")?;
        let mut described = HashMap::new();
        for &(element, id, ref value) in &self.others {
            let ty = value.get_type();
            if ty.getattr("__repr__")?.is(&default_repr) {
                return Err(PyTypeError::new_err(format!(
                    "{element} key of type {ty} has the default repr, which does not order it \
                     canonically; {}",
                    hint(element, "ordered")
                )));
            }
            let description = (
                ty.fully_qualified_name()?.to_string(),
                value.repr()?.to_string(),
            );
            match described.entry(id) {
                Entry::Vacant(entry) => {
                    entry.insert((description, element));
                }
                Entry::Occupied(entry) if entry.get().0 != description => {
                    return Err(PyTypeError::new_err(format!(
                        "equal {element} keys {} and {} differ in type or repr, which does not \
                         order them canonically; {}",
                        entry.get().0 .1,
                        description.1,
                        hint(element, "ordered")
                    )));
                }
                Entry::Occupied(_) => {}
            }
        }
        let mut order = described
            .into_iter()
            .map(|(id, (description, element))| (description, element, id))
            .collect::<Vec<_>>();
        order.sort_unstable();
        if let Some([_, (description, element, _)]) =
            order.windows(2).find(|pair| pair[0].0 == pair[1].0)
        {
            return Err(PyTypeError::new_err(format!(
                "distinct {element} keys share the type and repr {}, which does not order them \
                 canonically; {}",
                description.1,
                hint(element, "ordered")
            )));
        }
        Ok(order.into_iter().map(|(.., id)| id).collect())
    }

    /// A total order on distinct key ids, given the position of each Python value in its kind.
    fn cmp(&self, position: &[u32], a: u32, b: u32) -> Ordering {
        match (&self.terms[a as usize], &self.terms[b as usize]) {
            (Term::Tuple(x), Term::Tuple(y)) => x
                .iter()
                .zip(y.iter())
                .map(|(&a, &b)| self.cmp(position, a, b))
                .find(|order| order.is_ne())
                .unwrap_or_else(|| x.len().cmp(&y.len())),
            #[cfg(feature = "community")]
            (Term::Expression(x), Term::Expression(y)) => x.cmp(y),
            (x, y) => x
                .kind()
                .cmp(&y.kind())
                .then(position[a as usize].cmp(&position[b as usize])),
        }
    }

    /// Ranks of the key ids in their value order, so that they never depend on the element
    /// order. Only canonical forms need this; equal ids for equal keys suffice to compare graphs.
    pub(crate) fn rank(&self, keys: ElementKeys<u32>) -> PyResult<ElementKeys<u32>> {
        let mut position = vec![0; self.terms.len()];
        for objects in self.objects[..Kind::Other as usize].iter().flatten() {
            // Distinct values never tie on `==`, so pairs sort by value; Python's sort, unlike
            // Rust's, tolerates an inconsistent `<`.
            let items = objects.items();
            items.sort()?;
            for (p, item) in items.iter().enumerate() {
                let (_, id) = item.extract::<(Bound<'_, PyAny>, u32)>()?;
                position[id as usize] = p as u32;
            }
        }
        for (p, id) in self.order_others()?.into_iter().enumerate() {
            position[id as usize] = p as u32;
        }
        let mut order = (0..self.terms.len() as u32).collect::<Vec<_>>();
        order.sort_unstable_by(|&a, &b| self.cmp(&position, a, b));
        let mut rank = vec![0; order.len()];
        for (r, term) in order.into_iter().enumerate() {
            rank[term as usize] = r as u32;
        }
        let ranked = |terms: Vec<u32>| terms.into_iter().map(|t| rank[t as usize]).collect();
        Ok(ElementKeys {
            nodes: ranked(keys.nodes),
            edges: ranked(keys.edges),
            hedges: ranked(keys.hedges),
        })
    }
}

fn hint(element: &str, keys: &str) -> String {
    format!("pass {element}_key= to map payloads to {keys} keys")
}

/// `error` from interning `value`, explained if the value's type is unhashable.
fn unhashable(element: &str, value: &Bound<'_, PyAny>, error: PyErr) -> PyErr {
    let (py, ty) = (value.py(), value.get_type());
    let unhashable = ty.getattr("__hash__").is_ok_and(|hash| hash.is_none());
    if !unhashable || !error.is_instance_of::<PyTypeError>(py) {
        return error;
    }
    let explained = PyTypeError::new_err(format!(
        "unhashable {element} key of type {ty}; {}",
        hint(element, "hashable")
    ));
    explained.set_cause(py, Some(error));
    explained
}
