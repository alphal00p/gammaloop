//! Finite diagnostics through existing public diagram and tensor owners.
use std::cell::RefCell;
use std::env;
use std::error::Error;
use std::fs;
use std::io::Error as IoError;
use std::iter;
use std::ptr;
use std::result::Result as StdResult;
use std::sync::Arc;

use feynkit_graph::expressions::GraphExpressions;
use feynkit_graph::{FeynmanDiagram, NumeratorInterface, VertexId};
use feynkit_model::Model;
use idenso::tensor::SymbolicTensor;
use idenso::tensor::composition;
use idenso::tensor::composition::ProductPlan;
use idenso::{CookMode, CookSettings};
use linnet::half_edge::subgraph::SuBitGraph;
use serde_json::{Value, json};
use spenso::network::Network;
use spenso::network::library::panicing::ErroringLibrary;
use spenso::network::library::{DummyKey, DummyLibrary};
use spenso::network::parsing::{ParseSettings, ShadowedStructure};
use spenso::network::store::NetworkStore;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{PartialIndex, PartialStructure, PartialStructureExt};
use spenso::structure::slot::IsAbstractSlot;
use spenso::structure::{Canonicalized, TensorStructure};
use spenso::tensors::parametric::ParamTensor;
use symbolica::atom::{Atom, AtomCore, Symbol};

type Result<T> = StdResult<T, Box<dyn Error>>;
type Tensor = SymbolicTensor<PartialStructure>;
type ComponentTensor = ParamTensor<ShadowedStructure<AbstractIndex>>;
type ComponentNet = Network<NetworkStore<ComponentTensor, Atom>, DummyKey, Symbol, AbstractIndex>;

fn check(ok: bool, stage: &str) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(IoError::other(stage.to_owned()).into())
    }
}
fn cooking() -> CookSettings {
    CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true)
}
fn parse(source: &str) -> Result<Atom> {
    let mut values = NumeratorInterface::parse(&serde_json::to_string(&[source])?)?.0;
    Ok(values.remove(0))
}
fn admit(source: &str, ports: Option<&Value>) -> Result<Tensor> {
    let interface = ports
        .map(|v| NumeratorInterface::parse(&v.to_string()))
        .transpose()?;
    Ok(NumeratorInterface::admit(parse(source)?, interface)?)
}
fn slots(t: &Tensor) -> Vec<String> {
    t.structure()
        .logical_slots()
        .into_iter()
        .map(|s| {
            cooking()
                .uncook(composition::port_atom(s).as_view())
                .to_canonical_string()
        })
        .collect()
}
fn record(t: &Tensor) -> Value {
    json!({"ports":slots(t), "is_metric":t.is_metric(), "is_composite":t.is_composite(),
        "canonical_to_logical":(0..t.rank()).map(|i| t.structure().layout().canonical_to_logical_position(i)).collect::<Vec<_>>(),
        "canonical_ports":t.structure().canonical().external_structure().into_iter().map(|s|cooking().uncook(composition::port_atom(s).as_view()).to_canonical_string()).collect::<Vec<_>>()})
}
fn owner<'a>(g: &'a FeynmanDiagram, name: &str) -> Result<&'a Tensor> {
    if let Some(n) = name.strip_prefix("vertex:") {
        return Ok(&g
            .vertex(Some(VertexId(n.parse()?)))
            .ok_or_else(|| IoError::other("missing vertex"))?
            .numerator);
    }
    if let Some(n) = name.strip_prefix("edge:") {
        let i: usize = n.parse()?;
        return Ok(&g
            .edges()
            .find(|(id, _, _)| id.0 == i)
            .ok_or_else(|| IoError::other("missing edge"))?
            .2
            .numerator);
    }
    match name {
        "numerator" => Ok(g.numerator()?),
        "numerator_prefactor" => Ok(g.numerator_prefactor()),
        "projector" => Ok(g.projector()),
        _ => Err(IoError::other("bad owner").into()),
    }
}
fn named_plan(plan: &ProductPlan, left: &Tensor, right: &Tensor) -> Value {
    let l = slots(left);
    let r = slots(right);
    match plan {
        ProductPlan::Scalar => json!({"kind":"Scalar"}),
        ProductPlan::Outer => json!({"kind":"Outer"}),
        ProductPlan::Contract(pairs) => {
            let mut names = pairs
                .iter()
                .map(|p| json!([l[p.left], r[p.right]]))
                .collect::<Vec<_>>();
            names.sort_by_key(Value::to_string);
            json!({"kind":"Contract","names":names})
        }
        ProductPlan::Compose(a, b) => {
            json!({"kind":"Compose","left":[a.input,a.output],"right":[b.input,b.output]})
        }
    }
}
fn lowered(t: &Tensor) -> Result<Value> {
    // Exact owning lowering only; no graph contraction, Dirac trace or component substitution.
    let library = DummyLibrary::<ComponentTensor>::new();
    let functions = ErroringLibrary::<Symbol>::new();
    let network: Canonicalized<ComponentNet> =
        t.clone()
            .into_network(&library, &functions, &ParseSettings::default(), None)?;
    let mut dangling = network
        .canonical()
        .graph
        .dangling_indices()
        .into_iter()
        .map(|s| {
            cooking()
                .uncook(s.to_atom().as_view())
                .to_canonical_string()
        })
        .collect::<Vec<_>>();
    dangling.sort();
    let mut expected = slots(t);
    expected.sort();
    check(dangling == expected, "component-lowering boundary")?;
    Ok(
        json!({"dangling":dangling,"layout":format!("{:?}",network.layout()),"execution":"lowering only, no component execution"}),
    )
}
fn probe_owner(g: &FeynmanDiagram, h: &FeynmanDiagram, r: &Value) -> Result<Value> {
    let name = r["owner"]
        .as_str()
        .ok_or_else(|| IoError::other("owner schema"))?;
    let c = owner(g, name)?;
    let old = owner(h, name)?;
    check(
        slots(c) == serde_json::from_value::<Vec<String>>(r["current"].clone())?,
        "current stored ports",
    )?;
    check(
        slots(old) == serde_json::from_value::<Vec<String>>(r["historical"].clone())?,
        "historical inferred ports",
    )?;
    check(
        c.expression() == old.expression(),
        "owner mapped expression",
    )?;
    check(
        c.is_metric() == r["is_metric"].as_bool().unwrap() && old.is_metric() == c.is_metric(),
        "metric owner binding",
    )?;
    check(
        c.is_composite() == r["is_composite"].as_bool().unwrap()
            && old.is_composite() == c.is_composite(),
        "composite owner binding",
    )?;
    check(
        !c.structure().has_unresolved_ports()
            && c.structure()
                .logical_slots()
                .iter()
                .all(|s| matches!(s.aind, PartialIndex::Explicit(_))),
        "named ports",
    )?;
    let pi = serde_json::from_value::<Vec<usize>>(r["permutation"].clone())?;
    let transported = c.permuted(&pi)?;
    check(
        transported.same_value(old),
        "justified actual-owner permutation",
    )?;
    let (expr, interface) = NumeratorInterface::export(c);
    let restored = NumeratorInterface::admit(expr, interface)?;
    check(c.same_value(&restored), "owner export/admit same_value")?;
    let declared = admit(r["expression"].as_str().unwrap(), Some(&r["current"]))?;
    check(c.same_value(&declared), "explicit declared admission")?;
    let old_declared = admit(r["expression"].as_str().unwrap(), Some(&r["historical"]))?;
    check(
        old.same_value(&old_declared),
        "historical observer declaration binding",
    )?;
    let affected = pi.iter().copied().ne(0..pi.len());
    let mut result = json!({"owner":name,"affected":affected,"current":record(c),"historical":record(old),"same_value_without_transport":c.same_value(old),"transport":pi,"transport_same_value":true,"declared_admission":true,"owner_roundtrip":true});
    if affected {
        let mut wrong = (0..pi.len()).collect::<Vec<_>>();
        wrong.rotate_left(1);
        let wrong_reject = !c.permuted(&wrong)?.same_value(old);
        let mut inverse = vec![0; pi.len()];
        for (i, &j) in pi.iter().enumerate() {
            inverse[j] = i;
        }
        let mut duplicate = r["current"].as_array().unwrap().clone();
        duplicate[0] = duplicate[1].clone();
        let mut omitted = r["current"].as_array().unwrap().clone();
        omitted.remove(0);
        let mut wrong_owner = r["current"].as_array().unwrap().clone();
        wrong_owner[0] = Value::String(
            wrong_owner[0]
                .as_str()
                .unwrap()
                .replace("::hedge(", "::wrong_interface_owner("),
        );
        let mut wrong_type = r["current"].as_array().unwrap().clone();
        let first = wrong_type[0].as_str().unwrap();
        // All affected arrays begin in the self-dual bis representation. Switch its dimension.
        wrong_type[0] = Value::String(first.replacen("::bis(4,", "::bis(5,", 1));
        let source = r["expression"].as_str().unwrap();
        let reject = |v: &Vec<Value>| -> Value {
            match admit(source, Some(&json!(v))) {
                Ok(_) => json!({"accepted":true}),
                Err(e) => json!({"accepted":false,"error":e.to_string()}),
            }
        };
        result["controls"] = json!({"wrong_permutation_reject":wrong_reject,"inverse_accepted":c.permuted(&inverse)?.same_value(old),"inverse_is_same_permutation":inverse==pi,"duplicate":reject(&duplicate),"omitted":reject(&omitted),"wrong_owner":reject(&wrong_owner),"wrong_type":reject(&wrong_type),"altered_metric_claim_rejected":c.is_metric()!=!r["is_metric"].as_bool().unwrap(),"altered_composite_claim_rejected":c.is_composite()!=!r["is_composite"].as_bool().unwrap()});
        result["current_lowering"] = lowered(c)?;
        result["historical_lowering"] = lowered(old)?;
    }
    Ok(result)
}
fn public_json(data: &Value, key: &str) -> Result<String> {
    // census.rs deliberately removes id and name. Rehydrate from the separately
    // typed, exact hex ID, never through serde_json::Value's f64 number carrier.
    let id = u128::from_str_radix(data["id"].as_str().unwrap(), 16)?;
    let payload = data[key].to_string();
    Ok(format!(
        "{{\"id\":{id},\"name\":{},{}",
        data["name"],
        &payload[1..]
    ))
}
fn probe_graph(cell: &Value, data: &Value) -> Result<Value> {
    let model = Arc::new(Model::from_json(&cell["model"].to_string())?);
    let g = FeynmanDiagram::from_json(Arc::clone(&model), &public_json(data, "current_payload")?)?;
    let h = FeynmanDiagram::from_json(
        Arc::clone(&model),
        &public_json(data, "mapped_historical_payload")?,
    )?;
    check(g.id()? == h.id()?, "graph id equality")?;
    let round = FeynmanDiagram::from_json(Arc::clone(&model), &g.to_json()?)?;
    let mut owners = Vec::new();
    for r in data["owners"].as_array().unwrap() {
        let name = r["owner"].as_str().unwrap();
        check(
            owner(&g, name)?.same_value(owner(&round, name)?),
            "whole diagram roundtrip owner",
        )?;
        owners.push(
            probe_owner(&g, &h, r).map_err(|e| IoError::other(format!("owner {name}: {e}")))?,
        );
    }
    // Enumerate actual owner-pair plans, translating every positional edge to typed field identities.
    let rs = data["owners"].as_array().unwrap();
    let mut plans = Vec::new();
    for a in rs.iter().filter(|r| r["explicit"] == true) {
        let name = a["owner"].as_str().unwrap();
        let l = owner(&g, name)?;
        let hl = owner(&h, name)?;
        for b in rs.iter().filter(|r| {
            r["owner"].as_str().unwrap().starts_with("edge:")
                || r["owner"].as_str().unwrap().starts_with("vertex:")
        }) {
            let partner = b["owner"].as_str().unwrap();
            if name == partner {
                continue;
            }
            let r = owner(&g, partner)?;
            let hr = owner(&h, partner)?;
            let cp = l.product_plan(r)?;
            let hp = hl.product_plan(hr)?;
            let cn = named_plan(&cp, l, r);
            let hn = named_plan(&hp, hl, hr);
            check(cn == hn, "named plan covariance")?;
            let stale = named_plan(&hp, l, r);
            let wiring_changed = stale != hn;
            plans.push(json!({"left":name,"right":partner,"current_plan":format!("{cp:?}"),"historical_plan":format!("{hp:?}"),"physical_plan":cn,"stale_positional_plan_changes_wiring":wiring_changed}));
        }
    }
    let derived_g = g
        .clone()
        .with_numerator_fragments(iter::empty(), iter::empty())?;
    let derived_h = h
        .clone()
        .with_numerator_fragments(iter::empty(), iter::empty())?;
    check(
        derived_g.numerator()?.same_value(derived_h.numerator()?),
        "actual fragment_numerator multiply_many",
    )?;
    let graph = g.underlying();
    let full: SuBitGraph = graph.full_filter();
    let empty: SuBitGraph = graph.empty_subgraph();
    let stages = RefCell::new(Vec::<Value>::new());
    let labels = RefCell::new(Vec::<String>::new());
    let (ht,ct)=graph.try_fold_numerator(&full,&empty,(SymbolicTensor::infer(Atom::one())?,SymbolicTensor::infer(Atom::one())?),
      |id,v| {labels.borrow_mut().push(format!("vertex:{}",id.0));(h.vertex(Some(VertexId(id.0))).unwrap().numerator.clone(),v.numerator.clone())},
      |e| {let id=g.edges().find(|(_,_,candidate)|ptr::eq(*candidate,e)).unwrap().0.0;labels.borrow_mut().push(format!("edge:{id}"));(h.edges().find(|(i,_,_)|i.0==id).unwrap().2.numerator.clone(),e.numerator.clone())},
      |e|e.is_dummy,
      |(ha,ca),(hb,cb)|->Result<(Tensor,Tensor)>{
          let hp=ha.product_plan(&hb)?;let cp=ca.product_plan(&cb)?;
          let hn=named_plan(&hp,&ha,&hb);let cn=named_plan(&cp,&ca,&cb);
          check(hn==cn,"actual fold physical plan")?;
          let hv=ha.multiply(&hb)?;let cv=ca.multiply(&cb)?;
          check(hv.same_value(&cv),"actual fold output")?;
          stages.borrow_mut().push(json!({"owner":labels.borrow().last(),"historical_plan":format!("{hp:?}"),"current_plan":format!("{cp:?}"),"physical_plan":cn,"output":record(&cv),"same_value":true}));
          Ok((hv,cv))
      })?;
    check(ht.same_value(&ct), "fold result")?;
    let public = g.numerator_of(&full, &empty)?;
    let historical_public = h.numerator_of(&full, &empty)?;
    check(
        public.same_value(&ct) && public.same_value(&historical_public),
        "public graph numerator_of",
    )?;
    Ok(
        json!({"id":data["id"],"name":data["name"],"whole_graph_import":true,"fragment_numerator_multiply_many_same_value":true,"whole_graph_roundtrip":true,"owners":owners,"plans":plans,"actual_fold":stages.into_inner(),"public_numerator_of":record(&public),"public_numerator_of_same_value":true}),
    )
}
fn main() -> Result<()> {
    let args = env::args().collect::<Vec<_>>();
    check(args.len() == 3, "argv")?;
    let inputs: Value = serde_json::from_slice(&fs::read(&args[1])?)?;
    let mut result = json!({"schema":"raw-three-native-interface/1","cells":[]});
    for cell in inputs["cells"].as_array().unwrap() {
        let mut out = json!({"case":cell["case"],"mode":cell["mode"],"graphs":[]});
        for g in cell["graphs"].as_array().unwrap() {
            match probe_graph(cell, g) {
                Ok(v) => out["graphs"].as_array_mut().unwrap().push(v),
                Err(e) => {
                    out["error"] = json!({"graph":g["id"],"error":e.to_string(),"classification":"probe failure; inspect stage, not automatically semantic rejection"});
                    result["cells"].as_array_mut().unwrap().push(out);
                    fs::write(&args[2], serde_json::to_vec_pretty(&result)?)?;
                    return Err(e);
                }
            }
        }
        result["cells"].as_array_mut().unwrap().push(out);
        fs::write(&args[2], serde_json::to_vec_pretty(&result)?)?;
    }
    Ok(())
}
