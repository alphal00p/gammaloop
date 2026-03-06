use std::{
    fs::{self, File},
    path::Path,
};

use color_eyre::Result;
use eyre::{Context, eyre};
use linnet::half_edge::involution::Orientation;
use rand::Rng;
use symbolica::{
    atom::{Atom, AtomCore},
    printer::PrintOptions,
    state::State,
};

use crate::{
    integrands::process::{
        ProcessIntegrandImpl, GenericEvaluator, GenericEvaluatorFloat,
        amplitude::{
            AmplitudeIntegrand,
            load::{
                STANDALONE_EVALUATORS_VERSION, StandaloneEvaluatorArchive,
                StandaloneEvaluatorStackArchive, StandaloneGenericEvaluatorArchive,
                StandaloneGraphTermArchive,
            },
        },
    },
    momentum::ThreeMomentum,
    momentum::sample::{LoopMomenta, MomentumSample},
    processes::{StandaloneDataFormat, StandaloneExportMode, StandaloneExportSettings},
    utils::F,
};

const STANDALONE_DATA_FILE: &str = "standalone_evaluators";
const STANDALONE_RUST_SCRIPT_FILE: &str = "standalone_evaluators_rust.rs";

fn atom_to_bytes(atom: &Atom) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    atom.as_view()
        .write(&mut out)
        .with_context(|| "Failed to write atom to standalone binary stream")?;
    Ok(out)
}

fn export_generic_evaluator<T: ExportAtomTo>(
    evaluator: &GenericEvaluator,
) -> Result<StandaloneGenericEvaluatorArchive<T>> {
    let exprs = evaluator
        .exprs
        .as_ref()
        .ok_or_else(|| {
            eyre!("Standalone export requires stored evaluator atoms (`store_atom=true`)")
        })?
        .iter()
        .map(|a| T::export_atom_to(a))
        .collect::<Result<Vec<_>>>()?;

    let additional_fn_map_entries = evaluator
        .fn_map_entries
        .iter()
        .map(|e| e.archive())
        .collect::<Result<_>>()?;

    Ok(StandaloneGenericEvaluatorArchive {
        exprs,
        additional_fn_map_entries,
        dual_shape: evaluator.dual_shape.clone(),
    })
}

fn save_atom_export(atom: &Atom, path: &Path) -> Result<()> {
    let mut file = File::create(path)
        .with_context(|| format!("Failed to create atom export file {}", path.display()))?;
    atom.as_view()
        .export(&mut file)
        .with_context(|| format!("Failed to export atom to {}", path.display()))?;
    Ok(())
}

fn standalone_python_script() -> &'static str {
    r#"#!/usr/bin/env -S uv run --script
# /// script
# dependencies = ["symbolica"]
# ///
import json
from pathlib import Path
from symbolica import Expression, Replacement


def _load_expr(data_root: Path, rel: str):
    return Expression.load(str(data_root / rel))


def _make_evaluator(data_root: Path, evaluator_manifest: dict, params, replacements):
    exprs = [
        _load_expr(data_root, rel).replace_multiple(replacements)
        for rel in evaluator_manifest["exprs"]
    ]
    constants = {}
    functions = {}
    return Expression.evaluator_multiple(
        exprs,
        constants,
        functions,
        params,
        iterations=evaluator_manifest.get("horner_iterations", 100),
        n_cores=evaluator_manifest.get("n_cores", 4),
        verbose=evaluator_manifest.get("verbose", False),
    )


def load_standalone(data_root: str | Path | None = None):
    data_root = Path(data_root) if data_root is not None else (Path(__file__).resolve().parent / "data")
    manifest = json.loads((data_root / "manifest.json").read_text())
    loaded = []

    for graph in manifest["graphs"]:
        params = [_load_expr(data_root, rel) for rel in graph["params"]]
        replacements = [
            Replacement(_load_expr(data_root, rep["lhs"]), _load_expr(data_root, rep["rhs"]))
            for rep in graph["replacements"]
        ]

        original = graph["original_integrand"]
        loaded_graph = {
            "graph_name": graph["graph_name"],
            "single_parametric": _make_evaluator(data_root, original["single_parametric"], params, replacements),
            "iterative": _make_evaluator(data_root, original["iterative"], params, replacements) if original["iterative"] else None,
            "summed_function_map": _make_evaluator(data_root, original["summed_function_map"], params, replacements) if original["summed_function_map"] else None,
            "summed": _make_evaluator(data_root, original["summed"], params, replacements) if original["summed"] else None,
            "threshold_counterterms": [],
        }

        for ct in graph["threshold_counterterms"]:
            loaded_graph["threshold_counterterms"].append(
                {
                    "parametric": _make_evaluator(data_root, ct["parametric"], params, replacements),
                    "iterative": _make_evaluator(data_root, ct["iterative"], params, replacements) if ct["iterative"] else None,
                }
            )

        loaded.append(loaded_graph)

    return loaded

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Load standalone Symbolica evaluators")
    parser.add_argument("data_root", nargs="?", default=None, help="Optional path to data directory")
    args = parser.parse_args()

    data_root = Path(args.data_root) if args.data_root is not None else None
    graphs = load_standalone(data_root)
    print(f"Loaded {len(graphs)} graph evaluator sets")
"#
}

fn standalone_rust_script() -> String {
    let mut script = include_str!("load.rs").to_string();
    if let Some(rest) = script.strip_prefix("//#!/usr/bin/env -S rust-script\n") {
        script = format!("#!/usr/bin/env -S rust-script\n{rest}");
    }

    script
}

pub trait CreateArchive<T, S> {
    fn create_archive(&self, state: S) -> Result<StandaloneEvaluatorArchive<S, T>>;
}

pub trait ExportAtomTo: Sized {
    fn export_atom_to(atom: &Atom) -> Result<Self>;
}

impl ExportAtomTo for Vec<u8> {
    fn export_atom_to(atom: &Atom) -> Result<Self> {
        atom_to_bytes(atom)
    }
}

impl ExportAtomTo for String {
    fn export_atom_to(atom: &Atom) -> Result<Self> {
        Ok(atom.printer(PrintOptions::file()).to_string())
    }
}

impl<T: ExportAtomTo, S> CreateArchive<T, S> for AmplitudeIntegrand {
    fn create_archive(&self, symbolica_state: S) -> Result<StandaloneEvaluatorArchive<S, T>> {
        let sample_inputs = self.representative_input()?;

        let graph_terms = self
            .data
            .graph_terms
            .iter()
            .zip(sample_inputs)
            .map(
                |(term, (representative_input, mult_offset, start, override_pos))| {
                    let param_builder_params = (&term.param_builder.pairs)
                        .into_iter()
                        .flat_map(|pair| pair.params.iter())
                        .map(|a| T::export_atom_to(a))
                        .collect::<Result<Vec<_>>>()?;

                    let orientations: Vec<Vec<_>> = term
                        .orientations
                        .iter()
                        .map(|a| {
                            a.iter()
                                .map(|(a, o)| match o {
                                    Orientation::Default => 1i8,
                                    Orientation::Reversed => -1i8,
                                    Orientation::Undirected => 0,
                                })
                                .collect()
                        })
                        .collect();

                    let fn_map_entries = term
                        .param_builder
                        .reps
                        .iter()
                        .map(|entry| entry.archive())
                        .collect::<Result<Vec<_>>>()?;

                    let original_integrand = StandaloneEvaluatorStackArchive {
                        start,
                        mult_offset,
                        override_pos,
                        representative_input: representative_input.clone(),
                        single_parametric: export_generic_evaluator(
                            &term.original_integrand.single_parametric,
                        )?,
                        iterative: term
                            .original_integrand
                            .iterative
                            .as_ref()
                            .map(|(eval, len)| export_generic_evaluator(eval))
                            .transpose()?,
                        summed_function_map: term
                            .original_integrand
                            .summed_function_map
                            .as_ref()
                            .map(export_generic_evaluator)
                            .transpose()?,
                        summed: term
                            .original_integrand
                            .summed
                            .as_ref()
                            .map(export_generic_evaluator)
                            .transpose()?,
                    };

                    let threshold_counterterms = term
                        .threshold_counterterm
                        .evaluators
                        .iter()
                        .map(|counterterm| {
                            let evaluator_stack = counterterm.evaluator_stack.borrow();
                            Ok(StandaloneEvaluatorStackArchive {
                                start,
                                override_pos,
                                mult_offset,
                                representative_input: representative_input.clone(),
                                single_parametric: export_generic_evaluator(
                                    &evaluator_stack.single_parametric,
                                )?,
                                iterative: evaluator_stack
                                    .iterative
                                    .as_ref()
                                    .map(|(eval, len)| export_generic_evaluator(eval))
                                    .transpose()?,
                                summed_function_map: evaluator_stack
                                    .summed_function_map
                                    .as_ref()
                                    .map(export_generic_evaluator)
                                    .transpose()?,
                                summed: evaluator_stack
                                    .summed
                                    .as_ref()
                                    .map(export_generic_evaluator)
                                    .transpose()?,
                            })
                        })
                        .collect::<Result<Vec<_>>>()?;

                    Ok(StandaloneGraphTermArchive {
                        graph_name: term.graph.name.clone(),
                        orientations,
                        param_builder_params,
                        fn_map_entries,
                        original_integrand,
                        threshold_counterterms,
                    })
                },
            )
            .collect::<Result<Vec<_>>>()?;
        Ok(StandaloneEvaluatorArchive {
            version: STANDALONE_EVALUATORS_VERSION,
            symbolica_state,
            graph_terms,
        })
    }
}

impl AmplitudeIntegrand {
    pub(crate) fn representative_input(
        &self,
    ) -> Result<
        Vec<(
            Vec<symbolica::domains::float::Complex<f64>>,
            usize,
            usize,
            usize,
        )>,
    > {
        let hel = self.settings.kinematics.externals.get_helicities();
        let mut rng = rand::rng();

        let mut mom_samples = Vec::new();
        {
            let dependent_momenta_constructor = self.get_dependent_momenta_constructor();
            for t in self.data.graph_terms.iter() {
                let loop_moms = LoopMomenta::from_iter(
                    (0..t.graph.loop_momentum_basis.loop_edges.len()).map(|_| ThreeMomentum {
                        px: F(rng.random_range(-1.0..1.0)),
                        py: F(rng.random_range(-1.0..1.0)),
                        pz: F(rng.random_range(-1.0..1.0)),
                    }),
                );
                let momentum_sample = MomentumSample::new(
                    loop_moms,
                    0,
                    &self.settings.kinematics.externals,
                    0,
                    F(1.0),
                    dependent_momenta_constructor,
                    None,
                )?;
                mom_samples.push(momentum_sample);
            }
        }
        let mut inputs = Vec::new();

        for (t, momentum_sample) in self.data.graph_terms.iter().zip(&mom_samples) {
            let mut p_build = t.param_builder.clone();
            let input = f64::get_parameters(
                &mut p_build,
                (false, false),
                &t.graph,
                momentum_sample,
                hel,
                &self.settings.additional_params(),
                None,
                None,
                None,
            );

            inputs.push((
                input
                    .as_slice()
                    .iter()
                    .map(|c| symbolica::domains::float::Complex::new(c.re.0, c.im.0))
                    .collect(),
                input.multiplicative_offset,
                input.orientations_start,
                input.override_pos,
            ));
        }

        Ok(inputs)
    }

    // pub(crate) fn create_archive<

    pub(crate) fn export_standalone(
        &self,
        path: impl AsRef<Path>,
        settings: &StandaloneExportSettings,
    ) -> Result<()> {
        let mut symbolica_state = Vec::new();
        State::export(&mut symbolica_state)
            .with_context(|| "Failed to export Symbolica state for standalone evaluators")?;

        let mut standalone_path = path.as_ref().join(STANDALONE_DATA_FILE);

        match settings.format {
            StandaloneDataFormat::Binary => {
                let standalone: StandaloneEvaluatorArchive<Vec<u8>, Vec<u8>> =
                    self.create_archive(symbolica_state)?;
                let binary = bincode::encode_to_vec(&standalone, bincode::config::standard())?;
                standalone_path.add_extension("bin");
                fs::write(&standalone_path, binary).with_context(|| {
                    format!(
                        "Failed writing standalone evaluator binary to {}",
                        standalone_path.display()
                    )
                })?;
            }
            StandaloneDataFormat::Json => {
                let standalone: StandaloneEvaluatorArchive<(), String> = self.create_archive(())?;
                let json = serde_json::to_vec_pretty(&standalone)?;
                standalone_path.add_extension("json");
                fs::write(&standalone_path, json).with_context(|| {
                    format!(
                        "Failed writing standalone evaluator JSON to {}",
                        standalone_path.display()
                    )
                })?
            }
        }

        match settings.mode {
            StandaloneExportMode::Python => {
                return Err(eyre!("Python Export mode not implemented"));
            }
            StandaloneExportMode::Rust => {
                let script = standalone_rust_script();
                let script_path = path.as_ref().join(STANDALONE_RUST_SCRIPT_FILE);
                fs::write(&script_path, script).with_context(|| {
                    format!(
                        "Failed writing standalone rust-script loader to {}",
                        script_path.display()
                    )
                })?;
            }
        }

        Ok(())
    }

    fn export_standalone_python_directory(
        &self,
        path: &Path,
        archive: &StandaloneEvaluatorArchive,
    ) -> Result<()> {
        Err(eyre!("Not yet implemented"))
        // let root = path.join(STANDALONE_PYTHON_DIR);
        // let data_root = root.join("data");
        // fs::create_dir_all(&data_root)
        //     .with_context(|| format!("Failed creating {}", data_root.display()))?;
        // fs::write(
        //     root.join(STANDALONE_PYTHON_SCRIPT_FILE),
        //     standalone_python_script(),
        // )
        // .with_context(|| "Failed writing python standalone script")?;

        // let graphs = archive
        //     .graph_terms
        //     .iter()
        //     .enumerate()
        //     .map(|(graph_index, graph)| {
        //         let graph_dir = data_root.join(format!("graph_{graph_index}"));
        //         fs::create_dir_all(&graph_dir).with_context(|| {
        //             format!(
        //                 "Failed creating graph standalone directory {}",
        //                 graph_dir.display()
        //             )
        //         })?;

        //         let params_dir = graph_dir.join("params");
        //         fs::create_dir_all(&params_dir)?;
        //         let mut params = Vec::new();
        //         for (i, bytes) in graph.param_builder_params.iter().enumerate() {
        //             let rel = format!("graph_{graph_index}/params/param_{i}.bin");
        //             fs::write(data_root.join(&rel), bytes).with_context(|| {
        //                 format!("Failed writing {}", data_root.join(&rel).display())
        //             })?;
        //             params.push(rel);
        //         }

        //         let reps_dir = graph_dir.join("replacements");
        //         fs::create_dir_all(&reps_dir)?;
        //         let mut replacements = Vec::new();
        //         for (i, (lhs, rhs)) in graph.param_builder_replacements.iter().enumerate() {
        //             let lhs_rel = format!("graph_{graph_index}/replacements/lhs_{i}.bin");
        //             let rhs_rel = format!("graph_{graph_index}/replacements/rhs_{i}.bin");
        //             fs::write(data_root.join(&lhs_rel), lhs).with_context(|| {
        //                 format!("Failed writing {}", data_root.join(&lhs_rel).display())
        //             })?;
        //             fs::write(data_root.join(&rhs_rel), rhs).with_context(|| {
        //                 format!("Failed writing {}", data_root.join(&rhs_rel).display())
        //             })?;
        //             replacements.push(PythonReplacementManifest {
        //                 lhs: lhs_rel,
        //                 rhs: rhs_rel,
        //             });
        //         }

        //         let write_eval =
        //             |name: &str,
        //              eval: &StandaloneGenericEvaluatorArchive|
        //              -> Result<(Vec<String>, PythonEvaluatorManifest)> {
        //                 let eval_dir = graph_dir.join(name);
        //                 fs::create_dir_all(&eval_dir)?;
        //                 let mut expr_files = Vec::new();
        //                 for (i, expr) in eval.exprs.iter().enumerate() {
        //                     let rel = format!("graph_{graph_index}/{name}/expr_{i}.bin");
        //                     fs::write(data_root.join(&rel), expr).with_context(|| {
        //                         format!("Failed writing {}", data_root.join(&rel).display())
        //                     })?;
        //                     expr_files.push(rel);
        //                 }
        //                 Ok((expr_files.clone(), python_eval_manifest(eval, expr_files)))
        //             };

        //         let (_, single_parametric) = write_eval(
        //             "single_parametric",
        //             &graph.original_integrand.single_parametric,
        //         )?;
        //         let iterative = graph
        //             .original_integrand
        //             .iterative
        //             .as_ref()
        //             .map(|e| write_eval("iterative", e).map(|(_, m)| m))
        //             .transpose()?;
        //         let summed_function_map = graph
        //             .original_integrand
        //             .summed_function_map
        //             .as_ref()
        //             .map(|e| write_eval("summed_function_map", e).map(|(_, m)| m))
        //             .transpose()?;
        //         let summed = graph
        //             .original_integrand
        //             .summed
        //             .as_ref()
        //             .map(|e| write_eval("summed", e).map(|(_, m)| m))
        //             .transpose()?;

        //         let mut threshold_counterterms = Vec::new();
        //         for (ct_id, ct) in graph.threshold_counterterms.iter().enumerate() {
        //             let (_, parametric) =
        //                 write_eval(&format!("ct_{ct_id}_parametric"), &ct.parametric)?;
        //             let iterative = ct
        //                 .iterative
        //                 .as_ref()
        //                 .map(|e| write_eval(&format!("ct_{ct_id}_iterative"), e).map(|(_, m)| m))
        //                 .transpose()?;
        //             threshold_counterterms.push(PythonCountertermManifest {
        //                 parametric,
        //                 iterative,
        //             });
        //         }

        //         Ok(PythonStandaloneGraphManifest {
        //             graph_name: graph.graph_name.clone(),
        //             params,
        //             replacements,
        //             original_integrand: PythonEvaluatorStackManifest {
        //                 single_parametric,
        //                 iterative,
        //                 summed_function_map,
        //                 summed,
        //             },
        //             threshold_counterterms,
        //         })
        //     })
        //     .collect::<Result<Vec<_>>>()?;

        // let manifest = PythonStandaloneManifest { graphs };
        // let manifest_path = data_root.join(STANDALONE_PYTHON_MANIFEST_FILE);
        // let manifest_json = serde_json::to_vec_pretty(&manifest)?;
        // fs::write(&manifest_path, manifest_json)
        //     .with_context(|| format!("Failed writing {}", manifest_path.display()))?;

        // Ok(())
    }
}
