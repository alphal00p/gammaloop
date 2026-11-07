use std::{
    env, fs,
    path::{Path, PathBuf},
};

use clap::Subcommand;
#[cfg(test)]
use feynkit_graph::FeynmanDiagram;
use gammalooprs::graph::Graph;
use model::ImportModel;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::sync::Arc;
use tracing::info;

use crate::{
    commands::generate::parse_process_spec_string,
    completion::CompletionArgExt,
    state::{GraphImportOptions, ProcessRef, State},
    CLISettings,
};
use color_eyre::Result;
use eyre::{eyre, Context};

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Import {
    /// Load a UFO model and make it the active model for subsequent generation.
    Model(ImportModel),
    /// Load serialized graph data into a new or existing process.
    Graphs {
        // #[arg(short = 'p')]
        /// Graph file to load, resolved against the active state root and current directory.
        #[arg(value_name = "PATH", value_hint = clap::ValueHint::FilePath)]
        source: Option<String>,

        /// Inline DOT graph content to import instead of reading from a file path.
        #[arg(long = "inline-dot", value_name = "DOT")]
        inline_dot: Option<String>,

        /// Process definition used to select physical Cutkosky cuts for imported graphs.
        #[arg(long = "process-spec", value_name = "SPEC")]
        process_spec: Option<String>,

        /// Process reference: `#<id>`, `name:<name>`, or `<id>/<name>`
        #[arg(
            short = 'p',
            long = "process",
            value_name = "PROCESS",
            completion_process_selector(crate::completion::SelectorKind::Any)
        )]
        process: Option<ProcessRef>,

        /// Name assigned to the imported integrand.
        #[arg(short = 'i', completion_disable_special_value())]
        integrand_name: Option<String>,

        /// Replace an existing process or integrand with the same target name.
        #[arg(short = 'o', default_value_t = false, conflicts_with = "append")]
        overwrite: bool,

        /// Append imported graphs to the selected target instead of replacing it.
        #[arg(short = 'a', default_value_t = false, conflicts_with = "overwrite")]
        append: bool,
    },
}

pub mod model;

impl Import {
    pub fn run(self, state: &mut State, cli_settings: &CLISettings) -> Result<()> {
        match self {
            Import::Graphs {
                source,
                inline_dot,
                process_spec,
                process,
                integrand_name,
                overwrite,
                append,
            } => {
                let source = GraphImportSource::resolve(
                    source.as_deref(),
                    inline_dot.as_deref(),
                    &cli_settings.state.folder,
                )?;
                let default_process_name = source.default_process_name()?;
                let (process_name, process_id) = match process {
                    Some(ProcessRef::Id(id)) => (None, Some(id)),
                    Some(ProcessRef::Name(name)) => (Some(name), None),
                    Some(ProcessRef::Unqualified(value)) => {
                        let name_match = state
                            .process_list
                            .processes
                            .iter()
                            .position(|p| p.definition.folder_name == value);
                        if let Ok(id) = value.parse::<usize>() {
                            if name_match.is_some() {
                                return Err(eyre!(
                                    "Ambiguous process reference '{}'. Use '#{}' or 'name:{}' to disambiguate.",
                                    value,
                                    id,
                                    value
                                ));
                            }
                            (None, Some(id))
                        } else {
                            (Some(value), None)
                        }
                    }
                    None => (Some(default_process_name), None),
                };

                info!("Loading graphs from {}", source.display_name());
                // One model handle for every source file, shared with the graphs already held.
                let graphs = source.load(&state.shared_graph_model())?;
                let generation_type = State::infer_graph_list_generation_type(&graphs)?;
                let process_definition = process_spec
                    .as_deref()
                    .map(|spec| {
                        parse_process_spec_string(spec, generation_type, &state.model)
                            .map(|spec| spec.process_definition)
                            .wrap_err_with(|| format!("Failed to parse --process-spec `{spec}`"))
                    })
                    .transpose()?;
                state.import_graphs(
                    graphs,
                    GraphImportOptions {
                        process_name,
                        process_id,
                        process_definition,
                        integrand_name,
                        overwrite,
                        append,
                    },
                )
            }
            Import::Model(im) => im.run(state),
        }
    }

    fn load_graphs(path: &Path, model: &Arc<gammalooprs::model::Model>) -> Result<Vec<Graph>> {
        if path.is_dir() {
            let mut files = fs::read_dir(path)
                .with_context(|| format!("Could not read graph directory '{}'.", path.display()))?
                .filter_map(|entry| entry.ok().map(|entry| entry.path()))
                .filter(|path| path.extension().is_some_and(|extension| extension == "dot"))
                .collect::<Vec<_>>();
            files.sort();
            if files.is_empty() {
                return Err(eyre!(
                    "No .dot files found in directory: {}",
                    path.display()
                ));
            }
            return files
                .iter()
                .map(|path| Self::load_graph_file(path, model))
                .collect::<Result<Vec<_>>>()
                .map(|sets| sets.into_iter().flatten().collect());
        }
        Self::load_graph_file(path, model)
    }

    fn load_graph_file(path: &Path, model: &Arc<gammalooprs::model::Model>) -> Result<Vec<Graph>> {
        let input = fs::read_to_string(path)
            .with_context(|| format!("Could not read graph file '{}'.", path.display()))?;
        Self::load_graph_string(&input, model)
            .with_context(|| format!("Could not import graphs from '{}'.", path.display()))
    }

    fn load_graph_string(
        input: &str,
        model: &Arc<gammalooprs::model::Model>,
    ) -> Result<Vec<Graph>> {
        let graphs = Graph::import_dot_shared(input, model)?;
        if graphs.is_empty() {
            return Err(eyre!("No diagrams found."));
        }
        Ok(graphs)
    }

    fn resolve_graph_import_path(path: &Path, state_folder: &Path) -> Result<PathBuf> {
        let cwd = env::current_dir().wrap_err(
            "Failed to query the current working directory while resolving graph import path",
        )?;

        if path.is_absolute() {
            let normalized = Self::normalize_path_lexically(path);
            return normalized
                .exists()
                .then_some(normalized)
                .ok_or_else(|| eyre!("Graph file '{}' does not exist.", path.display()));
        }

        let state_root = state_folder.parent().unwrap_or(state_folder);
        let state_root_candidate;
        let absolute_state_root = Self::normalize_path_lexically(if state_root.is_absolute() {
            state_root
        } else {
            state_root_candidate = cwd.join(state_root);
            &state_root_candidate
        });
        let state_candidate = Self::normalize_path_lexically(&absolute_state_root.join(path));
        if state_candidate.exists() {
            return Ok(state_candidate);
        }

        let cwd_candidate = Self::normalize_path_lexically(&cwd.join(path));
        if cwd_candidate.exists() {
            return Ok(cwd_candidate);
        }

        Err(eyre!(
            "Could not find graph file '{}'. Tried '{}' (active state root) and '{}' (current working directory).",
            path.display(),
            state_candidate.display(),
            cwd_candidate.display()
        ))
    }

    fn display_graph_import_path(path: &Path) -> PathBuf {
        path.canonicalize()
            .unwrap_or_else(|_| Self::normalize_path_lexically(path))
    }

    fn normalize_path_lexically(path: &Path) -> PathBuf {
        use std::path::Component;

        let mut normalized = PathBuf::new();
        for component in path.components() {
            match component {
                Component::CurDir => {}
                Component::ParentDir => {
                    if normalized.components().next_back().is_some_and(|last| {
                        !matches!(last, Component::RootDir | Component::Prefix(_))
                    }) {
                        normalized.pop();
                    } else if !path.is_absolute() {
                        normalized.push(component.as_os_str());
                    }
                }
                Component::Normal(part) => normalized.push(part),
                Component::RootDir | Component::Prefix(_) => {
                    normalized.push(component.as_os_str());
                }
            }
        }
        normalized
    }
}

/// Graph input shared by the CLI and Python imports. Loading accepts canonical,
/// runtime and compact DOT through [`Graph::import_dot`].
pub(crate) enum GraphImportSource {
    Path(PathBuf),
    String(String),
}

impl GraphImportSource {
    fn resolve(
        source: Option<&str>,
        inline_dot: Option<&str>,
        state_folder: &Path,
    ) -> Result<Self> {
        if let Some(inline_dot) = inline_dot {
            if let Some(source) = source {
                return Err(eyre!(
                    "Cannot combine graph path '{}' with --inline-dot. Use either a path or inline DOT content.",
                    source
                ));
            }
            return Ok(Self::String(inline_dot.to_string()));
        }

        let source = source.ok_or_else(|| {
            eyre!("`import graphs` requires either a graph path or --inline-dot <DOT>.")
        })?;

        Ok(Self::Path(Import::resolve_graph_import_path(
            Path::new(source),
            state_folder,
        )?))
    }

    pub(crate) fn default_process_name(&self) -> Result<String> {
        match self {
            Self::Path(path) => path
                .file_stem()
                .ok_or_else(|| {
                    eyre!(
                        "Could not derive a process name from graph path '{}'",
                        path.display()
                    )
                })
                .map(|stem| stem.to_string_lossy().into_owned()),
            Self::String(_) => Ok("inline_graphs".to_string()),
        }
    }

    fn display_name(&self) -> String {
        match self {
            Self::Path(path) => format!("'{}'", Import::display_graph_import_path(path).display()),
            Self::String(_) => "inline DOT string".to_string(),
        }
    }

    pub(crate) fn load(self, model: &Arc<gammalooprs::model::Model>) -> Result<Vec<Graph>> {
        match self {
            Self::Path(path) => Import::load_graphs(&path, model),
            Self::String(dot_string) => Import::load_graph_string(&dot_string, model),
        }
    }
}

#[cfg(test)]
mod tests {
    use feynkit_generator::{GenerationOptions, Process};
    use gammalooprs::processes::ProcessCollection;

    use super::*;

    #[test]
    fn imports_canonical_feynkit_dot_sets_with_typed_cuts() -> Result<()> {
        let model = gammalooprs::model::Model::from_json(include_str!(
            "../../../../../assets/models/json/scalars/scalars_2p_3p.json"
        ))?;
        let generated = Process::new(["scalar_1"], ["scalar_1", "scalar_1"])
            .generate_cross_section(
                Arc::new(model.clone()),
                &GenerationOptions::default()
                    .with_loop_count(1, 1)?
                    .threads(1)
                    .max_vertices(4),
            )?;
        assert!(!generated.diagrams.is_empty());
        let dot = generated
            .diagrams
            .iter()
            .map(FeynmanDiagram::to_dot)
            .collect::<Result<Vec<_>, _>>()?
            .join("\n");
        let directory = tempfile::tempdir()?;
        let path = directory.path().join("canonical.dot");
        fs::write(&path, &dot)?;

        let shared_model = Arc::new(model);
        let from_path = GraphImportSource::Path(path).load(&shared_model)?;
        let from_string = GraphImportSource::String(dot).load(&shared_model)?;

        for imported in [from_path, from_string] {
            assert_eq!(imported.len(), generated.diagrams.len());
            assert!(imported
                .iter()
                .all(|graph| !graph.finalized_cuts.is_empty()));
        }
        Ok(())
    }

    #[test]
    fn graph_import_source_accepts_compact_amplitude_dot() -> Result<()> {
        let model = gammalooprs::model::Model::from_json(include_str!(
            "../../../../../assets/models/json/scalars/scalars_2p_3p.json"
        ))?;
        let directory = tempfile::tempdir()?;
        let path = directory.path().join("scalar_box_compact.dot");
        fs::write(&path, COMPACT_SCALAR_BOX)?;

        let source = GraphImportSource::Path(path);
        assert_eq!(source.default_process_name()?, "scalar_box_compact");
        let shared_model = Arc::new(model);
        let from_path = source.load(&shared_model)?;
        let from_string =
            GraphImportSource::String(COMPACT_SCALAR_BOX.to_owned()).load(&shared_model)?;

        for imported in [from_path, from_string] {
            assert_eq!(imported.len(), 1);
            assert_eq!(imported[0].name, "scalar_box");
        }
        Ok(())
    }

    const COMPACT_SCALAR_BOX: &str =
        include_str!("../../../../../tests/resources/graphs/scalar_box_scalar_compact.dot");
    #[test]
    fn directory_import_shares_one_model_across_files_and_reload() -> Result<()> {
        gammalooprs::initialisation::test_initialise()?;
        let temp = tempfile::tempdir()?;
        let mut state = State::new_test();
        let mut cli_settings = CLISettings::default();
        cli_settings.state.folder = temp.path().join("state");
        ImportModel {
            path: "scalars".into(),
            simplify_model: true,
        }
        .run(&mut state)?;
        let directory = temp.path().join("graphs");
        fs::create_dir(&directory)?;
        let source = fs::read_to_string(
            crate::test_workspace_root().join("tests/resources/graphs/scalar_bubble.dot"),
        )?;
        for name in ["bubble_a", "bubble_b"] {
            fs::write(
                directory.join(format!("{name}.dot")),
                source.replacen("digraph bubble", &format!("digraph {name}"), 1),
            )?;
        }
        Import::Graphs {
            source: Some(directory.display().to_string()),
            inline_dot: None,
            process_spec: None,
            process: Some(ProcessRef::Name("bubbles".into())),
            integrand_name: Some("default".into()),
            overwrite: false,
            append: false,
        }
        .run(&mut state, &cli_settings)?;

        let models = |state: &State| {
            let ProcessCollection::Amplitudes(amplitudes) =
                &state.process_list.processes[0].collection
            else {
                unreachable!("bubble graphs are amplitudes")
            };
            amplitudes["default"]
                .graphs
                .iter()
                .map(|graph| (graph.graph.name.clone(), Arc::clone(&graph.graph.model)))
                .collect::<Vec<_>>()
        };
        let imported = models(&state);
        assert_eq!(
            imported
                .iter()
                .map(|(name, _)| name.as_str())
                .collect::<Vec<_>>(),
            ["bubble_a", "bubble_b"]
        );
        assert!(
            Arc::ptr_eq(&imported[0].1, &imported[1].1),
            "graphs from two files of one directory must share one model"
        );

        state.save(&cli_settings.state.folder, true, false)?;
        let loaded = State::load(cli_settings.state.folder.clone(), None, None)?;
        let reloaded = models(&loaded);
        assert_eq!(reloaded.len(), 2);
        assert!(
            Arc::ptr_eq(&reloaded[0].1, &reloaded[1].1),
            "reloaded graphs must share one model"
        );
        assert_eq!(reloaded[0].1.fingerprint(), loaded.model.fingerprint());
        Ok(())
    }
}
