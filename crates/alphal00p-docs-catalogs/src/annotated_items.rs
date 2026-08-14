//! Source-backed proc-macro descriptors for the supported Rust surfaces.
//!
//! Each marker is consumed by the attribute macro. At compile time the macro
//! parses the named product source item, captures its real signature, members,
//! and existing Rustdoc, and emits a static descriptor plus a dependency edge.
//! Product crates therefore remain independently publishable: only this
//! internal adapter depends on the documentation crates.

use alphal00p_docs_schema::DocItem;

#[alphal00p_docs::trait_item(
    id = "GammaLoopContext",
    title = "GammaLoopContext",
    summary = "Combines the state-map and model capabilities required by GammaLoop computations.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/lib.rs",
    source_id = "gammalooprs::GammaLoopContext"
)]
fn gammalooprs_gamma_loop_context() {}

#[alphal00p_docs::trait_item(
    id = "HasModel",
    title = "HasModel",
    summary = "Provides shared access to the physics model owned by a computation context.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/lib.rs",
    source_id = "gammalooprs::HasModel"
)]
fn gammalooprs_has_model() {}

#[alphal00p_docs::func(
    id = "set_interrupt_handler",
    title = "Install interrupt handling",
    summary = "Installs GammaLoop's process interrupt handler before a long-running integration.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/lib.rs",
    source_id = "gammalooprs::set_interrupt_handler"
)]
fn gammalooprs_set_interrupt_handler() {}

#[alphal00p_docs::ty(
    id = "StateLoadOption",
    title = "State loading options",
    summary = "Selects the state folder, boot commands, settings overrides, logging, and read-only mode.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/lib.rs",
    source_id = "gammaloop_api::StateLoadOption"
)]
fn gammaloop_api_state_load_option() {}

#[alphal00p_docs::func(
    id = "StateLoadOption::load",
    title = "Load a GammaLoop state",
    summary = "Loads the selected state and returns the supported session boundary.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/lib.rs",
    source_id = "gammaloop_api::StateLoadOption::load"
)]
fn gammaloop_api_state_load_option_load() {}

#[alphal00p_docs::func(
    id = "LoadedState::cli_session",
    title = "Create a CLI session",
    summary = "Borrows a loaded state as a command session while keeping state ownership explicit.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/lib.rs",
    source_id = "gammaloop_api::LoadedState::cli_session"
)]
fn gammaloop_api_loaded_state_cli_session() {}

#[alphal00p_docs::ty(
    id = "CLISettings",
    title = "Command-line settings",
    summary = "Owns the generated command-line, global, session, and state-loading settings.",
    kind = "setting",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/lib.rs",
    source_id = "gammaloop_api::CLISettings"
)]
fn gammaloop_api_cli_settings() {}

#[alphal00p_docs::ty(
    id = "gammaloop",
    name = "gammaloop",
    title = "GammaLoop command",
    summary = "Runs the Clap-derived command tree for state creation, inspection, generation, and evaluation.",
    kind = "command",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/lib.rs",
    source_id = "gammaloop_api::CLISettings"
)]
fn gammaloop_api_command() {}

#[alphal00p_docs::ty(
    id = "HedgeGraph",
    title = "Half-edge graph",
    summary = "Stores graph topology as paired half-edges with independent edge, vertex, and half-edge payloads.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge.rs",
    source_id = "linnet::half_edge::HedgeGraph"
)]
fn linnet_hedge_graph() {}

#[alphal00p_docs::ty(
    id = "HedgeGraphBuilder",
    title = "Graph builder",
    summary = "Builds a valid half-edge graph while assigning typed node, edge, and hedge identifiers.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge/builder.rs",
    source_id = "linnet::half_edge::builder::HedgeGraphBuilder"
)]
fn linnet_hedge_graph_builder() {}

#[alphal00p_docs::ty(
    id = "SuBitGraph",
    title = "Bit-vector subgraph",
    summary = "Represents a subgraph as an efficient half-edge inclusion filter.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge/subgraph/subset.rs",
    source_id = "linnet::half_edge::subgraph::subset::SuBitGraph"
)]
fn linnet_su_bit_graph() {}

#[alphal00p_docs::ty(
    id = "SimpleTraversalTree",
    title = "Traversal tree",
    summary = "Records parent, root, and traversal information produced by breadth-first or depth-first graph walks.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge/tree.rs",
    source_id = "linnet::half_edge::tree::SimpleTraversalTree"
)]
fn linnet_simple_traversal_tree() {}

#[alphal00p_docs::ty(
    id = "DotGraph",
    title = "DOT-backed graph",
    summary = "Combines a half-edge graph with DOT attributes and parser/serializer support.",
    format = "rust-markdown",
    source = "crates/linnet/src/parser/mod.rs",
    source_id = "linnet::parser::DotGraph"
)]
fn linnet_dot_graph() {}

#[alphal00p_docs::trait_item(
    id = "TensorStructure",
    title = "Tensor structure",
    summary = "Describes ordered slots, dimensions, representations, and contraction compatibility independently of data.",
    format = "rust-markdown",
    source = "crates/spenso/src/structure.rs",
    source_id = "spenso::structure::TensorStructure"
)]
fn spenso_tensor_structure() {}

#[alphal00p_docs::ty(
    id = "DenseTensor",
    title = "Dense tensor",
    summary = "Stores every tensor component in structure-defined flat order.",
    format = "rust-markdown",
    source = "crates/spenso/src/tensors/data/dense.rs",
    source_id = "spenso::tensors::data::dense::DenseTensor"
)]
fn spenso_dense_tensor() {}

#[alphal00p_docs::ty(
    id = "SparseTensor",
    title = "Sparse tensor",
    summary = "Stores only explicit tensor components keyed by flat index.",
    format = "rust-markdown",
    source = "crates/spenso/src/tensors/data/sparse.rs",
    source_id = "spenso::tensors::data::sparse::SparseTensor"
)]
fn spenso_sparse_tensor() {}

#[alphal00p_docs::trait_item(
    id = "Contract",
    title = "Pairwise contraction",
    summary = "Defines pairwise tensor contraction with an explicit settings type and output.",
    format = "rust-markdown",
    source = "crates/spenso/src/contraction.rs",
    source_id = "spenso::contraction::Contract"
)]
fn spenso_contract() {}

#[alphal00p_docs::ty(
    id = "Network",
    title = "Tensor network",
    summary = "Binds symbolic tensor expressions to libraries and an executable contraction graph.",
    format = "rust-markdown",
    source = "crates/spenso/src/network/mod.rs",
    source_id = "spenso::network::Network"
)]
fn spenso_network() {}

#[alphal00p_docs::ty(
    id = "ParamTensor",
    title = "Parametric tensor",
    summary = "Stores tensor components as parameter-dependent symbolic coefficients.",
    format = "rust-markdown",
    source = "crates/spenso/src/tensors/parametric.rs",
    source_id = "spenso::tensors::parametric::ParamTensor"
)]
fn spenso_param_tensor() {}

#[alphal00p_docs::macro_item(
    id = "SimpleRepresentation",
    title = "SimpleRepresentation derive",
    summary = "Derives representation naming and duality boilerplate from the representation helper attribute.",
    format = "rust-markdown",
    source = "crates/spenso-macros/src/lib.rs",
    source_id = "spenso_macros::SimpleRepresentation"
)]
fn spenso_macros_simple_representation() {}

#[alphal00p_docs::func(
    id = "hep_lib",
    title = "HEP tensor library",
    summary = "Constructs the standard high-energy-physics tensor library for a chosen scalar data type.",
    format = "rust-markdown",
    source = "crates/spenso-hep-lib/src/lib.rs",
    source_id = "spenso_hep_lib::hep_lib"
)]
fn spenso_hep_lib_hep_lib() {}

#[alphal00p_docs::func(
    id = "gamma_data_dirac",
    title = "Dirac gamma data",
    summary = "Creates sparse Dirac-basis gamma-matrix data for the supplied tensor structure.",
    format = "rust-markdown",
    source = "crates/spenso-hep-lib/src/lib.rs",
    source_id = "spenso_hep_lib::gamma_data_dirac"
)]
fn spenso_hep_lib_gamma_data_dirac() {}

#[alphal00p_docs::func(
    id = "su3_generator_data",
    title = "SU(3) generator data",
    summary = "Creates the standard sparse SU(3) fundamental-generator tensor.",
    format = "rust-markdown",
    source = "crates/spenso-hep-lib/src/lib.rs",
    source_id = "spenso_hep_lib::su3_generator_data"
)]
fn spenso_hep_lib_su3_generator_data() {}

#[alphal00p_docs::trait_item(
    id = "IndexTooling",
    title = "Index tooling",
    summary = "Adds canonicalization, wrapping, conjugation, and dangling-index inspection to Symbolica atoms.",
    format = "rust-markdown",
    source = "crates/idenso/src/lib.rs",
    source_id = "idenso::IndexTooling"
)]
fn idenso_index_tooling() {}

#[alphal00p_docs::ty(
    id = "CookSettings",
    title = "Index cooking settings",
    summary = "Controls reversible, flattened, and filtered encodings of tensor indices.",
    format = "rust-markdown",
    source = "crates/idenso/src/cook.rs",
    source_id = "idenso::cook::CookSettings"
)]
fn idenso_cook_settings() {}

#[alphal00p_docs::trait_item(
    id = "Cookable",
    title = "Cookable expressions",
    summary = "Applies or reverses index-cooking transformations on symbolic expressions.",
    format = "rust-markdown",
    source = "crates/idenso/src/cook.rs",
    source_id = "idenso::cook::Cookable"
)]
fn idenso_cookable() {}

#[alphal00p_docs::trait_item(
    id = "SelectiveExpand",
    title = "Selective expansion",
    summary = "Expands only tensor families selected by their representation.",
    format = "rust-markdown",
    source = "crates/idenso/src/selective_expand.rs",
    source_id = "idenso::selective_expand::SelectiveExpand"
)]
fn idenso_selective_expand() {}

#[alphal00p_docs::macro_item(
    id = "bis",
    name = "bis!",
    title = "Bispinor representation macro",
    summary = "Builds stripped or indexed bispinor representation atoms in Spenso syntax.",
    format = "rust-markdown",
    source = "crates/idenso/src/lib.rs",
    source_id = "idenso::bis"
)]
fn idenso_bis() {}

#[alphal00p_docs::macro_item(
    id = "cof",
    name = "cof!",
    title = "Color-fundamental macro",
    summary = "Builds stripped or indexed color-fundamental representation atoms.",
    format = "rust-markdown",
    source = "crates/idenso/src/lib.rs",
    source_id = "idenso::cof"
)]
fn idenso_cof() {}

#[alphal00p_docs::macro_item(
    id = "coad",
    name = "coad!",
    title = "Color-adjoint representation macro",
    summary = "Builds stripped or indexed color-adjoint representation atoms.",
    format = "rust-markdown",
    source = "crates/idenso/src/lib.rs",
    source_id = "idenso::coad"
)]
fn idenso_coad() {}

#[alphal00p_docs::macro_item(
    id = "gamma",
    name = "gamma!",
    title = "Dirac gamma macro",
    summary = "Builds a Dirac gamma chain factor or an explicitly indexed gamma tensor.",
    format = "rust-markdown",
    source = "crates/idenso/src/dirac/macros.rs",
    source_id = "idenso::gamma"
)]
fn idenso_gamma() {}

#[alphal00p_docs::macro_item(
    id = "gamma0",
    name = "gamma0!",
    title = "Gamma-zero macro",
    summary = "Builds a gamma-zero chain factor or an explicitly indexed gamma-zero tensor.",
    format = "rust-markdown",
    source = "crates/idenso/src/dirac/macros.rs",
    source_id = "idenso::gamma0"
)]
fn idenso_gamma0() {}

#[alphal00p_docs::macro_item(
    id = "gamma5",
    name = "gamma5!",
    title = "Gamma-five macro",
    summary = "Builds a gamma-five chain factor or an explicitly indexed gamma-five tensor.",
    format = "rust-markdown",
    source = "crates/idenso/src/dirac/macros.rs",
    source_id = "idenso::gamma5"
)]
fn idenso_gamma5() {}

#[alphal00p_docs::macro_item(
    id = "epsilon",
    name = "epsilon!",
    title = "Levi-Civita tensor macro",
    summary = "Builds an arbitrary-rank antisymmetric epsilon tensor from atom-like arguments.",
    format = "rust-markdown",
    source = "crates/idenso/src/epsilon.rs",
    source_id = "idenso::epsilon"
)]
fn idenso_epsilon() {}

#[alphal00p_docs::macro_item(
    id = "color_t",
    name = "color_t!",
    title = "Color-generator macro",
    summary = "Builds a color-generator chain factor or an explicitly indexed generator tensor.",
    format = "rust-markdown",
    source = "crates/idenso/src/color/macros.rs",
    source_id = "idenso::color_t"
)]
fn idenso_color_t() {}

#[alphal00p_docs::macro_item(
    id = "color_f",
    name = "color_f!",
    title = "Color structure-constant macro",
    summary = "Builds an explicitly indexed adjoint color structure-constant tensor.",
    format = "rust-markdown",
    source = "crates/idenso/src/color/macros.rs",
    source_id = "idenso::color_f"
)]
fn idenso_color_f() {}

#[alphal00p_docs::func(
    id = "representations::initialize",
    title = "Initialize representations",
    summary = "Initializes Idenso's standard representations and tensor symbols.",
    format = "rust-markdown",
    source = "crates/idenso/src/representations.rs",
    source_id = "idenso::representations::initialize"
)]
fn idenso_initialize() {}

#[alphal00p_docs::macro_item(
    id = "vakint_parse",
    name = "vakint_parse!",
    title = "Vakint expression parser",
    summary = "Parses a Symbolica expression using Vakint's namespace by default or a caller-supplied namespace.",
    format = "rust-markdown",
    source = "crates/vakint/src/utils.rs",
    source_id = "vakint::vakint_parse"
)]
fn vakint_parse_macro() {}

#[alphal00p_docs::macro_item(
    id = "vakint_symbol",
    name = "vakint_symbol!",
    title = "Vakint symbol parser",
    summary = "Resolves a symbol name in Vakint's namespace while preserving an explicit namespace.",
    format = "rust-markdown",
    source = "crates/vakint/src/utils.rs",
    source_id = "vakint::vakint_symbol"
)]
fn vakint_symbol_macro() {}

#[alphal00p_docs::ty(
    id = "Vakint",
    title = "Vakint engine",
    summary = "Owns supported topology matching, canonicalization, reduction, and evaluation operations.",
    format = "rust-markdown",
    source = "crates/vakint/src/lib.rs",
    source_id = "vakint::Vakint"
)]
fn vakint_engine() {}

#[alphal00p_docs::ty(
    id = "VakintSettings",
    title = "Vakint settings",
    summary = "Controls symbols, precision, normalization, external tools, and backend evaluation order.",
    kind = "setting",
    format = "rust-markdown",
    source = "crates/vakint/src/lib.rs",
    source_id = "vakint::VakintSettings"
)]
fn vakint_settings() {}

#[alphal00p_docs::ty(
    id = "VakintExpression",
    title = "Vakint expression",
    summary = "Stores a sum of matched vacuum-integral terms with separate numerators.",
    format = "rust-markdown",
    source = "crates/vakint/src/lib.rs",
    source_id = "vakint::VakintExpression"
)]
fn vakint_expression() {}

#[alphal00p_docs::ty(
    id = "EvaluationOrder",
    title = "Evaluation order",
    summary = "Defines the ordered analytic and numerical backends available to an evaluation.",
    format = "rust-markdown",
    source = "crates/vakint/src/lib.rs",
    source_id = "vakint::EvaluationOrder"
)]
fn vakint_evaluation_order() {}

#[alphal00p_docs::func(
    id = "VakintExpression::canonicalize",
    title = "Canonicalize an expression",
    summary = "Matches and rewrites each integral into Vakint's canonical topology representation.",
    format = "rust-markdown",
    source = "crates/vakint/src/lib.rs",
    source_id = "vakint::VakintExpression::canonicalize"
)]
fn vakint_expression_canonicalize() {}

#[alphal00p_docs::func(
    id = "VakintExpression::tensor_reduce",
    title = "Tensor-reduce an expression",
    summary = "Reduces tensor numerators to scalar-integral combinations.",
    format = "rust-markdown",
    source = "crates/vakint/src/lib.rs",
    source_id = "vakint::VakintExpression::tensor_reduce"
)]
fn vakint_expression_tensor_reduce() {}

#[alphal00p_docs::func(
    id = "VakintExpression::evaluate_integral",
    title = "Evaluate an expression",
    summary = "Evaluates canonical scalar integrals using the configured ordered backend policy.",
    format = "rust-markdown",
    source = "crates/vakint/src/lib.rs",
    source_id = "vakint::VakintExpression::evaluate_integral"
)]
fn vakint_expression_evaluate_integral() {}

#[alphal00p_docs::ty(
    id = "NumericalEvaluationResult",
    title = "Numerical Laurent result",
    summary = "Stores numerical Laurent coefficients by power of the dimensional regulator.",
    format = "rust-markdown",
    source = "crates/vakint/src/lib.rs",
    source_id = "vakint::NumericalEvaluationResult"
)]
fn vakint_numerical_evaluation_result() {}

pub(super) fn for_component(component: &str) -> Option<Vec<DocItem>> {
    Some(match component {
        "gammalooprs" => vec![
            __alphal00p_docs_trait_gammalooprs_gamma_loop_context(),
            __alphal00p_docs_trait_gammalooprs_has_model(),
            __alphal00p_docs_func_gammalooprs_set_interrupt_handler(),
        ],
        "gammaloop-api" => vec![
            __alphal00p_docs_ty_gammaloop_api_state_load_option(),
            __alphal00p_docs_func_gammaloop_api_state_load_option_load(),
            __alphal00p_docs_func_gammaloop_api_loaded_state_cli_session(),
            __alphal00p_docs_ty_gammaloop_api_cli_settings(),
            __alphal00p_docs_ty_gammaloop_api_command(),
        ],
        "linnet" => vec![
            __alphal00p_docs_ty_linnet_hedge_graph(),
            __alphal00p_docs_ty_linnet_hedge_graph_builder(),
            __alphal00p_docs_ty_linnet_su_bit_graph(),
            __alphal00p_docs_ty_linnet_simple_traversal_tree(),
            __alphal00p_docs_ty_linnet_dot_graph(),
        ],
        "spenso" => vec![
            __alphal00p_docs_trait_spenso_tensor_structure(),
            __alphal00p_docs_ty_spenso_dense_tensor(),
            __alphal00p_docs_ty_spenso_sparse_tensor(),
            __alphal00p_docs_trait_spenso_contract(),
            __alphal00p_docs_ty_spenso_network(),
            __alphal00p_docs_ty_spenso_param_tensor(),
        ],
        "spenso-macros" => {
            vec![__alphal00p_docs_macro_spenso_macros_simple_representation()]
        }
        "spenso-hep-lib" => vec![
            __alphal00p_docs_func_spenso_hep_lib_hep_lib(),
            __alphal00p_docs_func_spenso_hep_lib_gamma_data_dirac(),
            __alphal00p_docs_func_spenso_hep_lib_su3_generator_data(),
        ],
        "idenso" => vec![
            __alphal00p_docs_trait_idenso_index_tooling(),
            __alphal00p_docs_ty_idenso_cook_settings(),
            __alphal00p_docs_trait_idenso_cookable(),
            __alphal00p_docs_trait_idenso_selective_expand(),
            __alphal00p_docs_macro_idenso_bis(),
            __alphal00p_docs_macro_idenso_cof(),
            __alphal00p_docs_macro_idenso_coad(),
            __alphal00p_docs_macro_idenso_gamma(),
            __alphal00p_docs_macro_idenso_gamma0(),
            __alphal00p_docs_macro_idenso_gamma5(),
            __alphal00p_docs_macro_idenso_epsilon(),
            __alphal00p_docs_macro_idenso_color_t(),
            __alphal00p_docs_macro_idenso_color_f(),
            __alphal00p_docs_func_idenso_initialize(),
        ],
        "vakint" => vec![
            __alphal00p_docs_macro_vakint_parse_macro(),
            __alphal00p_docs_macro_vakint_symbol_macro(),
            __alphal00p_docs_ty_vakint_engine(),
            __alphal00p_docs_ty_vakint_settings(),
            __alphal00p_docs_ty_vakint_expression(),
            __alphal00p_docs_ty_vakint_evaluation_order(),
            __alphal00p_docs_func_vakint_expression_canonicalize(),
            __alphal00p_docs_func_vakint_expression_tensor_reduce(),
            __alphal00p_docs_func_vakint_expression_evaluate_integral(),
            __alphal00p_docs_ty_vakint_numerical_evaluation_result(),
        ],
        _ => return None,
    })
}
