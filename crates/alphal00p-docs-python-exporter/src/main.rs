//! Emit one pyo3-stub-gen inventory as a deterministic checked-in stub.

use std::{env, error::Error, fs, path::PathBuf};

#[cfg(any(feature = "gammaloop", feature = "idenso", feature = "vakint"))]
use pyo3::types::{PyAnyMethods as _, PyDictMethods as _, PyModuleMethods as _};
#[cfg(any(
    feature = "gammaloop",
    feature = "idenso",
    feature = "spenso",
    feature = "vakint"
))]
use std::collections::BTreeSet;

fn main() -> Result<(), Box<dyn Error>> {
    let mut arguments = env::args().skip(1);
    let component = arguments
        .next()
        .ok_or("usage: alphal00p-docs-python-exporter <component> <output.pyi>")?;
    let output = PathBuf::from(
        arguments
            .next()
            .ok_or("usage: alphal00p-docs-python-exporter <component> <output.pyi>")?,
    );
    let check = match arguments.next().as_deref() {
        None => false,
        Some("--check") => true,
        Some(_) => return Err("the optional third argument must be --check".into()),
    };
    if arguments.next().is_some() {
        return Err("expected a component, output, and optional --check".into());
    }

    let (module_name, stub_info) = gather(&component)?;
    let module = stub_info
        .modules
        .get(module_name)
        .ok_or_else(|| format!("{component} inventory has no module {module_name}"))?;
    #[cfg(feature = "gammaloop")]
    if component == "gammaloop-python" {
        validate_gammaloop_stub_surface(module)?;
    }
    #[cfg(feature = "spenso")]
    if component == "spynso3" {
        validate_spenso_stub_surface(module)?;
    }
    #[cfg(feature = "idenso")]
    if component == "idenso-community" {
        validate_idenso_stub_surface(module)?;
    }
    #[cfg(feature = "vakint")]
    if component == "vakint-community" {
        validate_vakint_stub_surface(module)?;
    }
    if let Some(parent) = output.parent() {
        fs::create_dir_all(parent)?;
    }
    let rendered = render(&component, module)?;
    let mut normalized = rendered
        .lines()
        .map(str::trim_end)
        .collect::<Vec<_>>()
        .join("\n");
    normalized.push('\n');
    let mut outputs = vec![output];
    if component == "linnet-py" {
        outputs.push(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../linnet-py/linnet_py.pyi"));
    }
    outputs.sort();
    outputs.dedup();
    for output in outputs {
        if check {
            let checked_in = fs::read_to_string(&output)?;
            if checked_in != normalized {
                let hint = if component == "linnet-py" {
                    "regenerate the shared Linnet package/docs surface"
                } else {
                    "regenerate the checked-in snapshot"
                };
                return Err(format!(
                    "generated Python stub drifted: {}; {hint}",
                    output.display(),
                )
                .into());
            }
        } else {
            if let Some(parent) = output.parent() {
                fs::create_dir_all(parent)?;
            }
            fs::write(output, &normalized)?;
        }
    }
    Ok(())
}

#[cfg(feature = "gammaloop")]
fn validate_gammaloop_stub_surface(
    module: &pyo3_stub_gen::generate::Module,
) -> Result<(), Box<dyn Error>> {
    let runtime = pyo3::Python::attach(|py| {
        let module = pyo3::types::PyModule::new(py, "gammaloop._gammaloop")?;
        gammaloop_api::python::register_python_api_for_docs(&module)?;
        public_module_names(&module)
    })?;
    let stub = stub_names(module);
    if stub == runtime {
        return Ok(());
    }
    Err(format!(
        "GammaLoop StubInfo does not match runtime registration: missing {:?}; unreachable {:?}",
        runtime.difference(&stub).collect::<Vec<_>>(),
        stub.difference(&runtime).collect::<Vec<_>>()
    )
    .into())
}

fn render(
    component: &str,
    module: &pyo3_stub_gen::generate::Module,
) -> Result<String, Box<dyn Error>> {
    match component {
        #[cfg(feature = "linnet")]
        "linnet-py" => Ok(linnet_py::canonical_stub()?),
        _ => Ok(module.to_string()),
    }
}

#[cfg(feature = "spenso")]
fn validate_spenso_stub_surface(
    module: &pyo3_stub_gen::generate::Module,
) -> Result<(), Box<dyn Error>> {
    let surface = &spynso3::PYTHON_STUB_SURFACE;
    let expected = surface
        .registered
        .iter()
        .chain(surface.returned_opaque)
        .copied()
        .collect::<BTreeSet<_>>();
    let actual = stub_names(module);
    let expected = expected
        .into_iter()
        .map(str::to_owned)
        .collect::<BTreeSet<_>>();

    let missing = expected.difference(&actual).cloned().collect::<Vec<_>>();
    let unreachable = actual.difference(&expected).cloned().collect::<Vec<_>>();
    if missing.is_empty() && unreachable.is_empty() {
        return Ok(());
    }

    Err(format!(
        "Spenso StubInfo does not match the community-module surface: missing {missing:?}; unreachable {unreachable:?}"
    )
    .into())
}

#[cfg(feature = "idenso")]
fn validate_idenso_stub_surface(
    module: &pyo3_stub_gen::generate::Module,
) -> Result<(), Box<dyn Error>> {
    use symbolica::api::python::SymbolicaCommunityModule;

    let expected = idenso::python::PYTHON_STUB_SURFACE
        .iter()
        .map(|name| (*name).to_owned())
        .collect::<BTreeSet<_>>();
    let runtime = pyo3::Python::attach(|py| {
        let module = pyo3::types::PyModule::new(py, "symbolica.community.idenso")?;
        idenso::python::IdensoModule::register_module(&module)?;
        public_module_names(&module)
    })?;
    validate_exact_surface("Idenso", &expected, &stub_names(module), &runtime)
}

#[cfg(feature = "vakint")]
fn validate_vakint_stub_surface(
    module: &pyo3_stub_gen::generate::Module,
) -> Result<(), Box<dyn Error>> {
    use symbolica::api::python::SymbolicaCommunityModule;

    let expected = vakint::symbolica_community_module::PYTHON_STUB_SURFACE
        .iter()
        .map(|name| (*name).to_owned())
        .collect::<BTreeSet<_>>();
    let runtime = pyo3::Python::attach(|py| {
        let module = pyo3::types::PyModule::new(py, "symbolica.community.vakint")?;
        vakint::symbolica_community_module::VakintWrapper::register_module(&module)?;
        public_module_names(&module)
    })?;
    validate_exact_surface("Vakint", &expected, &stub_names(module), &runtime)
}

#[cfg(any(
    feature = "gammaloop",
    feature = "idenso",
    feature = "spenso",
    feature = "vakint"
))]
fn stub_names(module: &pyo3_stub_gen::generate::Module) -> BTreeSet<String> {
    module
        .class
        .values()
        .map(|class| class.name.to_owned())
        .chain(module.enum_.values().map(|enum_| enum_.name.to_owned()))
        .chain(module.function.keys().map(|name| (*name).to_owned()))
        .chain(module.variables.keys().map(|name| (*name).to_owned()))
        .collect()
}

#[cfg(any(feature = "gammaloop", feature = "idenso", feature = "vakint"))]
fn public_module_names(
    module: &pyo3::Bound<'_, pyo3::types::PyModule>,
) -> pyo3::PyResult<BTreeSet<String>> {
    module
        .dict()
        .iter()
        .filter_map(|(name, _)| name.extract::<String>().ok())
        .filter(|name| !name.starts_with('_'))
        .map(Ok)
        .collect()
}

#[cfg(any(feature = "idenso", feature = "vakint"))]
fn validate_exact_surface(
    product: &str,
    expected: &BTreeSet<String>,
    stub: &BTreeSet<String>,
    runtime: &BTreeSet<String>,
) -> Result<(), Box<dyn Error>> {
    if stub == expected && runtime == expected {
        return Ok(());
    }
    Err(format!(
        "{product} Python surface mismatch: declared {expected:?}; StubInfo {stub:?}; registered {runtime:?}"
    )
    .into())
}

fn gather(component: &str) -> Result<(&'static str, pyo3_stub_gen::StubInfo), Box<dyn Error>> {
    match component {
        #[cfg(feature = "gammaloop")]
        "gammaloop-python" => Ok(("gammaloop._gammaloop", gammaloop_api::python::stub_info()?)),
        #[cfg(feature = "linnet")]
        "linnet-py" => Ok(("linnet_py", linnet_py::stub_info()?)),
        #[cfg(feature = "spenso")]
        "spynso3" => Ok(("symbolica.community.spenso", spynso3::stub_info()?)),
        #[cfg(feature = "idenso")]
        "idenso-community" => Ok(("symbolica.community.idenso", idenso::stub_info()?)),
        #[cfg(feature = "vakint")]
        "vakint-community" => Ok(("symbolica.community.vakint", vakint::stub_info()?)),
        _ => Err(format!(
            "component {component} is not enabled; select its matching Cargo feature"
        )
        .into()),
    }
}

#[cfg(test)]
mod tests {
    #[cfg(feature = "gammaloop")]
    use std::ffi::CString;

    #[cfg(feature = "gammaloop")]
    const GAMMALOOP_STUB: &str = include_str!("../../../docs/api/python/gammaloop-python.pyi");

    #[cfg(feature = "gammaloop")]
    const GAMMALOOP_RUNTIME_VALIDATOR: &str = r#"
import ast
import inspect

_EMPTY = inspect.Parameter.empty
_PROTOCOL_SLOT_ONLY = {"__getattr__"}


def _fail(message):
    raise AssertionError(message)


def _stub_parameters(function, drop_receiver=False):
    arguments = function.args
    positional = [
        (argument, inspect.Parameter.POSITIONAL_ONLY)
        for argument in arguments.posonlyargs
    ] + [
        (argument, inspect.Parameter.POSITIONAL_OR_KEYWORD)
        for argument in arguments.args
    ]
    defaults = [None] * (len(positional) - len(arguments.defaults)) + list(arguments.defaults)
    parameters = []
    for (argument, kind), default in zip(positional, defaults):
        parameters.append((argument.arg, kind, default))
    if arguments.vararg is not None:
        parameters.append((arguments.vararg.arg, inspect.Parameter.VAR_POSITIONAL, None))
    parameters.extend(
        (argument.arg, inspect.Parameter.KEYWORD_ONLY, default)
        for argument, default in zip(arguments.kwonlyargs, arguments.kw_defaults)
    )
    if arguments.kwarg is not None:
        parameters.append((arguments.kwarg.arg, inspect.Parameter.VAR_KEYWORD, None))
    if drop_receiver and parameters and parameters[0][0] in {"self", "cls"}:
        parameters.pop(0)
    return parameters


def _runtime_parameters(callable_, drop_receiver=False):
    try:
        parameters = list(inspect.signature(callable_).parameters.values())
    except (TypeError, ValueError) as error:
        _fail(f"cannot inspect {callable_!r}: {error}")
    if drop_receiver and parameters and parameters[0].name in {"self", "cls"}:
        parameters.pop(0)
    return parameters


def _literal_default(node):
    if node is None:
        return _EMPTY
    if isinstance(node, ast.Constant) and node.value is Ellipsis:
        return Ellipsis
    try:
        return ast.literal_eval(node)
    except (TypeError, ValueError):
        return Ellipsis


def _check_signature(label, stub_function, runtime_callable, drop_receiver=False):
    expected = _stub_parameters(stub_function, drop_receiver)
    actual = _runtime_parameters(runtime_callable, drop_receiver)
    expected_shape = [(name, kind) for name, kind, _ in expected]
    actual_shape = [(parameter.name, parameter.kind) for parameter in actual]
    if actual_shape != expected_shape:
        _fail(f"{label} parameters differ: stub {expected_shape!r}; runtime {actual_shape!r}")
    for (name, _, stub_default), parameter in zip(expected, actual):
        expected_default = _literal_default(stub_default)
        runtime_default = parameter.default
        if (expected_default is _EMPTY) != (runtime_default is _EMPTY):
            _fail(
                f"{label}.{name} required/default status differs: "
                f"stub {expected_default!r}; runtime {runtime_default!r}"
            )
        if (
            expected_default is not _EMPTY
            and expected_default is not Ellipsis
            and runtime_default is not Ellipsis
            and runtime_default != expected_default
        ):
            _fail(
                f"{label}.{name} default differs: "
                f"stub {expected_default!r}; runtime {runtime_default!r}"
            )


def _stub_exports(tree):
    for statement in tree.body:
        if isinstance(statement, (ast.Assign, ast.AnnAssign)):
            targets = statement.targets if isinstance(statement, ast.Assign) else [statement.target]
            if any(isinstance(target, ast.Name) and target.id == "__all__" for target in targets):
                return set(ast.literal_eval(statement.value))
    return {
        statement.name
        for statement in tree.body
        if isinstance(statement, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
    }


def _check_class(stub_class, runtime_class):
    expected_members = {
        statement.target.id
        for statement in stub_class.body
        if isinstance(statement, ast.AnnAssign) and isinstance(statement.target, ast.Name)
    }
    expected_members.update(
        target.id
        for statement in stub_class.body
        if isinstance(statement, ast.Assign)
        for target in statement.targets
        if isinstance(target, ast.Name)
    )
    callables = {}
    for statement in stub_class.body:
        if not isinstance(statement, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        if any(
            isinstance(decorator, ast.Attribute) and decorator.attr == "setter"
            for decorator in statement.decorator_list
        ):
            expected_members.add(statement.name)
            continue
        callables[statement.name] = statement
    expected_members.update(callables)
    expected_runtime_members = expected_members - _PROTOCOL_SLOT_ONLY
    runtime_members = {
        name
        for name in vars(runtime_class)
        if not name.startswith("_") or name in expected_runtime_members
    }
    if runtime_members != expected_runtime_members:
        _fail(
            f"{stub_class.name} members differ: "
            f"missing {sorted(expected_runtime_members - runtime_members)!r}; "
            f"undocumented {sorted(runtime_members - expected_runtime_members)!r}"
        )
    for name, stub_function in callables.items():
        if name in _PROTOCOL_SLOT_ONLY:
            continue
        decorators = {
            decorator.id
            for decorator in stub_function.decorator_list
            if isinstance(decorator, ast.Name)
        }
        if "property" in decorators:
            if callable(getattr(runtime_class, name)):
                _fail(f"{stub_class.name}.{name} is documented as a property but is callable")
            continue
        if name == "__new__":
            runtime_callable = runtime_class
        else:
            runtime_callable = getattr(runtime_class, name)
        _check_signature(
            f"{stub_class.name}.{name}",
            stub_function,
            runtime_callable,
            drop_receiver=True,
        )


def validate(runtime_module, stub_source):
    tree = ast.parse(stub_source, filename="gammaloop-python.pyi", mode="exec")
    declarations = {
        statement.name: statement
        for statement in tree.body
        if isinstance(statement, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
    }
    expected_exports = _stub_exports(tree)
    runtime_exports = {name for name in vars(runtime_module) if not name.startswith("_")}
    if runtime_exports != expected_exports:
        _fail(
            "GammaLoop module exports differ: "
            f"missing {sorted(expected_exports - runtime_exports)!r}; "
            f"undocumented {sorted(runtime_exports - expected_exports)!r}"
        )
    if set(declarations) != expected_exports:
        _fail(
            "GammaLoop stub declarations and __all__ differ: "
            f"missing {sorted(expected_exports - set(declarations))!r}; "
            f"unexported {sorted(set(declarations) - expected_exports)!r}"
        )
    for name in sorted(expected_exports):
        stub_declaration = declarations[name]
        runtime_declaration = getattr(runtime_module, name)
        if isinstance(stub_declaration, ast.ClassDef):
            if not isinstance(runtime_declaration, type):
                _fail(f"{name} is documented as a class but runtime exposes {runtime_declaration!r}")
            _check_class(stub_declaration, runtime_declaration)
        else:
            if not callable(runtime_declaration):
                _fail(f"{name} is documented as callable but runtime exposes {runtime_declaration!r}")
            _check_signature(name, stub_declaration, runtime_declaration)
"#;

    #[cfg(feature = "linnet")]
    #[test]
    fn linnet_package_and_docs_share_the_typed_stub_info_surface() {
        let canonical = linnet_py::canonical_stub().unwrap();
        assert_eq!(canonical, include_str!("../../linnet-py/linnet_py.pyi"));
        assert_eq!(
            canonical,
            include_str!("../../../docs/api/python/linnet-py.pyi")
        );
        for declaration in [
            "_NativeValue: typing.TypeAlias = (",
            "def __new__(cls, *, encode_node: typing.Callable[[NodeValue], DotVertexData]",
            "def build(*items: _GraphItem, name: _OptionalString = None",
            "def map(self, *, node: typing.Callable[[Node], typing.Any] | None = None",
            "def render(self, output: builtins.str | os.PathLike[builtins.str], *, config: RenderConfig | None = None)",
            "node_selector: _NodeSelector = ...",
            "node_outset: _AutoNumber = ...",
            "_Padding: typing.TypeAlias = builtins.int | builtins.float | _NativeArray | _NativeDict | Insets",
            "show_node_index: _Boolean = ...",
            "def call(self, *args: _NativeValue, **kwargs: _NativeValue) -> TypstCall",
        ] {
            assert!(canonical.contains(declaration), "missing `{declaration}`");
        }
        for legacy_or_broad in ["**kwargs: typing.Any", "**drawing", "slf:"] {
            assert!(
                !canonical.contains(legacy_or_broad),
                "unexpected `{legacy_or_broad}`"
            );
        }
        assert_eq!(canonical.matches("def anchor(self)").count(), 1);
    }

    #[cfg(feature = "gammaloop")]
    #[test]
    fn gammaloop_stubinfo_tracks_runtime_registration() {
        let (module_name, stub_info) =
            super::gather("gammaloop-python").expect("GammaLoop StubInfo");
        let module = stub_info
            .modules
            .get(module_name)
            .expect("GammaLoop extension module");
        super::validate_gammaloop_stub_surface(module)
            .expect("StubInfo and runtime registration match");

        let mut rendered = module
            .to_string()
            .lines()
            .map(str::trim_end)
            .collect::<Vec<_>>()
            .join("\n");
        rendered.push('\n');
        assert_eq!(rendered, GAMMALOOP_STUB);
    }

    #[cfg(feature = "gammaloop")]
    #[test]
    fn gammaloop_runtime_surface_and_signatures_match_the_docs_stub() -> pyo3::PyResult<()> {
        use pyo3::{prelude::*, types::PyModule};

        Python::attach(|py| {
            let runtime = PyModule::new(py, "gammaloop._gammaloop")?;
            gammaloop_api::python::register_python_api_for_docs(&runtime)?;
            let validator = PyModule::from_code(
                py,
                CString::new(GAMMALOOP_RUNTIME_VALIDATOR)
                    .expect("validator has no NUL bytes")
                    .as_c_str(),
                c"gammaloop_docs_runtime_validator.py",
                c"gammaloop_docs_runtime_validator",
            )?;
            validator
                .getattr("validate")?
                .call1((runtime, GAMMALOOP_STUB))?;
            Ok(())
        })
    }

    #[cfg(feature = "gammaloop")]
    #[test]
    fn gammaloop_stub_preserves_supported_reference_prose() {
        let (module_name, stub_info) =
            super::gather("gammaloop-python").expect("GammaLoop StubInfo");
        let rendered = stub_info
            .modules
            .get(module_name)
            .expect("GammaLoop extension module")
            .to_string()
            .lines()
            .map(str::trim_end)
            .collect::<Vec<_>>()
            .join("\n");

        for excerpt in [
            "def e(self) -> builtins.float:\n        r\"\"\"\n        Energy component.",
            "def graph_groups(self) -> builtins.list[IntegrandGraphGroup]:\n        r\"\"\"\n        Per-group graph, orientation, basis, threshold, and cut metadata.",
            "def sample(self) -> SampleEvaluationResult:\n        r\"\"\"\n        Evaluation record for the requested sample.",
            "def combine_diagrams(self, value: builtins.bool) -> None:\n        r\"\"\"\n        Write all diagrams of an integrand to one file during filesystem export.",
            "Create an evenly binned continuous histogram accumulator.\n\n        Parameters\n        ----------\n        title : str",
            "Resolve a relative string path or array index.\n\n        Parameters\n        ----------\n        key : str or int",
            "Return the causal-flow orientations generated for one graph.\n\n        Each returned dictionary maps an edge id to ``1`` (default), ``-1``\n        (reversed), or ``0`` (undirected). Supply process and integrand selectors when\n        the active state does not identify a unique integrand.\n\n        Parameters\n        ----------\n        graph_name : str",
        ] {
            assert!(rendered.contains(excerpt), "missing `{excerpt}`");
        }
    }

    #[cfg(feature = "spenso")]
    #[test]
    fn spenso_stub_tracks_registered_and_returned_opaque_types() {
        let (module_name, stub_info) = super::gather("spynso3").expect("Spenso StubInfo");
        let module = stub_info
            .modules
            .get(module_name)
            .expect("Spenso community module");

        super::validate_spenso_stub_surface(module).expect("matching Spenso Python surface");
        let rendered = module.to_string();
        assert!(!rendered.contains("TensorFunctionLibrary"));
        assert!(!rendered.contains("TensorNamespace"));
    }

    #[cfg(feature = "spenso")]
    #[test]
    fn spenso_stub_preserves_optional_structure_defaults_and_execution_semantics() {
        let (module_name, stub_info) = super::gather("spynso3").expect("Spenso StubInfo");
        let module = stub_info
            .modules
            .get(module_name)
            .expect("Spenso community module");
        let rendered = module.to_string();

        for signature in [
            "def __call__(self, *args: builtins.int | Expression | str, extra_args: typing.Sequence[Expression | int | str | float | builtins.complex] | None = None)",
            "def symbolic(self, *args: builtins.int | Expression | str, extra_args: typing.Sequence[Expression | int | str | float | builtins.complex] | None = None)",
            "def index(self, *args: builtins.int | Expression | str, extra_args: typing.Sequence[Expression] | None = None, cook_indices: builtins.bool = False)",
        ] {
            assert!(rendered.contains(signature), "missing `{signature}`");
        }
        assert!(rendered.contains(
            "Single : Select one smallest-degree rewrite per step; without `n_steps`, continue until no work remains"
        ));
        assert!(!rendered.contains("Single : Execute one contraction at a time"));
    }
}
