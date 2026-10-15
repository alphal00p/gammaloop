use bincode_trait_derive::{Decode, Encode};
use eyre::{Result, eyre};
use serde::{Deserialize, Serialize};
use symbolica::{
    evaluate::{ExpressionEvaluator, Instruction, JITCompilationSettings, JITCompiledEvaluator},
    prelude::{Complex, Rational},
};

/// Compile the existing complex program without exceeding SymJIT's call scratch space.
/// Symbolica 3.0.1 does not reserve slots for retained-call arguments; SymJIT 2.27
/// defaults to 64 and silently clamps reservations to 2048. Its direct translator
/// also accepts oversized calls without validating their stores. Until upstream
/// fixes both boundaries, size the reservation from the complete instruction DAG
/// and reject unsupported calls before either translator runs. Hoisted constants
/// are actual call arguments and must be counted, unlike root input parameters.
pub(crate) fn compile_symjit(
    evaluator: &ExpressionEvaluator<Complex<f64>>,
    settings: JITCompilationSettings,
) -> Result<JITCompiledEvaluator<Complex<f64>>> {
    let instructions = evaluator.export_instructions();
    let mut pending = vec![&instructions];
    let mut visited = std::collections::HashSet::new();
    let mut capacity = 64;
    while let Some(program) = pending.pop() {
        for instruction in &program.instructions {
            if let Instruction::Fun(_, function, _) = instruction {
                capacity = capacity.max(function.2.len());
            }
        }
        for callee in &program.sub_evaluators {
            if visited.insert(std::sync::Arc::as_ptr(&callee.instructions)) {
                pending.push(callee.instructions.as_ref());
            }
        }
    }
    if capacity > 2048 {
        return Err(eyre!(
            "SymJIT 2.27 supports at most 2048 call arguments, but this evaluator requires {capacity} including hoisted constants; refusing unsafe JIT compilation"
        ));
    }
    evaluator
        .jit_compile(settings.with_option("num_args", capacity.to_string()))
        .map_err(|error| eyre!(error))
}

/// Physical inputs precede any program-specific row or auxiliary ports. Their
/// archived constants are parameter-major jets in the same shape for every layout.
/// Cross-section archives may omit representative inputs and require caller inputs.
pub(crate) fn validate_standalone_input_layout<'a>(
    parameter_count: usize,
    input_count: Option<usize>,
    mult_offset: usize,
    starts: [usize; 2],
    dual_shapes: impl IntoIterator<Item = Option<&'a [Vec<usize>]>>,
) -> Result<()> {
    let mut shapes = dual_shapes.into_iter();
    let shape = shapes
        .next()
        .ok_or_else(|| eyre!("Standalone stack has no evaluator"))?;
    let width = shape.map_or(1, <[Vec<usize>]>::len);
    if width == 0
        || mult_offset != width
        || input_count.is_some_and(|count| parameter_count.checked_mul(width) != Some(count))
        // Empty parameter groups can start immediately after the last input.
        || starts.into_iter().any(|start| start > parameter_count)
        || shapes.any(|other| other != shape)
    {
        return Err(eyre!(
            "Standalone stack has inconsistent physical input or dual layout"
        ));
    }
    Ok(())
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneParametricResidueRows<A> {
    pub(crate) parameters: Vec<A>,
    pub(crate) rows: Vec<Vec<Rational>>,
    /// Reachable native inputs, indexed into the unchanged coefficient catalog.
    pub(crate) input_columns: Vec<usize>,
    pub(crate) auxiliary_parameters: Vec<A>,
    pub(crate) auxiliary_values: Vec<Rational>,
    pub(crate) auxiliary_shared_outputs: usize,
    /// Dynamic carrier ports gather complete output jets; None supplies zero.
    pub(crate) auxiliary_gathers: Vec<Vec<Option<usize>>>,
}

impl<A> StandaloneParametricResidueRows<A> {
    pub(crate) fn validate(&self, auxiliary_layout: Option<(usize, usize)>) -> Result<()> {
        if self.rows.is_empty() {
            return Err(eyre!(
                "Standalone parametric residue catalog must contain at least one row"
            ));
        }
        if self
            .rows
            .iter()
            .any(|row| row.len() != self.parameters.len())
        {
            return Err(eyre!(
                "Standalone parametric residue row width differs from its parameter count"
            ));
        }
        if self
            .input_columns
            .iter()
            .any(|&index| index >= self.parameters.len())
            || self.input_columns.windows(2).any(|pair| pair[0] >= pair[1])
        {
            return Err(eyre!(
                "Standalone coefficient input columns must be ordered, unique, and in range"
            ));
        }
        if self.auxiliary_parameters.len() != self.auxiliary_values.len()
            || self.auxiliary_gathers.len() != self.rows.len()
        {
            return Err(eyre!(
                "Standalone auxiliary coefficient or gather counts are inconsistent"
            ));
        }
        let (ports, outputs) = auxiliary_layout.unwrap_or((0, 0));
        if self.auxiliary_shared_outputs > ports
            || self.auxiliary_shared_outputs > outputs
            || (auxiliary_layout.is_none() && !self.auxiliary_parameters.is_empty())
        {
            return Err(eyre!(
                "Standalone auxiliary prefix or coefficient inputs require a matching payload"
            ));
        }
        let width = ports - self.auxiliary_shared_outputs;
        if self
            .auxiliary_gathers
            .iter()
            .any(|row| row.len() != width || row.iter().flatten().any(|&index| index >= outputs))
        {
            return Err(eyre!(
                "Standalone auxiliary gather dimensions or output indices are invalid"
            ));
        }
        Ok(())
    }
}
