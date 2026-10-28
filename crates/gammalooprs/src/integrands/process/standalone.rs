use eyre::{Result, eyre};

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
