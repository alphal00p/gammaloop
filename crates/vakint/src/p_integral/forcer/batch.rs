//! Independent exact requests share one linear FORCER invocation. Only raw
//! scalar channels enter FORM; their signs, scales and spectators stay owned
//! by the original request and its existing cache key.

use super::*;

const MAX_REQUESTS: usize = 32;
const MAX_CHANNELS: usize = 256;
const MAX_PROGRAM_BYTES: usize = 1 << 20;

fn combined_program(inputs: &[&ForcerInput]) -> Result<(String, Vec<String>), VakintError> {
    let mut names = Vec::new();
    let mut locals = Vec::new();
    let mut substitutions = vec![".sort".to_owned()];
    for (request, input) in inputs.iter().enumerate() {
        let mut request_names = Vec::new();
        for (part, polynomial) in input.polynomials.iter().enumerate() {
            let name = format!("B{request}F{part}");
            locals.push(format!("Local {name} = {polynomial}*{};", input.topology));
            request_names.push(name.clone());
            names.push(name);
        }
        substitutions.extend([
            "Skip;".to_owned(),
            format!("NSkip {};", request_names.join(",")),
            input.substitutions.clone(),
            ".sort".to_owned(),
        ]);
    }
    let channels = (0..names.len())
        .map(|index| format!("vakintchannel{index}"))
        .collect::<Vec<_>>();
    locals.insert(0, format!("S {};", channels.join(",")));
    if inputs.iter().any(|input| {
        input
            .polynomials
            .iter()
            .any(|polynomial| polynomial.contains("rat("))
    }) {
        locals.insert(0, "PolyRatFun rat;".to_owned());
    }
    substitutions.push(format!("Drop {};", names.join(",")));
    substitutions.push(format!(
        "Local VakintChannels = {};",
        channels
            .iter()
            .zip(&names)
            .map(|(channel, name)| format!("{channel}*{name}"))
            .collect::<Vec<_>>()
            .join("+")
    ));
    let mut writes = vec![format!(
        "Bracket {};\n.sort\nDrop VakintChannels;",
        channels.join(",")
    )];
    writes.extend(
        names
            .iter()
            .zip(&channels)
            .map(|(name, channel)| format!("Local {name} = VakintChannels[{channel}];")),
    );
    writes.push(".sort".into());
    writes.extend(
        names
            .iter()
            .map(|name| format!("#write <{name}.out> \"%E\", {name}")),
    );
    let loops = inputs.iter().map(|input| input.loops).max().unwrap_or(0);
    let template = Template::parse_template(TEMPLATES.get("run_forcer.txt").unwrap()).unwrap();
    let program = template
        .render(&RenderOptions {
            variables: HashMap::from([
                (
                    "vectors".into(),
                    (1..=loops)
                        .map(|index| format!("k{index}"))
                        .collect::<Vec<_>>()
                        .join(","),
                ),
                ("locals".into(), locals.join("\n")),
                ("substitutions".into(), substitutions.join("\n")),
                ("writes".into(), writes.join("\n")),
            ]),
            ..Default::default()
        })
        .map_err(|error| VakintError::InvalidGenericExpression(error.to_string()))?;
    Ok((program, names))
}

/// Bounded batches of cache misses; completed and pending requests retain the
/// same individual keys as `reduce_p`, independently of order or batch peers.
pub(crate) fn reduce_inputs(
    inputs: &[ForcerInput],
    settings: &VakintSettings,
) -> Result<Vec<Vec<ReducedPropagatorValue>>, VakintError> {
    let _stage = StageSpan::new("p_integral.forcer.batch_dispatch");
    let started = Instant::now();
    if inputs.iter().any(|input| !input.parts.is_empty()) {
        crate::citations::USED_CITATIONS.record(crate::citations::CitationSource::Form);
        crate::citations::USED_CITATIONS.record(crate::citations::CitationSource::Forcer);
    }
    let keys = (settings.forcer_cache && settings.clean_tmp_dir).then(|| {
        let _stage = StageSpan::new("p_integral.forcer.cache_key");
        inputs
            .iter()
            .map(|input| {
                let names = input
                    .parts
                    .iter()
                    .map(|part| part.name.as_str())
                    .collect::<Vec<_>>();
                cache::ForcerRunKey::new(&form_command(settings), &input.program, &names)
            })
            .collect::<Vec<_>>()
    });
    let mut jobs = 0_usize;
    let mut misses = 0_usize;
    let run_started = Instant::now();
    let mut compute =
        |missing: &[usize],
         publish: &mut dyn FnMut(usize, Vec<String>) -> Result<(), VakintError>| {
            let _stage = StageSpan::new("p_integral.forcer.batch_prepare");
            misses += missing.len();
            let mut begin = 0;
            while begin < missing.len() {
                if inputs[missing[begin]].parts.is_empty() {
                    publish(missing[begin], Vec::new())?;
                    begin += 1;
                    continue;
                }
                let mut end = begin;
                let mut channels = 0;
                let mut bytes = 0;
                while end < missing.len() {
                    let input = &inputs[missing[end]];
                    if end > begin
                        && (input.parts.is_empty()
                            || end - begin >= MAX_REQUESTS
                            || channels + input.parts.len() > MAX_CHANNELS
                            || bytes + input.program.len() > MAX_PROGRAM_BYTES)
                    {
                        break;
                    }
                    channels += input.parts.len();
                    bytes += input.program.len();
                    end += 1;
                }
                let (group, program, names) = loop {
                    let group = missing[begin..end]
                        .iter()
                        .map(|index| &inputs[*index])
                        .collect::<Vec<_>>();
                    let (program, names) = if group.len() == 1 {
                        (
                            group[0].program.clone(),
                            group[0]
                                .parts
                                .iter()
                                .map(|part| part.name.clone())
                                .collect(),
                        )
                    } else {
                        combined_program(&group)?
                    };
                    if group.len() == 1 || program.len() <= MAX_PROGRAM_BYTES {
                        break (group, program, names);
                    }
                    end = begin + (end - begin) / 2;
                };
                let names = names.iter().map(String::as_str).collect::<Vec<_>>();
                let raw = run_forcer_uncached(settings, &program, &names)?;
                jobs += 1;
                if raw.len() != names.len() {
                    return Err(VakintError::InvalidGenericExpression(
                        "FORCER group returned the wrong number of scalar channels".into(),
                    ));
                }
                let mut offset = 0;
                // Peers can use this group's results while later groups execute.
                for (&index, input) in missing[begin..end].iter().zip(group) {
                    publish(index, raw[offset..offset + input.parts.len()].to_vec())?;
                    offset += input.parts.len();
                }
                begin = end;
            }
            Ok(())
        };
    let outputs = if let Some(keys) = keys {
        let _stage = StageSpan::new("p_integral.forcer.cache_lookup_wait");
        FORCER_CACHE.run_batch(&keys, &mut compute)?
    } else {
        let _stage = StageSpan::new("p_integral.forcer.cache.bypass");
        let mut outputs = vec![Vec::new(); inputs.len()];
        compute(
            &(0..inputs.len()).collect::<Vec<_>>(),
            &mut |index, output| {
                outputs[index] = output;
                Ok(())
            },
        )?;
        outputs
    };
    let run_elapsed = run_started.elapsed();
    let parsing_started = Instant::now();
    let epsilon = vk_symbol!(&settings.epsilon_symbol);
    let values = inputs
        .iter()
        .zip(outputs)
        .map(|(input, output)| combine_parts(&input.parts, &output, epsilon))
        .collect::<Result<Vec<_>, _>>()?;
    tracing::info!(
        uv = true,
        rstar = true,
        timing = true,
        requests = inputs.len(),
        computed_requests = misses,
        reused_requests = inputs.len() - misses,
        form_jobs = jobs,
        run_or_cache_seconds = run_elapsed.as_secs_f64(),
        output_parsing_seconds = parsing_started.elapsed().as_secs_f64(),
        seconds = started.elapsed().as_secs_f64(),
        "Reduced p-integral request batch with FORCER"
    );
    Ok(values)
}
