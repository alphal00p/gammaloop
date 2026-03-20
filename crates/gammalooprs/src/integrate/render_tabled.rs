use colored::Colorize;
use itertools::Itertools;
use tabled::{
    Table,
    builder::Builder,
    settings::{
        Alignment, Modify, Panel, Span,
        object::{Cell, Rows},
        style::Style,
        themes::BorderCorrection,
        width::Width,
    },
};

use crate::utils::normalize_tabled_separator_rows;

use super::{
    display::{StyledText, TextColor},
    status_update::{
        DiscreteMaxWeightDetailsSection, MainResultsRowGroupKind, MainResultsSection,
        MaxWeightDetailsSection, StatusUpdate,
    },
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct TabledRenderOptions {
    pub max_table_width: usize,
    pub show_statistics: bool,
    pub show_max_weight_details: bool,
    pub show_top_discrete_grid: bool,
    pub show_discrete_contributions_sum: bool,
    pub show_max_weight_info_for_discrete_bins: bool,
}

impl Default for TabledRenderOptions {
    fn default() -> Self {
        Self {
            max_table_width: 250,
            show_statistics: true,
            show_max_weight_details: true,
            show_top_discrete_grid: false,
            show_discrete_contributions_sum: false,
            show_max_weight_info_for_discrete_bins: false,
        }
    }
}

struct StatusTable {
    table: Table,
    separator_after_rows: Vec<usize>,
    hidden_vertical_boundaries: Vec<usize>,
    full_row_vertical_count: usize,
    suppress_header_middle_separator: bool,
    suppress_header_tail_separator: bool,
}

fn render_styled_text(text: &StyledText) -> String {
    text.spans
        .iter()
        .map(|span| {
            let mut styled = match span.style.color {
                Some(TextColor::Green) => span.text.green(),
                Some(TextColor::Blue) => span.text.blue(),
                Some(TextColor::Red) => span.text.red(),
                None => span.text.normal(),
            };
            if span.style.dimmed {
                styled = styled.dimmed();
            }
            if span.style.bold {
                styled = styled.bold();
            }
            styled.to_string()
        })
        .collect()
}

fn maybe_render<T>(value: &Option<super::display::DisplayField<T>>) -> String {
    value
        .as_ref()
        .map(|value| render_styled_text(&value.display))
        .unwrap_or_default()
}

fn suppress_iteration_header_separators(
    rendered: &str,
    suppress_middle_separator: bool,
    suppress_tail_separator: bool,
) -> String {
    rendered
        .lines()
        .enumerate()
        .map(|(line_index, line)| {
            if line_index != 1 {
                return line.to_string();
            }

            let vertical_positions = line.match_indices('│').map(|(idx, _)| idx).collect_vec();
            if vertical_positions.len() < 3 {
                return line.to_string();
            }

            let mut updated = line.to_string();
            if suppress_tail_separator {
                let suppressed_index = vertical_positions[vertical_positions.len() - 2];
                updated.replace_range(suppressed_index..suppressed_index + '│'.len_utf8(), " ");
            }
            if suppress_middle_separator && vertical_positions.len() >= 4 {
                let suppressed_index = vertical_positions[1];
                updated.replace_range(suppressed_index..suppressed_index + '│'.len_utf8(), " ");
            }
            updated
        })
        .join("\n")
}

fn hide_hidden_vertical_boundaries(
    rendered: &str,
    hidden_vertical_boundaries: &[usize],
    full_row_vertical_count: usize,
) -> String {
    rendered
        .lines()
        .map(|line| {
            let vertical_positions = line.match_indices('│').map(|(idx, _)| idx).collect_vec();
            if vertical_positions.len() != full_row_vertical_count {
                return line.to_string();
            }

            let mut updated = line.to_string();
            for boundary in hidden_vertical_boundaries.iter().rev() {
                let Some(position) = vertical_positions.get(boundary + 1).copied() else {
                    continue;
                };
                updated.replace_range(position..position + '│'.len_utf8(), " ");
            }
            updated
        })
        .join("\n")
}

fn insert_separator_rows(rendered: &str, separator_after_rows: &[usize]) -> String {
    if separator_after_rows.is_empty() {
        return rendered.to_string();
    }

    let mut lines = rendered.lines().map(str::to_string).collect_vec();
    let width = lines
        .first()
        .map(|line| line.chars().count())
        .unwrap_or_default();
    if width < 2 {
        return rendered.to_string();
    }
    let separator = format!("├{}┤", "─".repeat(width - 2));

    for (inserted, row) in separator_after_rows.iter().enumerate() {
        let insert_at = *row + 2 + inserted;
        lines.insert(insert_at, separator.clone());
    }

    lines.join("\n")
}

fn suppress_spanned_metadata_header_separators(rendered: &str) -> String {
    rendered
        .lines()
        .map(|line| {
            if !line.contains("χ²/dof") || !line.contains("mwi") {
                return line.to_string();
            }

            let mut updated = line.to_string();
            let separators_to_remove = [("χ²/dof", "mwi"), ("Δ [σ]", "Δ [%]")];

            for (left_label, right_label) in separators_to_remove {
                let Some(left_start) = updated.find(left_label) else {
                    continue;
                };
                let Some(right_start) = updated.find(right_label) else {
                    continue;
                };
                let left_end = left_start + left_label.len();
                if let Some((separator_index, _)) = updated[left_end..right_start]
                    .match_indices('│')
                    .next_back()
                {
                    let separator_index = left_end + separator_index;
                    updated.replace_range(separator_index..separator_index + '│'.len_utf8(), " ");
                }
            }

            updated
        })
        .join("\n")
}

fn render_tables_with_shared_width(mut tables: Vec<StatusTable>, max_table_width: usize) -> String {
    let max_width = tables
        .iter()
        .map(|table| table.table.total_width())
        .max()
        .unwrap_or(0)
        .min(max_table_width.max(1));

    for table in &mut tables {
        if table.table.total_width() < max_width {
            table.table.with(Width::increase(max_width));
        }
    }

    tables
        .into_iter()
        .map(|table| {
            let mut rendered = table.table.to_string();
            rendered = hide_hidden_vertical_boundaries(
                &rendered,
                &table.hidden_vertical_boundaries,
                table.full_row_vertical_count,
            );
            if table.suppress_header_middle_separator || table.suppress_header_tail_separator {
                rendered = suppress_iteration_header_separators(
                    &rendered,
                    table.suppress_header_middle_separator,
                    table.suppress_header_tail_separator,
                );
            }
            rendered = suppress_spanned_metadata_header_separators(&rendered);
            rendered = insert_separator_rows(&rendered, &table.separator_after_rows);
            normalize_tabled_separator_rows(&rendered)
        })
        .collect_vec()
        .join("\n")
}

fn build_main_results_table(
    section: &MainResultsSection,
    options: &TabledRenderOptions,
) -> StatusTable {
    let visible_row_groups = section
        .row_groups
        .iter()
        .filter(|group| match group.kind {
            MainResultsRowGroupKind::All => true,
            MainResultsRowGroupKind::Sum => options.show_discrete_contributions_sum,
            MainResultsRowGroupKind::Bins => options.show_top_discrete_grid,
        })
        .collect_vec();
    let show_discrete_columns = section.has_discrete_columns
        && (options.show_top_discrete_grid || options.show_discrete_contributions_sum);
    let slot_block_width = if show_discrete_columns { 5 } else { 2 };
    let metadata_columns = if section.has_target_columns { 4 } else { 2 };
    let n_columns = 2 + section.slot_headers.len() * slot_block_width + metadata_columns;

    let mut builder = Builder::new();
    let mut first_row = vec![
        render_styled_text(&section.header_left),
        String::new(),
        render_styled_text(&section.header_middle),
    ];
    first_row.resize(n_columns, String::new());
    first_row[n_columns - 2] = render_styled_text(&section.header_tail);
    builder.push_record(first_row);

    let mut header_row = vec![
        render_styled_text(&section.contribution_header),
        String::new(),
    ];
    for slot_header in &section.slot_headers {
        header_row.push(render_styled_text(slot_header));
        header_row.extend(std::iter::repeat_n(String::new(), slot_block_width - 1));
    }
    header_row.push("χ²/dof".blue().bold().to_string());
    header_row.push("mwi".blue().bold().to_string());
    if section.has_target_columns {
        header_row.push("Δ [σ]".blue().bold().to_string());
        header_row.push("Δ [%]".blue().bold().to_string());
    }
    builder.push_record(header_row);

    for group in &visible_row_groups {
        for row in &group.rows {
            let mut record = vec![
                render_styled_text(&row.contribution.display),
                render_styled_text(&row.component.display),
            ];
            for slot_cells in &row.slot_cells {
                record.push(maybe_render(&slot_cells.value));
                record.push(maybe_render(&slot_cells.relative_error));
                if show_discrete_columns {
                    record.push(maybe_render(&slot_cells.sample_fraction));
                    record.push(maybe_render(&slot_cells.sample_count));
                    record.push(maybe_render(&slot_cells.target_pdf));
                }
            }
            record.push(maybe_render(&row.chi_sq));
            record.push(maybe_render(&row.max_weight_impact));
            if section.has_target_columns {
                record.push(maybe_render(&row.delta_sigma));
                record.push(maybe_render(&row.delta_percent));
            }
            builder.push_record(record);
        }
    }

    let mut table = builder.build();
    table.modify((0, 0), Span::column(2));
    table.modify((0, 2), Span::column((n_columns - 4) as isize));
    table.modify((0, n_columns - 2), Span::column(2));
    table.modify((1, 0), Span::column(2));
    for (slot_index, _) in section.slot_headers.iter().enumerate() {
        table.modify(
            (1, 2 + slot_index * slot_block_width),
            Span::column(slot_block_width as isize),
        );
    }

    let first_metadata_column = 2 + section.slot_headers.len() * slot_block_width;

    let mut separator_rows = vec![1usize, 2usize];
    let mut row_offset = 2usize;
    for (group_index, group) in visible_row_groups.iter().enumerate() {
        row_offset += group.rows.len();
        if group_index + 1 < visible_row_groups.len() {
            separator_rows.push(row_offset);
        }
    }
    separator_rows.sort_unstable();
    separator_rows.dedup();

    table.with(Style::rounded().remove_horizontals());
    table.with(BorderCorrection::span());
    table.with(Modify::new(Rows::new(0..)).with(Alignment::left()));
    table.with(Modify::new(Cell::new(0, 2)).with(Alignment::center()));
    table.with(Modify::new(Cell::new(0, n_columns - 2)).with(Alignment::center()));
    table.with(Modify::new(Rows::new(1..2)).with(Alignment::center()));

    let mut hidden_vertical_boundaries = vec![0usize];
    for slot_index in 0..section.slot_headers.len() {
        let block_start = 2 + slot_index * slot_block_width;
        hidden_vertical_boundaries.extend(block_start..(block_start + slot_block_width - 1));
    }
    hidden_vertical_boundaries.push(first_metadata_column);
    if section.has_target_columns {
        hidden_vertical_boundaries.push(first_metadata_column + 2);
    }

    StatusTable {
        table,
        separator_after_rows: separator_rows
            .into_iter()
            .map(|row| row.saturating_sub(1))
            .collect(),
        hidden_vertical_boundaries,
        full_row_vertical_count: n_columns + 1,
        suppress_header_middle_separator: true,
        suppress_header_tail_separator: true,
    }
}

fn build_max_weight_details_table(section: &MaxWeightDetailsSection) -> StatusTable {
    let mut builder = Builder::new();
    builder.push_record([
        "Integrand".blue().bold().to_string(),
        String::new(),
        "Max eval".blue().bold().to_string(),
        "Max eval coordinates".blue().bold().to_string(),
    ]);

    for group in &section.rows_by_slot {
        for row in group {
            builder.push_record([
                render_styled_text(&row.slot.display),
                render_styled_text(&row.component_sign.display),
                render_styled_text(&row.max_eval.display),
                render_styled_text(&row.coordinates.display),
            ]);
        }
    }

    let mut table = builder.build();
    table.modify((0, 0), Span::column(2));
    table.with(Panel::header(
        "Maximum weight details".bold().green().to_string(),
    ));
    table.with(Style::rounded().remove_horizontals());
    table.with(BorderCorrection::span());
    table.with(Modify::new(Rows::new(0..2)).with(Alignment::center()));
    table.with(Modify::new(Rows::new(2..)).with(Alignment::left()));

    let mut separator_rows = vec![0usize, 1usize];
    let mut row_offset = 1usize;
    for (group_index, group) in section.rows_by_slot.iter().enumerate() {
        row_offset += group.len();
        if group_index + 1 < section.rows_by_slot.len() {
            separator_rows.push(row_offset);
        }
    }

    StatusTable {
        table,
        separator_after_rows: separator_rows,
        hidden_vertical_boundaries: vec![0],
        full_row_vertical_count: 5,
        suppress_header_middle_separator: false,
        suppress_header_tail_separator: false,
    }
}

fn build_discrete_max_weight_details_table(
    section: &DiscreteMaxWeightDetailsSection,
) -> StatusTable {
    let mut builder = Builder::new();
    let mut header = vec![
        render_styled_text(&section.contribution_header),
        String::new(),
    ];
    for slot_header in &section.slot_headers {
        header.push(render_styled_text(slot_header));
    }
    header.push("Max eval coordinates".blue().bold().to_string());
    builder.push_record(header);

    for group in &section.row_groups {
        for row in group {
            let mut record = vec![
                render_styled_text(&row.contribution.display),
                render_styled_text(&row.component_sign.display),
            ];
            record.extend(row.slot_values.iter().map(|value| {
                value
                    .as_ref()
                    .map(|value| render_styled_text(&value.display))
                    .unwrap_or_default()
            }));
            record.push(
                row.slot_coordinates
                    .iter()
                    .map(|entry| {
                        format!(
                            "{}: {}",
                            render_styled_text(&entry.slot.display),
                            render_styled_text(&entry.coordinates.display),
                        )
                    })
                    .join("\n"),
            );
            builder.push_record(record);
        }
    }

    let mut table = builder.build();
    table.modify((0, 0), Span::column(2));
    table.with(Panel::header(
        "Maximum weight details by discrete bin"
            .bold()
            .green()
            .to_string(),
    ));

    let mut separator_rows = vec![1usize, 2usize];
    let mut row_offset = 1usize;
    for (group_index, group) in section.row_groups.iter().enumerate() {
        row_offset += group.len();
        if group_index + 1 < section.row_groups.len() {
            separator_rows.push(row_offset + 1);
        }
    }

    table.with(Style::rounded().remove_horizontals());
    table.with(BorderCorrection::span());
    table.with(Modify::new(Rows::new(0..2)).with(Alignment::center()));
    table.with(Modify::new(Rows::new(2..)).with(Alignment::left()));

    StatusTable {
        table,
        separator_after_rows: separator_rows
            .into_iter()
            .map(|row| row.saturating_sub(1))
            .collect(),
        hidden_vertical_boundaries: vec![0],
        full_row_vertical_count: section.slot_headers.len() + 4,
        suppress_header_middle_separator: false,
        suppress_header_tail_separator: false,
    }
}

pub(crate) fn render_status_update(update: &StatusUpdate, options: &TabledRenderOptions) -> String {
    let mut tables = vec![build_main_results_table(&update.main_results, options)];
    if options.show_max_weight_details
        && let Some(section) = update.max_weight_details.as_ref()
    {
        tables.push(build_max_weight_details_table(section));
    }
    if options.show_max_weight_details
        && options.show_max_weight_info_for_discrete_bins
        && let Some(section) = update.discrete_max_weight_details.as_ref()
    {
        tables.push(build_discrete_max_weight_details_table(section));
    }
    if options.show_statistics
        && let Some(section) = update.statistics.as_ref()
    {
        tables.push(StatusTable {
            table: section.raw.build_status_table(),
            separator_after_rows: Vec::new(),
            hidden_vertical_boundaries: Vec::new(),
            full_row_vertical_count: 0,
            suppress_header_middle_separator: false,
            suppress_header_tail_separator: false,
        });
    }

    render_tables_with_shared_width(tables, options.max_table_width)
}
