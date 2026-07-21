//! End-of-run summaries shared by the alignment and collapse commands.

use std::fs::File;
use std::io::{BufWriter, Result, Write};
use std::path::Path;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SummaryTable {
    name: String,
    headers: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl SummaryTable {
    pub fn new(name: impl Into<String>, headers: &[&str]) -> Self {
        assert!(
            headers.len() >= 2,
            "Summary tables require at least two columns"
        );
        Self {
            name: name.into(),
            headers: headers.iter().map(|header| header.to_string()).collect(),
            rows: Vec::new(),
        }
    }

    pub fn push_row<I, S>(&mut self, values: I)
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        let row = values.into_iter().map(Into::into).collect::<Vec<_>>();
        assert_eq!(
            row.len(),
            self.headers.len(),
            "Summary row does not match its table header"
        );
        self.rows.push(row);
    }

    fn render(&self) -> String {
        let mut widths = self.headers.iter().map(String::len).collect::<Vec<_>>();
        for row in &self.rows {
            for (index, value) in row.iter().enumerate() {
                widths[index] = widths[index].max(value.len());
            }
        }

        let separator = format!(
            "+{}+",
            widths
                .iter()
                .map(|width| "-".repeat(width + 2))
                .collect::<Vec<_>>()
                .join("+")
        );
        let render_row = |row: &[String]| {
            format!(
                "| {} |",
                row.iter()
                    .enumerate()
                    .map(|(index, value)| format!("{:<width$}", value, width = widths[index]))
                    .collect::<Vec<_>>()
                    .join(" | ")
            )
        };

        let mut lines = vec![
            self.name.clone(),
            separator.clone(),
            render_row(&self.headers),
            separator.clone(),
        ];
        lines.extend(self.rows.iter().map(|row| render_row(row)));
        lines.push(separator);
        lines.join("\n")
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RunSummary {
    title: String,
    tables: Vec<SummaryTable>,
}

impl RunSummary {
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            tables: Vec::new(),
        }
    }

    pub fn add_table(&mut self, table: SummaryTable) {
        self.tables.push(table);
    }

    pub fn render(&self) -> String {
        let tables = self
            .tables
            .iter()
            .map(SummaryTable::render)
            .collect::<Vec<_>>()
            .join("\n\n");
        format!("{}\n{}", self.title, tables)
    }

    /// Print the summary to stderr and optionally write the same values in tidy TSV form.
    pub fn emit(&self, output: Option<&Path>) -> Result<()> {
        eprintln!("\n{}", self.render());
        if let Some(path) = output {
            self.write_tsv(path)?;
            info!("Wrote run summary to {}", path.display());
        }
        Ok(())
    }

    pub fn write_tsv(&self, path: &Path) -> Result<()> {
        let mut writer = BufWriter::new(File::create(path)?);
        writeln!(writer, "section\tscope\tmetric\tvalue")?;
        for table in &self.tables {
            for row in &table.rows {
                let scope = sanitize_tsv(&row[0]);
                for (header, value) in table.headers.iter().skip(1).zip(row.iter().skip(1)) {
                    if !value.is_empty() {
                        writeln!(
                            writer,
                            "{}\t{}\t{}\t{}",
                            sanitize_tsv(&table.name),
                            scope,
                            sanitize_tsv(header),
                            sanitize_tsv(value)
                        )?;
                    }
                }
            }
        }
        writer.flush()
    }
}

fn sanitize_tsv(value: &str) -> String {
    value.replace(['\t', '\n', '\r'], " ")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn renders_aligned_table_and_writes_tidy_tsv() {
        let mut table = SummaryTable::new("Results", &["Scope", "Input", "Output"]);
        table.push_row(["All", "10", "8"]);
        table.push_row(["reference_a", "6", "5"]);
        let mut summary = RunSummary::new("Run summary");
        summary.add_table(table);

        let rendered = summary.render();
        assert!(rendered.contains("| Scope       | Input | Output |"));
        assert!(rendered.contains("| reference_a | 6     | 5      |"));

        let directory = tempfile::tempdir().unwrap();
        let output = directory.path().join("summary.tsv");
        summary.write_tsv(&output).unwrap();
        let contents = std::fs::read_to_string(output).unwrap();
        assert_eq!(
            contents,
            "section\tscope\tmetric\tvalue\nResults\tAll\tInput\t10\nResults\tAll\tOutput\t8\nResults\treference_a\tInput\t6\nResults\treference_a\tOutput\t5\n"
        );
    }

    #[test]
    fn sanitizes_tabs_and_newlines_in_tsv_values() {
        let mut table = SummaryTable::new("Results", &["Scope", "Value"]);
        table.push_row(["ref\tone", "line\none"]);
        let mut summary = RunSummary::new("Run summary");
        summary.add_table(table);

        let directory = tempfile::tempdir().unwrap();
        let output = directory.path().join("summary.tsv");
        summary.write_tsv(&output).unwrap();
        assert!(std::fs::read_to_string(output)
            .unwrap()
            .contains("ref one\tValue\tline one"));
    }
}
