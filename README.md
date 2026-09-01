# Meta-analysis through random effects

MATLAB reproduction code for Jan R. Magnus and Andrey L. Vasnev, *More information, less precision: Meta-analysis through random effects*.

## Repository contents

- `Simple_figure_v4_github.m` produces manuscript Figures 1 and 2.
- `Meta_analysis_v10_github.m` produces manuscript Figures 3–6 and prints manuscript Table 1. Unrounded table values remain available as `table1_results`.
- `Menkveld_table2_v2_github.m` reproduces all 24 rows of manuscript Table 2. Unrounded values remain available as `table2_results`.
- `RT_research_results.csv` contains the Menkveld et al. (2024) data used by the case-12 figure and Table 2.

## Requirements

- MATLAB with the Statistics and Machine Learning Toolbox
- MATLAB Optimization Toolbox

The current programs were checked and run successfully in MATLAB R2025b.

## Data location

The two programs that use the Menkveld data expect this file in the
repository root:

```text
RT_research_results.csv
```

## Running the programs

Set the MATLAB Current Folder to the repository root, then run:

```matlab
Simple_figure_v4_github
Meta_analysis_v10_github
Menkveld_table2_v2_github
```

By default, the figure programs create figures without writing files. To export vector PDFs:

- Select the desired entries in `export_figures` in `Simple_figure_v4_github.m`.
- Select the desired entries in `print_figures` in `Meta_analysis_v10_github.m`.

Exported files are written to `output/pdf`.

Both table programs perform numerical checks before printing their results in the MATLAB Command Window.

## Table 2 reproduction notes

- Quantiles use the Excel `PERCENTILE.INC` definition used for the submitted table.
- Quantiles and IQRs use all 164 observations for each hypothesis and stage.
- The meta-analysis removes the highest and lowest 5% of estimates and reported standard errors.
- The submitted Stage-1 results also remove the lowest 5% of peer ratings.
- The archived H6 Stage-1 result is preserved at the original optimizer's four-iteration stopping point so that it remains stable across MATLAB releases.
