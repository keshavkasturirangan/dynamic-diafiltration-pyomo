# Visual Validation Share Index

Updated: 2026-04-01

## Core artifacts

- Paper PNG extract root:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/reproduction/20260306-023356-paper-pdf-extract/pdf_extract`
- Unified-vs-paper composite PNGs:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/simulation_validation/composites`
- Validation report (CSV):
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/simulation_validation/simulation_validation_figure_report.csv`
- Validation report (Markdown):
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/simulation_validation/simulation_validation_figure_report.md`

## Current status snapshot

- Figure targets evaluated: 23
- PASS: 23
- FAIL: 0
- Comparator: composite-image similarity vs extracted paper panel images (per target threshold).

## Table side-by-side references

- DATA1 tables (latest full run snapshot used in validation docs):
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/reproduction/20260221-011609/tables/paper_vs_unified_side_by_side.csv`
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/reproduction/20260221-011609/tables/data1_table1_side_by_side.csv`
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/reproduction/20260221-011609/tables/data1_table2_side_by_side.csv`
- DATA2 main paper baseline tables:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/reproduction/20260306-mainpaper-table-baseline/tables/data2_table3_side_by_side.csv`
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/reproduction/20260306-mainpaper-table-baseline/tables/data2_table4_side_by_side.csv`
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/reproduction/20260306-mainpaper-table-baseline/tables/data2_table5_side_by_side.csv`
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/reproduction/20260306-mainpaper-table-baseline/tables/data2_table6_side_by_side.csv`

## Suggested share package

Share these three locations together:

1. `.../pdf_extract/` (published-paper evidence)
2. `.../results/simulation_validation/composites/` (unified-generated comparison visuals)
3. `.../results/simulation_validation/simulation_validation_figure_report.csv` (per-target PASS/FAIL)
