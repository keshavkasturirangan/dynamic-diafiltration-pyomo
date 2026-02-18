# Validation Docs

This folder tracks DATA1/DATA2 published-result reproduction for unified-code validation.

## Files

- `paper_target_matrix.md`: required figure/table targets and stage plan.
- `paper_reference_values_template.csv`: side-by-side paper/unified reference template.

## Reproduction scaffold

Run canonical datasets:

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/scripts/reproduce/reproduce_data1_data2.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --solver ipopt \
  --nfe 30
```

Outputs are written under:

- `results/reproduction/<run_id>/summaries/*.json`
- `results/reproduction/<run_id>/tables/unified_metrics.csv`
- `results/reproduction/<run_id>/tables/paper_vs_unified_side_by_side.csv`
