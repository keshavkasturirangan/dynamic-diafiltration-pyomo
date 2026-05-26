# DATA2 Workflow for DATA3

This folder holds a narrow compatibility workflow for checking whether the
legacy `utility.py` pipeline can be pointed at the NF270 / DATA3 Excel sheets
without changing the old solver logic.

What this workflow does:

1. Loads NF270 `.xlsx` sheets through the refactored Excel loader.
2. Uses the refactored conductivity-to-concentration conversion upstream.
3. Hands the converted `data_stru` to the legacy `utility.py` solver.
4. Runs the old plots so the mass-vs-time and concentration-vs-time figures
   can be compared with the refactored DATA3 outputs.

Why this exists:

- It gives us a direct apples-to-apples test of whether the older DATA2-style
  estimation path can reproduce the single-salt DATA3 campaign once the data
  are loaded and normalized the same way.
- It also makes any parameter drift easier to point at: if `Lp` stays stable
  but `B` or `sigma` move, the source is more likely the concentration
  conversion or the B-form choice rather than the mass balance itself.

How to run:

```bash
python refactored_codes_v1/DATA2_workflow_for_DATA3/run_data2_workflow_for_data3.py
```

Outputs are written inside this directory so they stay separate from the
main refactored campaign artifacts. If a reference NF270 summary is found in
`UnifiedFramework/DATA3/results/paper_artifacts/nf270/`, the runner also
writes a side-by-side comparison CSV/JSON for the overlapping sheet IDs.
