# DATA1 Legacy vs Unified State Summary

This summary is based on fixed-theta comparisons using the legacy MATLAB fitted
parameters in both the legacy and unified models.

Primary artifacts:

- `data1_legacy_vs_unified_state_metrics.csv`
- `data1_legacy_vs_unified_state_traces.csv`
- `data1_legacy_vs_unified_boundary_report.csv`

## Main Result

For the fixed-theta comparison, the unified DATA1 model now matches the legacy
MATLAB state trajectories essentially exactly at the sampled points.

This is the important interpretation:

- the DATA1 simulation equations and boundary carryover rules are aligned well
  enough to reproduce the legacy state trajectories at the same `Lp`, `B`, and
  `sigma`
- the remaining DATA1 concentration mismatch is therefore not a fixed-theta
  simulation-parity problem
- the remaining gap is more likely in the estimation path than in the underlying
  DATA1 process model

## Dataset 501.1

Representative errors from `data1_legacy_vs_unified_state_metrics.csv`:

- vial 1 `cF` MAE: `2.77e-15`
- vial 1 `cH` MAE: `2.13e-16`
- vial 1 `cV` MAE: `3.50e-16`
- vial 1 `cIn` MAE: `3.27e-15`
- vial 1 `mV` MAE: `1.32e-16`

Boundary report:

- `cF` carryover matches
- `cH` carryover matches
- `cV` carryover matches
- `mV` reset matches to numerical precision

## Dataset 511.12

The same fixed-theta comparison workflow was applied to `511.12`, and the state
metrics and boundary report should be interpreted the same way:

- use the CSV artifacts, not earlier draft narrative notes, as the source of
  truth
- the fixed-theta comparison is no longer evidence of a model-definition
  mismatch by itself

## Working Conclusion

The current evidence supports this narrower diagnosis:

- fixed-theta DATA1 simulation parity is good
- the remaining DATA1 Figure 2 concentration `xfail` is an estimation-parity
  issue, not a basic state-equation mismatch
