# Repository Target Layout (Non-Breaking Reorg)

Updated: 2026-04-01

## Goals

- Keep `UnifiedFramework/DATA3` as the canonical execution root during transition.
- Remove root-level clutter by moving legacy notebooks/scripts/assets into explicit domains.
- Separate source code, validation assets, and generated artifacts.
- Preserve nightly validation stability while path migration is in progress.

## Target top-level structure

```text
.
├── UnifiedFramework/
│   └── DATA3/                          # canonical runnable framework (keep)
├── legacy/
│   ├── data1_matlab/                   # migrated DATA1_matlab assets + scripts
│   ├── notebooks/
│   │   ├── data1/
│   │   └── data2/
│   └── scripts/                        # root compatibility scripts migrated here
├── docs/
│   ├── validation/                     # high-level validation docs (keep)
│   └── repo_reorg/                     # reorg plans + reports
├── tools/
│   └── repo/                           # repository-management tooling
├── data/
│   ├── published/                      # paper PDFs / archives (optional phase)
│   └── extracted/                      # extracted paper assets (optional phase)
└── artifacts/
    ├── figures/                        # shareable generated figures (optional phase)
    └── reports/                        # shareable generated reports (optional phase)
```

## Non-breaking policy

- Do not relocate `UnifiedFramework/DATA3` in Phase 1–2.
- For every moved script/notebook in early phases, keep a compatibility shim at old path.
- Update tests and nightly scripts in slices; run nightly gate after each slice.

## Execution slices

1. Phase A: Planning + path audit + move map (this stage).
2. Phase B: Move root-level notebooks/scripts into `legacy/` with compatibility wrappers.
3. Phase C: Consolidate DATA1 legacy assets under `legacy/data1_matlab`.
4. Phase D: Optional migration of large generated/reproduction outputs to `artifacts/` + symlinks/wrappers.
5. Phase E: Remove deprecated paths once references hit zero.
