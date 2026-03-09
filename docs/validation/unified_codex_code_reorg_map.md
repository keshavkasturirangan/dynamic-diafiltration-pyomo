# `unified_codex_code` Reorg Map

Date: 2026-03-09
Scope: path migration reference for `/Users/kkasturi/GitHub/unified_codex_code`

## Purpose

This note captures the large subtree reorganization observed in `unified_codex_code`.
It is intended as a migration checklist and safety reference only.

Important:
- This document does **not** authorize deleting files from `codex/unified_code`.
- The current branch should treat the reorg as a path-mapping exercise until all references are intentionally migrated.
- Validation targets remain anchored to figures/tables that appear in `published_works/`.

## High-Level Mapping

| Old path root | New path root |
|---|---|
| `docs/validation/` | `UnifiedFramework/DATA3/docs/validation/` |
| `docs/unified_code/` | `UnifiedFramework/DATA3/docs/unified_code/` |
| `figures/` | `UnifiedFramework/DATA3/figures/` |
| `results/` | `UnifiedFramework/DATA3/results/` |
| `published_works/` | `UnifiedFramework/DATA3/published_works/` |
| `overleaf-latex/` | `UnifiedFramework/DATA3/overleaf-latex/` |
| `scripts/` | `UnifiedFramework/DATA3/scripts/` |
| `UnifiedFramework/DATA3/ExperimentalDataLoader/` | `UnifiedFramework/DATA3/ExperimentalDataLoader_ParamEst_n_DoEanalysis/` |

## Representative File Mapping

| Old path | New path |
|---|---|
| `docs/validation/README.md` | `UnifiedFramework/DATA3/docs/validation/README.md` |
| `docs/validation/developmentplan.md` | `UnifiedFramework/DATA3/docs/validation/developmentplan.md` |
| `docs/unified_code/unified_codebase_library.md` | `UnifiedFramework/DATA3/docs/unified_code/unified_codebase_library.md` |
| `figures/Bpervial.png` | `UnifiedFramework/DATA3/figures/Bpervial.png` |
| `results/reproduction/20260221-004009/tables/data1_table1_side_by_side.csv` | `UnifiedFramework/DATA3/results/reproduction/20260221-004009/tables/data1_table1_side_by_side.csv` |
| `scripts/reproduce/reproduce_data1_data2.py` | `UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py` |
| `published_works/DATA1_main.pdf` | `UnifiedFramework/DATA3/published_works/DATA1_main.pdf` |
| `UnifiedFramework/DATA3/ExperimentalDataLoader/UnifiedCode/unified_codebase_library.py` | `UnifiedFramework/DATA3/ExperimentalDataLoader_ParamEst_n_DoEanalysis/UnifiedCode/unified_codebase_library.py` |
| `UnifiedFramework/DATA3/ExperimentalDataLoader/UnifiedCode/experiment_load_v24.py` | `UnifiedFramework/DATA3/ExperimentalDataLoader_ParamEst_n_DoEanalysis/UnifiedCode/experiment_load_v24.py` |

## Migration Guidance

1. Update documentation paths first.
2. Update workflow and script references second.
3. Update imports and internal code references only after the target location is confirmed.
4. Do not treat old-path deletions as safe-to-commit removals unless the new-path file is verified and the branch intentionally adopts the new layout.

## Safety Rule For `codex/unified_code`

- Documentation commits may reference the `unified_codex_code` reorg.
- They must not delete current files from this repository solely because the external repo moved them.
- Any future structural migration should be done as an explicit, separately reviewed change set.
