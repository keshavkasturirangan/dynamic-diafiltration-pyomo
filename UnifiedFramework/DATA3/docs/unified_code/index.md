# Unified Code Documentation (v24)

> **Canonical code path (source of truth):** `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/`  
> Runner entrypoint: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py`  
> Core pipeline: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`  
> Legacy root-level scripts (for example `utility.py`, `run_*.py`) are compatibility/reference paths, not primary development paths.


This section documents the current unified experimental loader and modeling code under:

- `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`
- `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py`
- `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/conductivity_paper.py`

The goal is to provide:

- onboarding context for new users,
- implementation context for refactoring,
- packaging context for future ReadTheDocs docs.

## Module map

- `unified_codebase_library.py`: core data model, file ingestion, validation, conductivity->concentration utilities, dynamic Pyomo model construction, ParmEst integration, DoE/FIM helpers.
- `unified_codebase_runfile.py`: user-facing entry script that configures one run.
- `conductivity_paper.py`: conductivity physics helpers (Shedlovsky variant + MSA transport routines).

## Recommended read order

1. `unified_codebase_runfile.md` (how users run the workflow)
2. `workflow_algorithms.md` (step-by-step runtime algorithms)
3. `unified_codebase_library.md` (core system design and APIs)
4. `conductivity_paper.md` (measurement-model internals)

## Future RTD migration notes

This folder is written in Markdown so it can be migrated to RTD/MyST quickly.

For Sphinx/MyST, add a `toctree` similar to:

```md
```{toctree}
:maxdepth: 2

unified_codebase_runfile
unified_codebase_library
conductivity_paper
```
```

