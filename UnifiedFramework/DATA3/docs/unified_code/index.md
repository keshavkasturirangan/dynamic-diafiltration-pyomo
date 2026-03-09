# Unified Code Documentation (v24)

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
2. `unified_codebase_library.md` (core system design and APIs)
3. `conductivity_paper.md` (measurement-model internals)

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

