# Refactor Prompt for Unified Codebase (Option B)

I would like your help significantly refactoring this code base. This code base grew organically over several publications and involves several contributors. We want to unify the existing code into an easier-to-maintain foundation that can then be extended. We will validate the refactor by reproducing results (plots, tables) from the DATA1 and DATA2 papers.

Immediate challenge to be addressed: validation of DATA1 and DATA2 plots and tables.

Refactor the unified analysis workflow centered on `unified_codebase_library.py` and `unified_codebase_runfile.py` into a clean, extensible, object-oriented pipeline without changing scientific intent.

Do not start coding immediately. First inspect the repository and produce an implementation plan. Then execute the work step by step.

## Primary Goal

Create a unified, maintainable workflow that supports DATA1, DATA2, and DATA3, preserves existing working behavior where possible, validates DATA1/DATA2 plots and tables, and adds XLSX loading support for DATA3.

## Reference Files and Code to Inspect

- `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`
- `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py`
- `utility.py`
- `unified_codebase_v24`
- `experiment_load_v24.py`
- `conductivity_paper.py`
- Existing tests
- Existing nightly pytest configuration and helpers
- Existing parmest-related code paths

## Revised Refactor Organization

Use the following high-level organization instead of a long chain of thin top-level pipeline classes:

- `DataLoader`
- `ModelOptions`
- `DiafiltrationExperiment`
- `UQEngine`

Important design decision:

- `ModelComparator` should not be implemented as a separate top-level class.
- Model comparison should be folded into `UQEngine`.
- The same engine should handle parameter estimation, uncertainty analysis, FIM/DoE analysis, comparison across candidate model options, and recommendation of next experiments.

This organization is preferred because it better matches the intended workflow and avoids creating unnecessary interface complexity just for the sake of OOP.

## Core Responsibilities

### 1. `DataLoader`

Responsible for:

- loading MAT and XLSX experiment files,
- validating inputs,
- normalizing them into a common internal representation,
- optionally applying conductivity-to-concentration transformations where appropriate.

Requirements:

- Do not replace the existing MAT loader.
- Add XLSX support as a first-class loader under the same abstraction.
- The loader abstraction should be reusable and extensible to future formats.

### 2. `ModelOptions`

Responsible for:

- storing one candidate model configuration,
- capturing modeling choices such as transport formulation, thermodynamic options, sigma treatment, permeability formulation, and relevant estimation/DoE settings.

This should function as a lightweight configuration/data class rather than a large behavior-heavy class.

### 3. `DiafiltrationExperiment`

Responsible for:

- binding together one dataset and one `ModelOptions` instance,
- constructing the corresponding Pyomo model,
- preparing Parmest-compatible labeling,
- exposing methods needed for simulation, parameter estimation, and DoE/FIM-related analyses.

This should be the main experiment/model object representing one concrete dataset-model combination.

### 4. `UQEngine`

Responsible for orchestrating the main scientific workflow across one or more datasets and one or more candidate model options.

It should:

- build `DiafiltrationExperiment` instances,
- run parameter estimation,
- compute uncertainty summaries,
- compute FIM and DoE summaries,
- compare candidate model options,
- support model selection, including metrics such as AIC where appropriate,
- support planning the next experiment for either parameter estimation or model discrimination.

This is the key design change: model comparison is part of `UQEngine`, not a separate subsystem.

Additional required capabilities of `UQEngine`:

- Support simultaneous regression across multiple datasets, including a clear treatment of shared parameters versus dataset-specific parameters.
- Support model-selection workflows across candidate model formulations using metrics such as AIC and related fit-versus-complexity criteria.
- Support design of experiments for model discrimination, not just parameter estimation, so the framework can recommend experiments that best distinguish between competing model structures.
- Make the distinction between parameter-estimation workflows and model-discrimination workflows explicit in both interfaces and outputs.

## Core Workflow

The intended workflow should look conceptually like:

`DataLoader -> UQEngine`

Inside `UQEngine`, for each candidate `ModelOptions` object:

1. build the corresponding `DiafiltrationExperiment`,
2. run parameter estimation,
3. analyze uncertainty,
4. analyze FIM / DoE,
5. compare candidate models,
6. recommend next experiments.

## Multi-Dataset Support

The design should explicitly support simultaneous regression across multiple datasets.

At minimum, the design should handle these modes:

- one dataset + many model options,
- many datasets + one model option,
- many datasets + many model options.

This requirement should be discussed during planning before class interfaces are finalized.

The implementation plan should explicitly address:

- how parameters are shared or separated across datasets during joint regression,
- how dataset-specific metadata and experimental conditions are preserved inside the joint workflow,
- how the joint objective is constructed for simultaneous regression,
- how outputs are organized so per-dataset and global regression results are both easy to interpret.

## Plotting Integration

Integrate plotting into the refactored workflow rather than leaving it scattered across scripts.

Requirements:

- Preserve current plotting capability used for DATA1, DATA2, and relevant DATA3 workflows.
- Use current behavior in `unified_codebase_runfile.py`, `unified_codebase_library.py`, `experiment_load_v24.py`, and `conductivity_paper.py` as reference behavior.
- Plot generation should be callable from the refactored workflow or a clearly defined analysis/output layer.

## DATA1 and DATA2 Validation

Use tests and output comparison to validate that the refactor preserves existing DATA1 and DATA2 behavior.

Requirements:

- Validate generated plots and tables against current/public Diafiltration project outputs wherever practical.
- Reuse existing tests and validation utilities where possible.
- Add or update tests only where needed.
- Document any expected differences and explain why they are acceptable.

## Nightly Validation Support

Integrate the refactored workflow into the existing nightly validation framework that now uses pytest.

Requirements:

- Reuse or update current nightly test infrastructure instead of duplicating it.
- Ensure nightly validation exercises the refactored workflow.
- Ensure the nightly validation path is suitable for automated execution.
- Update documentation for how nightly validation works.

## XLSX Support for DATA3

Add `.xlsx` loading support as a first-class capability within the unified framework.

Requirements:

- Do not replace the existing MAT loader.
- Add XLSX support as an additional loader under the same abstraction.
- DATA3 must be able to use the XLSX loader in the refactored workflow.
- Validate that XLSX-loaded data works through downstream stages that DATA3 is expected to support.
- Ensure MAT loading still works after the refactor.

## Non-Negotiable Implementation Constraints

- Do not make unrelated changes.
- Do not remove currently working functionality unless necessary.
- If behavior must change, explain exactly what changed and why.
- Prefer simple, explicit designs over clever abstractions.
- Do not leave placeholder code unless absolutely unavoidable.
- Use atomic, reviewable commits or commit-like logical change groups.
- Add docstrings and type hints where they improve maintainability and clarity.

## Required Workflow for This Task

### Phase 1: Repository Inspection and Planning

Before editing code, inspect and summarize:

- current workflow and entry points,
- current data-loading paths,
- current parmest usage,
- current plotting flow,
- current test coverage,
- current nightly test setup,
- existing reusable unified-framework components,
- main architectural pain points blocking a clean design.

Output required before coding:

- a concrete implementation plan mapped to the requirements above,
- a proposed class/interface design for `DataLoader`, `ModelOptions`, `DiafiltrationExperiment`, and `UQEngine`,
- discussion of multi-dataset regression support,
- discussion of model-selection support, including AIC-style comparison,
- discussion of DoE support for model discrimination,
- discussion of how parameter-estimation DoE and model-discrimination DoE will coexist in the same framework,
- discussion of whether any additional helper classes are truly needed.

The object design may need a few iterations before being locked in.

### Phase 2: Implementation

Execute the refactor in logical steps. For each major step, report:

1. what changed,
2. why it changed,
3. which files were modified,
4. how it was tested or validated.

### Phase 3: Validation

Run and summarize relevant validation:

- unit/integration tests,
- refactor preservation checks for DATA1/DATA2 outputs,
- nightly-test-related validation,
- XLSX-loader validation for DATA3.

### Phase 4: Final Summary

Provide:

- final architecture summary,
- validation summary,
- known limitations or follow-up items.

## Success Criteria

The task is complete only if all of the following are true:

- The unified workflow is refactored into a clear OOP design built around `DataLoader`, `ModelOptions`, `DiafiltrationExperiment`, and `UQEngine`.
- `UQEngine` cleanly orchestrates estimation, uncertainty, FIM/DoE, model comparison, and next-experiment planning.
- The design supports one or more datasets and one or more candidate model options.
- `UQEngine` supports simultaneous regression across multiple datasets.
- `UQEngine` supports model-selection workflows such as AIC-based comparison.
- `UQEngine` supports DoE for model discrimination in addition to parameter-estimation-focused DoE.
- The separation between shared parameters, dataset-specific parameters, and model-option-specific settings is clear in the design.
- Existing parmest-based estimation is preserved or cleanly integrated.
- Plotting is integrated into the refactored design.
- DATA1 and DATA2 validation coverage is updated and runnable.
- Existing nightly validation support is updated to cover the refactor.
- DATA3 can load from XLSX through the new loader abstraction.
- MAT loading still works.
- Documentation is updated where behavior or usage changed.
- Final output includes a clear summary of architecture, modified files, and validation performed.

## Preferred Implementation Style

- Small, reviewable changes
- Minimal duplication
- Clear naming
- Explicit interfaces
- Testable components
- Backward-compatible behavior where practical
