# Pseudocode for Unified Refactor Plan (Option B)

This pseudocode reflects the revised design direction:

- `ModelComparator` is folded into `UQEngine`
- `UQEngine` is the main orchestration layer
- the top-level organization is:
  - `DataLoader`
  - `ModelOptions`
  - `DiafiltrationExperiment`
  - `UQEngine`

This is analysis only. It does not imply any code changes by itself.

```text
BEGIN REFACTOR PLAN

GOAL:
  create a unified workflow for DATA1 / DATA2 / DATA3
  preserve scientific behavior
  support MAT and XLSX loading
  validate against DATA1/DATA2 plots and tables
  support nightly pytest validation

--------------------------------------------------
PHASE 1: INSPECTION AND PLANNING
--------------------------------------------------

inspect repository:
  inspect current entry points
  inspect unified_codebase_library.py
  inspect unified_codebase_runfile.py
  inspect utility.py and legacy behavior
  inspect parmest-related code
  inspect plotting flow
  inspect existing tests
  inspect nightly pytest setup
  inspect MAT/XLSX loading paths
  identify reusable components
  identify architectural pain points

propose revised object organization:
  DataLoader
  ModelOptions
  DiafiltrationExperiment
  UQEngine

decide how to support:
  one dataset + many model options
  many datasets + one model option
  many datasets + many model options
  simultaneous regression with shared vs dataset-specific parameters
  model selection using criteria such as AIC
  DoE for model discrimination
  separation of global outputs vs per-dataset outputs

produce implementation plan before coding

--------------------------------------------------
CLASS: DataLoader
--------------------------------------------------

responsibility:
  load and normalize experimental data

inputs:
  MAT files
  XLSX files
  optional overrides / conversion settings

methods:
  load(input_source)
  load_many(input_sources)
  validate(raw_data)
  normalize(validated_data)
  apply_conductivity_conversion(data, settings)

outputs:
  ExperimentalData object
  or list of ExperimentalData objects

--------------------------------------------------
CLASS: ModelOptions
--------------------------------------------------

responsibility:
  store one candidate model configuration

possible contents:
  transport model choice
  thermodynamic option
  parameterization choice
  sigma treatment
  permeability formulation
  advanced XLSX transport/thermo flags
  estimation settings
  DoE settings

outputs:
  one structured configuration object used to build and analyze a model

--------------------------------------------------
CLASS: DiafiltrationExperiment
--------------------------------------------------

responsibility:
  represent one dataset-model combination

inputs:
  one ExperimentalData object
  one ModelOptions object

methods:
  build_model()
  discretize_model()
  label_for_parmest()
  simulate()
  estimate_parameters()
  compute_fim()
  run_doe()

outputs:
  Pyomo/Parmest/DoE-ready experiment object
  simulation / estimation / analysis outputs for one concrete case

--------------------------------------------------
CLASS: UQEngine
--------------------------------------------------

responsibility:
  serve as the main orchestration and decision-support layer

important design choice:
  model comparison is folded into UQEngine
  not a separate top-level class

inputs:
  one dataset or many datasets
  one ModelOptions or many ModelOptions

main algorithm:

  run(data_or_dataset_list, model_options_list):
    datasets = ensure_list(data_or_dataset_list)
    options_list = ensure_list(model_options_list)
    all_results = []

    FOR each model_options in options_list:
      experiment_list = []

      FOR each dataset in datasets:
        experiment = build_experiment(dataset, model_options)
        experiment_list.append(experiment)

      estimation_result = run_parameter_estimation(experiment_list)
      uq_result = analyze_parameter_uncertainty(experiment_list, estimation_result)
      fim_result = analyze_fim(experiment_list, estimation_result)
      doe_result = analyze_doe(experiment_list, estimation_result)

      candidate_result = assemble_candidate_result(
        model_options,
        experiment_list,
        estimation_result,
        uq_result,
        fim_result,
        doe_result
      )

      all_results.append(candidate_result)

    comparison_result = compare_candidates(all_results)
    recommendation_result = propose_next_experiments(all_results, comparison_result)

    RETURN final_result_bundle(
      all_results,
      comparison_result,
      recommendation_result
    )

methods:
  build_experiment(dataset, model_options)
  run_parameter_estimation(experiment_list)
  analyze_parameter_uncertainty(experiment_list, estimation_result)
  analyze_fim(experiment_list, estimation_result)
  analyze_doe(experiment_list, estimation_result)
  compare_candidates(all_results)
  propose_next_experiments(all_results, comparison_result)

--------------------------------------------------
DETAIL: UQEngine methods
--------------------------------------------------

build_experiment(dataset, model_options):
  create one DiafiltrationExperiment from dataset + options
  return experiment

run_parameter_estimation(experiment_list):
  perform parmest-based regression
  support simultaneous regression over multiple datasets where needed
  explicitly handle shared parameters vs dataset-specific parameters
  return fitted parameters, global objective summaries, per-dataset residual summaries, predictions

analyze_parameter_uncertainty(experiment_list, estimation_result):
  compute covariance / standard deviations / correlations / confidence intervals
  using parmest covariance and/or inverse FIM where available
  return uncertainty summary

analyze_fim(experiment_list, estimation_result):
  compute Jacobian, FIM, eigenstructure, identifiability metrics
  return FIM summary

analyze_doe(experiment_list, estimation_result):
  evaluate experimental informativeness
  support parameter-estimation-focused DoE
  support model-discrimination-focused DoE
  return DoE summary with the purpose of the recommended experiment made explicit

compare_candidates(all_results):
  compare candidate model options
  support model selection metrics such as AIC where appropriate
  compare fit quality, uncertainty, identifiability, and discrimination performance
  return comparison summary

propose_next_experiments(all_results, comparison_result):
  use FIM / DoE / model comparison outputs
  propose next experimental conditions
  distinguish between:
    next experiment for better parameter estimation
    next experiment for model discrimination
  return recommendation summary

--------------------------------------------------
PLOTTING INTEGRATION
--------------------------------------------------

plotting should not remain scattered

conceptual behavior:
  plots can be triggered from the workflow or a dedicated reporting/output layer

plots should include:
  simulation vs experimental data
  DATA1 / DATA2 reproduction plots
  candidate-model comparison plots where needed
  DoE / uncertainty summaries where useful

--------------------------------------------------
PHASE 2: IMPLEMENTATION
--------------------------------------------------

implement in small logical steps:
  step 1:
    establish target object structure
  step 2:
    refactor loading into DataLoader abstraction
  step 3:
    define ModelOptions as structured configuration
  step 4:
    define DiafiltrationExperiment around one dataset + one model option
  step 5:
    move estimation / FIM / DoE / comparison orchestration into UQEngine
  step 6:
    integrate plotting cleanly
  step 7:
    update tests and nightly validation
  step 8:
    add and validate DATA3 XLSX loading path
  step 9:
    update documentation

for each step, report:
  what changed
  why it changed
  which files changed
  how it was tested / validated

--------------------------------------------------
PHASE 3: VALIDATION
--------------------------------------------------

validate:
  unit tests
  integration tests
  DATA1 output reproduction
  DATA2 output reproduction
  nightly pytest workflow
  DATA3 MAT workflow
  DATA3 XLSX workflow

compare against:
  existing plots
  existing tables
  existing fit behavior
  existing scientific conclusions where practical

document any expected differences

--------------------------------------------------
PHASE 4: FINAL SUMMARY
--------------------------------------------------

provide:
  final architecture summary
  explanation of object responsibilities
  summary of modified files
  summary of validation performed
  known limitations or follow-up items

--------------------------------------------------
SUCCESS CONDITIONS
--------------------------------------------------

success if:
  DataLoader works for MAT and XLSX
  ModelOptions cleanly represents candidate model choices
  DiafiltrationExperiment cleanly represents one dataset-model pairing
  UQEngine orchestrates estimation, uncertainty, FIM/DoE, model comparison, and next-experiment planning
  simultaneous regression across multiple datasets is explicitly supported
  model selection using metrics such as AIC is supported
  DoE for model discrimination is supported
  the design clearly distinguishes shared parameters, dataset-specific parameters, and per-model-option settings
  DATA1 and DATA2 validation is runnable
  DATA3 XLSX workflow is supported
  nightly pytest validation is updated
  documentation reflects the new design

END REFACTOR PLAN
```
