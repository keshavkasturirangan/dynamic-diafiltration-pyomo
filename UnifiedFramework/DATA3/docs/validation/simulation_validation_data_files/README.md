# simulation_validation_data_files

One-time generated simulation-side CSV files used to validate legacy paper figures.

Source code paths:
- `utility.py`
- `DATA1_model_demo.ipynb`
- `DATA2_model_demo.ipynb`
- `DATA2_visualization.ipynb`
- `run_DATA2_model_variations.py`

These CSVs are intended to stay stable after generation so later pytest checks can compare unified-code outputs against the same published-paper baselines.

This folder is the current source of truth for simulation-side paper-validation CSVs.
