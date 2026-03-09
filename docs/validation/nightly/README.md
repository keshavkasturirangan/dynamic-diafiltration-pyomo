# Nightly Validation Folder

This folder contains all nightly-test-specific validation assets.

## Structure

- `config/known_nonpass.csv`: allowlisted non-pass targets (warning-only in gate).
- `logs/nightly_test_notes.md`: rolling notes from nightly runs.
- `digitized_baselines/`: digitized paper-vs-unified figure comparison inputs/outputs.
- `nightly_ops.md`: operational runbook and commands.

## Related scripts

- `scripts/validation/nightly/nightly_validation_gate.py`
- `scripts/validation/nightly/compare_digitized_baselines.py`
- `scripts/validation/nightly/append_nightly_note.py`

## Related workflow

- `.github/workflows/nightly-validation.yml`
