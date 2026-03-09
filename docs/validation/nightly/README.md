# Nightly Validation Folder

This folder contains all nightly-test-specific validation assets.

Nightly status should be interpreted using the main validation semantics:
- numeric `PASS` / `FAIL` apply only to published-paper targets with completed numeric comparisons;
- mapped/page-matched/locked figure targets are evidence of published-paper correspondence, not standalone numeric passes.

## Structure

- `config/known_nonpass.csv`: allowlisted non-pass targets (warning-only in gate).
- `logs/nightly_test_notes.md`: rolling notes from nightly runs.
- `digitized_baselines/`: digitized paper-vs-unified figure comparison inputs/outputs.
- `nightly_ops.md`: operational runbook and commands.
- `slides/nightly_ops_10_slide_deck.pptx`: optional briefing deck derived from `nightly_ops.md`.

## Related scripts

- `scripts/validation/nightly/nightly_validation_gate.py`
- `scripts/validation/nightly/compare_digitized_baselines.py`
- `scripts/validation/nightly/append_nightly_note.py`
- `scripts/validation/nightly/build_nightly_ops_deck.py` (regenerates the slide deck)

## Related workflow

- `.github/workflows/nightly-validation.yml`

## Related migration note

- `docs/validation/unified_codex_code_reorg_map.md`
