# conductivity_paper.py

This file contains the paper-based conductivity models.

## What it does

- Computes conductivity from concentration
- Provides the paper models used for conductivity-based conversion
- Supplies the physics used when a dataset gives conductivity instead of concentration

## Main functions

- `_shedlovsky(...)`
- `variant_shedlovsky(...)`
- `msa_transport(...)`

## How it is used

The refactored library calls this file when it needs to convert conductivity to concentration.

The simple flow is:

1. Read conductivity values from the experimental file
2. Use this file to evaluate the paper model
3. Numerically invert the model to recover concentration

## Notes

- This file is kept as the source of truth for conductivity physics
- The refactored library should call it directly when conversion is needed
- No old linear conversion is used

