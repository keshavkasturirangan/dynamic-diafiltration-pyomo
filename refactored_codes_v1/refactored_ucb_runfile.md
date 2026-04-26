# refactored_ucb_runfile.py

This is the user-facing entry point.

## What it does

- Asks the user to choose `DATA1`, `DATA2`, or `custom`
- Runs the correct workflow from `refactored_ucb_library.py`
- Keeps the menu simple for a non-coder

## User choices

- `DATA1`  
  Recreates the DATA1 paper plots

- `DATA2`  
  Recreates the DATA2 paper plots

- `custom`  
  Lets the user pick one or more `.mat` files and run the model on them

## How to run

In Spyder or from the command line, run:

```bash
python refactored_codes_v1/refactored_ucb_runfile.py
```

## Important behavior

- The runfile automatically calls the library
- If a custom file stores conductivity instead of concentration, the file is converted first
- The runfile does not contain the scientific model itself

## Simple rule

This file should stay very short and easy to read.
It is only the doorway into the library.

