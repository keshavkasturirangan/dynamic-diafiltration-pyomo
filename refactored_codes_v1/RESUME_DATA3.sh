#!/bin/bash
# ============================================================================
# RESUME_DATA3.sh — continue everything this session was running, in a new thread.
# IDEMPOTENT: every campaign skips already-completed work (WORKER_DONE / done nodes),
# so it is safe to run repeatedly and to interrupt/restart.
#   Usage:  bash RESUME_DATA3.sh            # run the remaining campaigns + rebuild decks
#           bash RESUME_DATA3.sh monitor    # just print a status snapshot and exit
# Run from refactored_codes_v1/.  Heavy: hours of DAE solves; uses ~8 cores.
# ============================================================================
cd "$(dirname "$0")" || exit 1
B="../UnifiedFramework/DATA3/results/paper_artifacts/nf270/bform_study"

status() {
  echo "=== DATA3 status $(date '+%H:%M') ==="
  echo "per-B-form slice : $(find $B/stack_contours -path '*/sigLp/slice/stacked.png' 2>/dev/null | wc -l | tr -d ' ')/24"
  echo "per-B-form profile: $(find $B/stack_contours -path '*/sigLp/profile/stacked.png' 2>/dev/null | wc -l | tr -d ' ')/24"
  echo "deck profile_contours: $(find $B/profile_contours -name WORKER_DONE 2>/dev/null | wc -l | tr -d ' ')/27"
  echo "running procs: bformsweep=$(pgrep -f bformsweep | wc -l | tr -d ' ') campaign=$(pgrep -f '_run_profile_contours.py run' | wc -l | tr -d ' ') ipopt=$(pgrep -f ipopt | wc -l | tr -d ' ')"
  echo "decks: $(ls -1 DATA3_STORY_deck.pptx DATA3_FITTING_4experiments.pptx 2>/dev/null | tr '\n' ' ')"
}

if [ "$1" = "monitor" ]; then status; exit 0; fi

# stop any stale leftovers from a previous session so we start clean
pkill -f "_run_profile_contours.py" 2>/dev/null; pkill -f "bformsweep" 2>/dev/null; pkill -f ipopt 2>/dev/null; sleep 2

echo ">>> [1/4] per-B-form σ×Lp stacking — 6 forms (const,linear,quad,cubic,sat,Donnan) × slice+profile × 4 groups"
for m in slice profile; do
  for g in CaCl2 NaCl_diluting NaCl_concentrating LaCl3; do
    python3 -u _run_profile_contours.py bformsweep "$g" 9 "$m"
  done
done

echo ">>> [2/4] constant-B 3-plane stacking (σ×Lp, B×Lp, σ×B) per group  [resumes / re-renders]"
for g in CaCl2 NaCl_diluting NaCl_concentrating LaCl3; do
  python3 -u _run_profile_contours.py stack "$g" single 11 profile
done

echo ">>> [3/4] per-experiment deck contours (4 headline + others), resumes from $(find $B/profile_contours -name WORKER_DONE 2>/dev/null | wc -l | tr -d ' ')/27"
python3 -u _run_profile_contours.py campaign deck 15 8

echo ">>> [4/4] rebuild both decks with all completed figures"
python3 _build_story_deck.py
python3 _build_fitting_deck.py

status
echo ">>> DATA3 resume complete. Decks: DATA3_STORY_deck.pptx (full story) + DATA3_FITTING_4experiments.pptx (4-experiment fitting)."
