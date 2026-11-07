#!/usr/bin/env bash
# Runs every v9 holdout once, in PREREGISTRATION_v9.md order; skips holdouts whose result already exists.
cd "$(dirname "$0")/.."
for h in derks2022 karagach2025_S1 px061710_GSTP1KO px061710_MDAfix px061710_RKOfix khan2023 brunner2022 karagach2025_S2 px077946_H2O2; do
  [ -f "results/confirm_v9/$h.json" ] && continue
  echo "=== $h start $(date -u +%FT%TZ)"
  (cd benchmarks && ../.venv/bin/python confirm_v9.py "$h" --jobs 3) 2>&1 | grep -v -i warning
  echo "=== $h end $(date -u +%FT%TZ)"
done
echo "=== ALL DONE"
