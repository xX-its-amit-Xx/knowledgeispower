#!/usr/bin/env bash
# Post-hoc v9 addendum: dev screen first, then each holdout whose preregistered result exists (waits for it).
cd "$(dirname "$0")/.."
export PIMMS_PYTHON="${PIMMS_PYTHON:-/tmp/claude-0/-home-user-knowledgeispower/19f0adc1-0e52-561e-b422-c31fec9b93e5/scratchpad/pimms/venv/bin/python}"
for t in petrosius2023 leduc2022 weiss2026_liver gregoire2023_mixCTRL px063590_ADbrain sadiku2025; do
  [ -f "results/v9/addendum/dev_$t.json" ] && continue
  (cd benchmarks && nice -n 10 ../.venv/bin/python addendum_v9.py --dev "$t") 2>&1 | grep -v -i warning
done
for h in derks2022 karagach2025_S1 px061710_GSTP1KO px061710_MDAfix px061710_RKOfix khan2023 brunner2022 karagach2025_S2 px077946_H2O2; do
  [ -f "results/v9/addendum/$h.json" ] && continue
  until [ -f "results/confirm_v9/$h.json" ]; do sleep 60; done
  echo "=== addendum $h start $(date -u +%FT%TZ)"
  (cd benchmarks && nice -n 10 ../.venv/bin/python addendum_v9.py "$h" --jobs ${JOBS:-1}) 2>&1 | grep -v -i warning
  echo "=== addendum $h end $(date -u +%FT%TZ)"
done
echo "=== ADDENDUM ALL DONE"
