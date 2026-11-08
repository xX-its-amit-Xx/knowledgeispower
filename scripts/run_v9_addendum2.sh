#!/usr/bin/env bash
# Post-hoc addendum 2 (R reference implementations + literature-standard ports) on the primary holdouts, 4 of the 10
# preregistered seeds per mask (runtime), after each holdout's preregistered result exists.
cd "$(dirname "$0")/.."
METHODS=scp_graphfs_mcp.addendum2.rbridge:r_bpca,scp_graphfs_mcp.addendum2.rbridge:r_missforest,scp_graphfs_mcp.addendum2.rbridge:r_knn,scp_graphfs_mcp.addendum2.seqknn:seqknn,scp_graphfs_mcp.addendum2.trknn:trknn,scp_graphfs_mcp.addendum2.impseq:impseq,scp_graphfs_mcp.addendum2.rem:rem_protein,scp_graphfs_mcp.addendum2.simple:median
for h in derks2022 karagach2025_S1 px061710_GSTP1KO px061710_MDAfix px061710_RKOfix; do
  [ -f "results/v9/addendum2/$h.json" ] && continue
  echo "=== addendum2 $h start $(date -u +%FT%TZ)"
  (cd benchmarks && nice -n 10 ../.venv/bin/python addendum2_v9.py "$h" --methods "$METHODS" --seeds 9000,9001,9002,9003) 2>&1 | grep -v -i warning
  echo "=== addendum2 $h end $(date -u +%FT%TZ)"
done
echo "=== ADDENDUM2 ALL DONE"
