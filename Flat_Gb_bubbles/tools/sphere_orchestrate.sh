#!/bin/bash
# For each dfBias: wait until the no-bubble sphere run has finished, then start the bubble run (3 threads).
ROOT=$(cd "$(dirname "$0")/.." && pwd); cd $ROOT
for df in 0 0.025 0.05 0.1; do
  ( while [ ! -f runs/sphere_nobub_df$df/diagnostics.txt ]; do sleep 60; done
    BIN=$ROOT/src/pf_mig_96 tools/run_case.sh sphere_bub_df$df 3 initIn=$ROOT/ic/sphere_R32_N25_c0.18.bin sphereR=32 bGB=2.0 dfBias=$df biasMode=1 outDt=0.1 snapDt=1 rawEvery=50 saveGrains=1 stopAfterVanish=10 tEnd=300 outDir=. ) &
done
wait
