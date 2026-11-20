#!/bin/bash
# When both theta104 runs are done, start the strong-GB-motion test (theta 120, N 33) in their 4 threads,
# then wait for every run and run the analysis.
ROOT=$(cd "$(dirname "$0")/.." && pwd); cd $ROOT
C="outDt=0.1 snapDt=1 rawEvery=50 saveGrains=1 sliceY=64 outDir=."
IC=$ROOT/ic/N33_z20_c0.18.bin
while [ ! -f runs/theta104_bubbles_only/diagnostics.txt ] || [ ! -f runs/theta104_coupled_df0.3/diagnostics.txt ]; do sleep 60; done
setsid nohup tools/run_case.sh coupled_df0.5      2 initIn=$IC bGB=2.0 dfBias=0.5 biasMode=1 $C tEnd=30 >/dev/null 2>&1 &
setsid nohup tools/run_case.sh L10_bubbles_only   1 initIn=$IC bGB=2.0 dfBias=0 LGB=10 $C tEnd=30 >/dev/null 2>&1 &
setsid nohup tools/run_case.sh L10_coupled_df0.1  1 initIn=$IC bGB=2.0 dfBias=0.1 biasMode=1 LGB=10 $C tEnd=25 >/dev/null 2>&1 &
sleep 60
for r in theta128_bubbles_only theta128_coupled_df0.1 N22_bubbles_only N22_coupled_df0.1 coupled_df0.5 L10_bubbles_only L10_coupled_df0.1; do
  while [ ! -f runs/$r/diagnostics.txt ]; do sleep 120; done; done
for d in runs/N22_* runs/theta1* runs/coupled_df0.5 runs/L10_*; do echo "== $d"; grep -E "end time|mass drift|free-energy|c range" $d/diagnostics.txt; done
python3 tools/analyze_robustness.py
