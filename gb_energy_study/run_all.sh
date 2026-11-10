#!/bin/bash
# Reproducible driver for the GB-energy study (serial).
# Safe to re-run at any time: finished stages are skipped and an interrupted
# 3-D run is resumed bit-exactly from its last checkpoint (metrics rows written
# after the checkpoint are discarded first, so nothing is duplicated).
set -e
ROOT=$(cd "$(dirname "$0")" && pwd); cd "$ROOT"
while pgrep -x phasefield_gb_b > /dev/null || pgrep -x pf2d > /dev/null; do sleep 5; done
CXX="g++ -O3 -march=native -fopenmp -std=c++17"
[ -x src/phasefield_gb_bubbles3d ] || $CXX src/phasefield_gb_bubbles3d.cpp -o src/phasefield_gb_bubbles3d
[ -x src/pf2d ] || $CXX -DGRID_NX=160 -DGRID_NY=1 -DGRID_NZ=96 src/phasefield_gb_bubbles3d.cpp -o src/pf2d
# common initial condition, generated ONCE (same seed-11 centres as the original code)
[ -f ic/common_ic.bin ] || src/phasefield_gb_bubbles3d initOut=ic/common_ic.bin initOnly=1 r0=4.5 cInit=0.12 outDir=ic
CASES="1.2 1.5 2.0 3.0"
T1=20; T2=${T2:-40}
endtime() { [ -f "$1/diagnostics.txt" ] && awk '/^end time/{print $4}' "$1/diagnostics.txt" || echo 0; }
ckpttime() { python3 -c "import struct;f=open('$1/checkpoint.bin','rb');f.read(20);print(struct.unpack('d',f.read(8))[0])"; }

# run_case DIR TEND ARGS... : start from the common IC, or resume from the checkpoint
run_case() {
  local d=$1 tend=$2; shift 2
  [ "$(endtime $d)" == "$tend" ] && return 0
  if [ -f $d/checkpoint.bin ]; then
    local tc=$(ckpttime $d)
    python3 - "$d" "$tc" <<'PY'
import sys, pandas as pd
d, tc = sys.argv[1], float(sys.argv[2])
for fn in ("bubble_metrics.csv", "snapshots.csv"):
    m = pd.read_csv(f"{d}/{fn}")
    m[m.time_dimensionless <= tc + 1e-9].to_csv(f"{d}/{fn}", index=False)
PY
    echo "./phasefield_gb_bubbles3d restartIn=checkpoint.bin $* tEnd=$tend outDir=." >> $d/run.sh
    (cd $d && ./phasefield_gb_bubbles3d restartIn=checkpoint.bin "$@" tEnd=$tend outDir=. >> run.log 2>&1)
  else
    echo "./phasefield_gb_bubbles3d initIn=../../ic/common_ic.bin $* tEnd=$tend outDir=." > $d/run.sh
    (cd $d && ./phasefield_gb_bubbles3d initIn=../../ic/common_ic.bin "$@" tEnd=$tend outDir=. > run.log 2>&1)
  fi
}

# stage 1: main experiment to T1 (identical IC and numerics, only bGB differs)
for b in $CASES; do
  d=cases/bGB_$b; mkdir -p $d; cp -n src/phasefield_gb_bubbles3d src/phasefield_gb_bubbles3d.cpp $d/
  [ "$(endtime $d)" != 0 ] && continue
  run_case $d $T1 bGB=$b outDt=0.1 snapDt=3
done
# stage 2: contact-angle equilibration, 2-D lens, two opposite starting shapes
for s in circle flat; do for b in $CASES; do
  d=angle/${s}_bGB$b; [ -f $d/diagnostics.txt ] && continue
  src/pf2d initIn=angle/ic_$s.bin bGB=$b tEnd=400 outDt=1 snapDt=50 outDir=$d > $d.log 2>&1
done; done
# stage 3: mobility sensitivity (Ms = Mgb = 0), extreme cases, horizon T1
for b in ${MOBCASES:-1.2 3.0}; do
  d=cases_mob0/bGB_$b; mkdir -p $d; cp -n src/phasefield_gb_bubbles3d $d/
  [ "$(endtime $d)" != 0 ] && continue
  run_case $d $T1 bGB=$b Ms=0 Mgb=0 outDt=0.1 snapDt=3 saveGrains=0
done
# stage 4: extend every main case to T2 (bit-exact restart)
for b in $CASES; do
  run_case cases/bGB_$b $T2 bGB=$b outDt=0.1 snapDt=3
done
echo ALLDONE > all_done.flag
