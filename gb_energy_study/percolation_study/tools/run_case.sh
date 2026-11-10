#!/bin/bash
# run_case.sh N bGB THREADS : start (or resume bit-exactly from checkpoint) one case, tEnd cap 300
N=$1; B=$2; TH=${3:-1}; ROOT=$(cd "$(dirname "$0")/.." && pwd)
D=$ROOT/runs/N${N}_bGB${B}; mkdir -p $D; cd $D
ARGS="bGB=$B outDt=0.1 snapDt=1 saveGrains=0 tEnd=300 outDir=."
if [ -f checkpoint.bin ]; then
  tc=$(python3 -c "import struct;f=open('checkpoint.bin','rb');f.read(20);print(struct.unpack('d',f.read(8))[0])")
  python3 - "$tc" <<'PY'
import sys
tc=float(sys.argv[1])
for fn in ("bubble_metrics.csv","snapshots.csv"):
    L=open(fn).read().splitlines(); keep=[L[0]]+[l for l in L[1:] if float(l.split(',')[1 if fn[0]=='b' else 0])<=tc+1e-9]
    open(fn,'w').write('\n'.join(keep)+'\n')
PY
  echo "OMP_NUM_THREADS=$TH $ROOT/bin/pf3d restartIn=checkpoint.bin $ARGS   # resumed at t=$tc" >> run.sh
  OMP_NUM_THREADS=$TH exec $ROOT/bin/pf3d restartIn=checkpoint.bin $ARGS >> run.log 2>&1
else
  echo "OMP_NUM_THREADS=$TH $ROOT/bin/pf3d initIn=$ROOT/ic/N$N/ic.bin $ARGS" > run.sh
  OMP_NUM_THREADS=$TH exec $ROOT/bin/pf3d initIn=$ROOT/ic/N$N/ic.bin $ARGS > run.log 2>&1
fi
