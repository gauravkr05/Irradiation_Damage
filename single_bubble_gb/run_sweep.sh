#!/usr/bin/env bash
# run_sweep.sh
#
# Parameter sweeps for the single-bubble growth study.
# Edit the lists below, then:   bash run_sweep.sh A     (or B, C, D, E)
# Each run writes its own folder; analyse with analyze_growth.py.
#
# Timing measured on 2 CPU cores, 48^3 box: ~2 ms per step,
# so 30000 steps take about 1 minute. A 64^3 box is ~2.4x slower.

set -e
EXE=./single_bubble_gb
COMMON="N=48 steps=30000 outputEvery=500"
STAGE=${1:-A}

case $STAGE in

  A)  # Sanity check: default paper values, conservation must be ~1e-14
      $EXE $COMMON out=A_default
      python3 analyze_growth.py A_default --out A.png
      ;;

  B)  # Supersaturation: how growth rate depends on matrix concentration
      for c0 in 0.12 0.15 0.20 0.25; do
        $EXE $COMMON Mb=1.0 R0=6 c0=$c0 out=B_c0_$c0
      done
      python3 analyze_growth.py B_c0_* --out B.png
      ;;

  C)  # Critical radius at fixed c0: smallest R0 that still grows
      C0=0.12
      for r in 3 4 5 6 8; do
        $EXE $COMMON Mb=1.0 c0=$C0 R0=$r out=C_R0_$r
      done
      python3 analyze_growth.py C_R0_* --out C.png
      ;;

  D)  # Mobility ratio: bulk vs GB/surface supply (changes mechanism)
      for mb in 0.01 0.1 1.0; do
        $EXE $COMMON c0=0.20 R0=6 Mb=$mb out=D_Mb_$mb
      done
      python3 analyze_growth.py D_Mb_* --out D.png
      python3 analyze_growth.py D_Mb_* --physical --out D_physical.png
      ;;

  E)  # Contact angle parameter sets from the paper
      for th in 40 50 60 80; do
        $EXE $COMMON Mb=1.0 c0=0.20 R0=6 theta=$th steps=60000 out=E_theta_$th
      done
      python3 analyze_growth.py E_theta_* --out E.png
      ;;

  *)  echo "usage: bash run_sweep.sh A|B|C|D|E" ;;
esac
