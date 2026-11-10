./phasefield_gb_bubbles3d initIn=../../ic/common_ic.bin bGB=1.5 tEnd=20 outDt=0.1 snapDt=3 outDir=.
./phasefield_gb_bubbles3d restartIn=checkpoint.bin bGB=1.5 outDt=0.1 snapDt=3 tEnd=20 outDir=.
./phasefield_gb_bubbles3d restartIn=checkpoint.bin bGB=1.5 outDt=0.1 snapDt=3 tEnd=40 outDir=.
