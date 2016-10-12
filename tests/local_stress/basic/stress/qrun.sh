#!/bin/sh

export GMX_MAXBACKUP=-1 # We do not need backup files.

grompp_LS -f ../lstress.mdp -c ../step7_700.gro -p ../topol.top -o ../lstress.tpr

# mdrun_LS -s ../lstress.tpr -rerun ../trj_centered10.trr -lsfd ccfd -lscont all -localsgrid 0.1 -lsgridx 1 -lsgridy 1 -lsgridz 83
# mdrun_LS -s ../lstress.tpr -rerun ../trj_centered10.trr -lsfd fcd -lscont all -localsgrid 0.1 -lsgridx 1 -lsgridy 1 -lsgridz 83
# mdrun_LS -s ../lstress.tpr -rerun ../trj_centered10.trr -lsfd hdlm -lscont angles -localsgrid 0.1 -lsgridx 1 -lsgridy 1 -lsgridz 83
# mdrun_LS -s ../lstress.tpr -rerun ../trj_centered10.trr -lsfd hdgm -lscont angles -localsgrid 0.1 -lsgridx 1 -lsgridy 1 -lsgridz 83
