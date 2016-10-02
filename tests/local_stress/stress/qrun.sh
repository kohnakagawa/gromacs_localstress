#!/bin/sh

mdrun_LS -s ../lstress.tpr -rerun ../trj_centered10.trr -lsfd ccfd -lscont all -localsgrid 0.1 -lsgridx 1 -lsgridy 1 -lsgridz 83
# mdrun_LS -s ../lstress.tpr -rerun ../trj_centered10.trr -lsfd fcd -lscont all -localsgrid 0.1 -lsgridx 1 -lsgridy 1 -lsgridz 83
# mdrun_LS -s ../lstress.tpr -rerun ../trj_centered10.trr -lsfd hdlm -lscont all -localsgrid 0.1 -lsgridx 1 -lsgridy 1 -lsgridz 83
# mdrun_LS -s ../lstress.tpr -rerun ../trj_centered10.trr -lsfd hdgm -lscont all -localsgrid 0.1 -lsgridx 1 -lsgridy 1 -lsgridz 83
