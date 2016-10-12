#!/bin/sh

# tensortools -f localstress.dat0 --integ x --oformat bin -o outputx.dat0
# tensortools -f outputx.dat0 --integ y --oformat bin -o outputxy.dat0
# tensortools -f outputxy.dat0 --integ z --oformat bin -o outputxyz.dat0
# tensortools -f outputxyz.dat0 --oformat txt -o out.txt

tensortools -f localstress.dat0 --oformat txt -o out.txt

# g_energy_LS -f ener.edr <<EOF
# 31 0
# EOF

# cat out.txt
