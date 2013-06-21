#!/bin/bash
for (( i=1;i<1000;i++)) ; do
time ./a.out
grep -v '^  0.0000' flux.dat>mag.dat.$i
grep -n -v '^  0.0000' flux.dat|gawk '{print $1-2047.5}' > pos.dat.$i
done

