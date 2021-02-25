#! /bin/bash

for i in $(ls 1eby_E2E)
do
  b=`basename $i .pdb`
  python split_docked.py 1eby_E2E/${i} 1eby_E2E/${b}
  rm 1eby_E2E/${i}
done

c=$(grep 'Affinity' 1eby_E2E/*/*.pdb | sort -nk3 | head -n 1 | awk -F ':' '{print $1}')
cp $c 1eby_E2E_top1.pdb

