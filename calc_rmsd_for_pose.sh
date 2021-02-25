#! /bin/bash

obabel -ipdb $1 -omol2 -Otest.mol2 

./DockRMSD test.mol2 $2 | grep 'RMSD' 
rm test.mol2
