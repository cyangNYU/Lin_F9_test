#! /bin/bash

mkdir -p 1eby_E2E

for i in $(ls 1eby_conformers)
do
  b=`basename $i .mol2`
  ./smina.static -r 1eby_protein.pdbqt -l 1eby_conformers/${i} \
   --autobox_ligand 1eby_ligand.pdbqt --scoring Lin_F9 -o 1eby_E2E/${b}.pdb 
done
