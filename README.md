# Lin_F9 tutorial

*Lin_F9 is a linear empirical scoring function for protein-ligand docking.*

<img src="/home/cyang/Lin_F9/tutorial/Lin_F9_test/plot.png" alt="plot" style="zoom:60;" />

### Preparation

Download a fork of Smina docking suite with Lin_F9 as an optional built-in scoring function at <https://github.com/cyangNYU/Lin_F9_test>

Take one protein-ligand complex (PDB_id: 1eby) as an example:

Initial files of 1eby_ligand.mol2 and 1eby_protein.pdb obtained from PDBbind database

MGLTools are used for preparing PDBQT files: 1eby_ligand.pdbqt and 1eby_protein.pdbqt.

(MGLTools can be downloaded at http://mgltools.scripps.edu/downloads)



### **Lin_F9 is a linear empirical scoring function for protein-ligand docking.Local optimization**

```shell
./smina.static -r 1eby_protein.pdbqt -l 1eby_ligand.pdbqt --local_only --scoring Lin_F9 -o 1eby_optimized.pdb 
```

Print the score of locally optimized pose.

```shell
grep 'Affinity' 1eby_optimized.pdb | awk '{print $3}'
```



### **Flexible re-docking**

```shell
./smina.static -r 1eby_protein.pdbqt -l 1eby_ligand.pdbqt --autobox_ligand 1eby_ligand.pdbqt --scoring Lin_F9 -o 1eby_flexRedock.pdb
```

Print the score of the best-scored pose

```shell
grep -m1 'Affinity' 1eby_flexRedock.pdb | awk '{print $3}'
```

Calculate the RMSD between the best-scored pose and crystal ligand pose 

 (DockRMSD is used to calculate RMSD, can be downloaded at https://zhanglab.ccmb.med.umich.edu/DockRMSD/)

```shell
python split_docked.py 1eby_flexRedock.pdb 1eby_flexRedock
./calc_rmsd_for_pose.sh 1eby_flexRedock/1eby_flexRedock_01.pdb 1eby_ligand.mol2
```



### **End-to-end (E2E) docking**

##### 1. OpenBabel 2.4.1 used to generate maximum 10 conformers per ligand.

```shell
babel -imol2 1eby_ligand.mol2 -osdf 1eby_ligand_2D.sdf --gen2D 

babel -isdf 1eby_ligand_2D.sdf -osdf 1eby_ligand_3D.sdf --gen3D 

obabel 1eby_ligand_3D.sdf -O conformers.sdf --conformer --nconf 10 --score rmsd --writeconformers
```

##### 2. All conformers are sequentially docked to protein receptor structure using flexible docking.

```shell
python split_sdf.py conformers.sdf 1eby_conformers/
./E2E_docking.sh
```

Print the score of the best-scored pose

```shell
grep -m1 'Affinity' 1eby_E2E/*.pdb | sort -nk3 | head -n 1 | awk '{print $3}'
```

Calculate the RMSD between the best-scored pose and crystal ligand pose

```shell
./GetSplit.sh
./calc_rmsd_for_pose.sh 1eby_E2E_top1.pdb 1eby_ligand.mol2
```

