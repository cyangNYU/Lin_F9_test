import os, sys
from rdkit import Chem
from rdkit.Chem import AllChem, TorsionFingerprints, PandasTools
from rdkit.ML.Cluster import Butina
import pandas as pd


#Based on https://gist.github.com/tdudgeon/b061dc67f9d879905b50118408c30aac
def gen_conformers(mol, numConfs=300, 
                   pruneRmsThresh=0.1, 
                   useExpTorsionAnglePrefs=True, 
                   useBasicKnowledge=True, 
                   enforceChirality=True):
    """Generate conformation with MMFF opt"""
    # generate conf


    ps = AllChem.ETKDG()
    ps.pruneRmsThresh = pruneRmsThresh
    #ps.useExpTorsionAnglePrefs = useExpTorsionAnglePrefs
    #ps.useBasicKnowledge = useBasicKnowledge
    #ps.enforceChirality = enforceChirality
    ps.numThreads = 1
    
    ids = AllChem.EmbedMultipleConfs(mol, numConfs, ps)
    
    # MMFF optimize
    for cid in ids:
        _ = AllChem.MMFFOptimizeMolecule(mol, confId = cid)
    return list(ids)


def calc_energy(mol, conformerId, minimizeIts=300):
    """
    Set minimizeIts to be 0 to turn off min
    since MMFF opt have been done before
    Here, it is used to get MMFF energy
    """
    mp = AllChem.MMFFGetMoleculeProperties(mol)
    ff = AllChem.MMFFGetMoleculeForceField(mol, mp, 
                                           confId=conformerId)

    results = {}
    
    ff.Minimize(maxIts=minimizeIts)
    results["energy_abs"] = ff.CalcEnergy()
    return (results, mol)
    
    
def cluster_conformers(mol, mode="RMSD", threshold=0.2):
    """
    Cluster conf based on heavy atom rmsd 
    Then Butina is used for clustering
    """
    # get heavy atom idx
    heavyatomidx = []
    for a in mol.GetAtoms():
        if a.GetAtomicNum() != 1:
            heavyatomidx.append(a.GetIdx())

    # align on heavy atom for each pair and get dmat
    n = mol.GetNumConformers()
    dmat = []
    for i in range(n-1):
        for j in range(i+1,n):
            dmat.append(Chem.rdMolAlign.AlignMol(mol, mol, i, j,
                                                 atomMap = [(k, k) for k in heavyatomidx]))
    # clustering         
    rms_clusters = Butina.ClusterData(dmat, mol.GetNumConformers(), 
                                      threshold, isDistData=True, reordering=True)
    return rms_clusters
    
    
def runGenerator(input_file, output_file ,numConfs, pruneRmsThresh, 
                 clusterMethod = 'RMSD', clusterThreshold = 2.0,
                 minimizeIterations = 1000):
    """
    Generate conformation as sdf for all smiles in input
    """
    w = Chem.SDWriter(output_file)

    #suppl = Chem.ForwardSDMolSupplier(input_file)
    i = 0
    suppl = Chem.SmilesMolSupplier(input_file, titleLine=False)
    
    for mol in suppl:
        
        i = i+1
        if mol is None: 
            continue
        m = Chem.AddHs(mol)
        # generate the confomers
        conformerIds = gen_conformers(m, numConfs, 
                                      pruneRmsThresh, True, True, True)
        conformerPropsDict = {}
        for conformerId in conformerIds:
            # energy minimise (optional) and energy calculation
            props, m2 = calc_energy(m, conformerId, minimizeIterations)
            conformerPropsDict[conformerId] = props
        
        # cluster the conformers
        rmsClusters = cluster_conformers(m, clusterMethod, clusterThreshold)
        #print(rmsClusters)
        #print ("Molecule", i, ": generated", len(conformerIds), 
        #       "conformers and", len(rmsClusters), "clusters")
        #rmsClustersPerCluster = []
        clusterNumber = 0
        #minEnergy = 9999999999999
        minClusters = []
        for cluster in rmsClusters:
            clusterNumber = clusterNumber+1
            #rmsWithinCluster = align_conformers(m, cluster)
            #print(cluster)
            minCluster = {}
            for conformerId in cluster:
                props = conformerPropsDict[conformerId]
                minCluster[conformerId] = props['energy_abs']
                props["molecule_id"] = i
                props["cluster_no"] = clusterNumber
                props["SMILE"] = Chem.MolToSmiles(mol)
            minClusters.append(min(minCluster, key = minCluster.get))            

        for idx, confId in enumerate(minClusters):
            for name in m.GetPropNames():
                m.ClearProp(name)
            conformerProps = conformerPropsDict[confId]
            m.SetIntProp("conformer_id", idx + 1)
            for key in conformerProps.keys():
                m.SetProp(key, str(conformerProps[key]))
            w.write(m, confId=confId)    
    w.flush()
    w.close()

def main():
    args = sys.argv[1:]

    if not args:
        print ('usage: python Get_confs.py crystal.smi outfile')

        sys.exit(1)

    elif sys.argv[1] == '--help':
        print ('usage: python Get_confs.py crystal.smi outfile')

        sys.exit(1)

    elif sys.argv[1].endswith('.smi'):

        runGenerator(sys.argv[1], sys.argv[2], 300, 0.1) 


if __name__ == '__main__':
    main()
