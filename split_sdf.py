### split conformers.sdf to individual.pdb ###

import sys, os
from openbabel import pybel

def split_conformers(filename, outdir='.'):
    '''
    split conformers into individual one

    '''
    mols = pybel.readfile('sdf', filename)
    os.makedirs(outdir, exist_ok=True)
    outname = os.path.basename(filename).split('.')[0]

    for count, mol in enumerate(mols, 1):
        output = pybel.Outputfile("mol2", "%s/%s_%02d.mol2"%(outdir, outname, count), overwrite='True')
        output.write(mol)
        output.close()

def main():
    args = sys.argv[1:]
    if not args:
        print ('usage: python split_sdf.py filename outdir')

        sys.exit(1)

    elif sys.argv[1] == '--help':
        print ('usage: python split_sdf.py filename outdir')

        sys.exit(1)

    elif sys.argv[1].endswith('.sdf'):
        if len(args) == 1:
            split_conformers(sys.argv[1])
        elif len(args) == 2:
            split_conformers(sys.argv[1], sys.argv[2])
        else:
            sys.exit(1)
            
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()
    

        
    
