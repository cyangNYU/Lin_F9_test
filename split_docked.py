import sys, os

def split_docked_file(filename, outdir='.'):
    f1 = open(filename, 'r')
    count = 1
    os.makedirs(outdir, exist_ok=True)
    outname = os.path.basename(filename).split('.')[0]
    for line in f1.readlines():
        if line.startswith('MODEL'):
            f2 = open('%s/%s_%02d.pdb'%(outdir,outname,count),'w')
            count = count+1
        if line.startswith(('ATOM','HETATM','REMARK','CONECT','MASTER')):
            f2.write(line)
        elif line.startswith('ENDMDL'):
            f2.close()

def main():
    args = sys.argv[1:]
    if not args:
        print ('usage: python split_docked.py filename outdir')

        sys.exit(1)

    elif sys.argv[1] == '--help':
        print ('usage: python split_docked.py filename outdir')

        sys.exit(1)

    elif sys.argv[1].endswith('.pdb'):
        if len(args) == 1:
            split_docked_file(sys.argv[1])
        elif len(args) == 2:
            split_docked_file(sys.argv[1], sys.argv[2])
        else:
            sys.exit(1)
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()
    
