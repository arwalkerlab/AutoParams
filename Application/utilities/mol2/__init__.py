from ..utilities import *

def GenerateMol2File(pdbfile,mol2file):
    S.call(f"antechamber -i {pdbfile} -fi pdb -at amber -o {mol2file} -fo mol2 -pf y 1> antechamber.log 2> antechamber.err",shell=True)

def writeMol2(charges, filename):
    lines = open(filename, 'r').readlines()
    i = 0
    start = 0
    end = 0
    for line in lines:
        if "@<TRIPOS>ATOM" in line:
            start = i + 1
        if "@<TRIPOS>BOND" in line:
            end = i
            break
        i = i + 1
    if end - start != len(charges):
        raise RuntimeError("Size of molecule in mol2 file differs from size of molecule in tcout file")
    with open("temp.txt", 'w') as w:
        for line in lines[:start]:
            w.write(line.replace("DU","NB"))
        for line, charge in zip(lines[start:end],charges):
            oldCharge = line.split()[8]
            w.write(line.replace(oldCharge,f"{float(charge):>9.06f}").replace("DU","NB"))
        for line in lines[end:]:
            w.write(line.replace("DU","NB"))
    S.call(f"mv temp.txt {filename}",shell=True)
    

