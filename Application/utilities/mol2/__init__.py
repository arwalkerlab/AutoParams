from ..utilities import *

def MAKE_MOL2_FILE(pdbfile,mol2file,respcharges,connections):
    ### Generate the mol2 file.
    GenerateMol2File(pdbfile,mol2file)
    ### Update the mol2 file with resp charges.
    writeMol2(respcharges, mol2file,connections)

def GenerateMol2File(pdbfile,mol2file):
    S.call(f"antechamber -i {pdbfile} -fi pdb -at amber -o {mol2file} -fo mol2 -pf y 1> antechamber.log 2> antechamber.err",shell=True)

def writeMol2(charges, filename,connections):
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
        if connections != ["0","0"]:
            w.write("@<TRIPOS>HEADTAIL\n")
            if connections[0] == "0":
                w.write("0 0\n")
            else:
                w.write(f"{connections[0]} 1\n")
            if connections[1] == "0":
                w.write("0 0\n")
            else:
                w.write(f"{connections[1]} 1\n")
            w.write("@<TRIPOS>RESIDUECONNECT\n")
            w.write(f"1 {connections[0]} {connections[1]} 0 0 0 0\n")
            w.write("\n")
    S.call(f"mv temp.txt {filename}",shell=True)
    

