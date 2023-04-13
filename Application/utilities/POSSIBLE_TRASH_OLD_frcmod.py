from .utilities import *
# Would it be better to parse through the structure and figure
# out all the parameters I actually need?  I think so, because
# the other option is to continue this huge dump of ALL params
# which could lead to bloat and parameter overwrites...

def GenerateFRCMODFile(mol2file,frcmodfile):
    S.call(f"parmchk2 -i {mol2file} -o {frcmodfile} -f mol2 -pf 2 -a Y -w Y 1> parmchk.log 2> parmchk.err",shell=True)

def ModifyFRCMOD(frcmod):
    bondfile = os.path.join(os.path.dirname(os.path.abspath(__file__)),"bonds.txt")
    newbonddata_lines   = open(bondfile).readlines()
    newbonddata = ""
    for line in newbonddata_lines:
        if line.strip() != "":
            newbonddata+=line

    anglefile = os.path.join(os.path.dirname(os.path.abspath(__file__)),"angles.txt")
    newangledata_lines  = open(anglefile).readlines()
    newangledata = ""
    for line in newangledata_lines:
        if line.strip() != "":
            newangledata+=line

    dihedralfile = os.path.join(os.path.dirname(os.path.abspath(__file__)),"dihedrals.txt")
    newdihedata_lines   = open(dihedralfile).readlines()
    newdihedata = ""
    for line in newdihedata_lines:
        if line.strip() != "":
            newdihedata += line

    torsionfile = os.path.join(os.path.dirname(os.path.abspath(__file__)),"torsions.txt")
    newimptordata_lines = open(torsionfile).readlines()
    newimptordata = ""
    for line in newimptordata_lines:
        if line.strip != "":
            newimptordata += line

    newfrcmod = os.path.join(os.path.dirname(frcmod),"temp_.frcmod")
    with open(frcmod,"r") as infrc:
        with open(newfrcmod,"w") as outfrc:
            for line in infrc.readlines():
                outfrc.write(line)
                if "BOND" in line:
                    if newbonddata != "":
                        outfrc.write(newbonddata+"\n")
                elif "ANGLE" in line:
                    if newangledata != "":
                        outfrc.write(newangledata)
                elif "DIHE" in line:
                    if newdihedata != "":
                        outfrc.write(newdihedata+"\n")
                elif "IMPROPER" in line:
                    if newimptordata != "":
                        outfrc.write(newimptordata+"\n")
    infrc.close()
    outfrc.close()
    C(f"mv {newfrcmod} {frcmod}")
    return

def HybridizeParams(frcmod,mol2,molecule_type):
    backbone_atoms_to_amberff = {"N2" : "NB","N1" : "N","C3" : "C"}
    ModifyFRCMOD(frcmod)
    newmol2 = os.path.join(os.path.dirname(mol2),"temp_.mol2")
    app_atoms=False
    app_bonds=False
    num_bonds_removed=0
    outmol=[]
    with open(mol2,"r") as inmol:
        for line in inmol.readlines():
            if app_bonds:
                if "@" not in line:
                    tmp = line.split()
                    bond_idx=tmp[0]
                    atom1=tmp[1]
                    atom2=tmp[2]
                    bondorder=tmp[3]
                    if atom1 == 1 and atom2 == 3:
                        line=""
                        num_bonds_removed+=1
                    else:
                        line = f"{int(bond_idx)-num_bonds_removed:>6} {atom1:>5} {atom2:>5} {bondorder}\n"
            if app_atoms:
                if "@" not in line:
                    tmp = line.split()
                    atomname=tmp[1]
                    if atomname in backbone_atoms_to_amberff.keys():
                        tmp[5] = backbone_atoms_to_amberff[atomname]
                    line = f"{tmp[0]:>7} {tmp[1]:<4} {tmp[2]:>14} {tmp[3]:>10} {tmp[4]:>10} {tmp[5]:<3} {tmp[6]:>9} {tmp[7]:>3} {tmp[8]:>15}\n"
            if "ATOM" in line:
                app_atoms=True
                app_bonds=False
            if "BOND"in line:
                app_atoms=False
                app_bonds=True
            outmol.append(line)
        ### Add connections to previous and next residues.
        if molecule_type == "PROTEIN":
            outmol.append("@<TRIPOS>HEADTAIL\nN1 1\nC3 1\n@<TRIPOS>RESIDUECONNECT\n1 N C3 0 0 0 0\n")
        elif molecule_type == "RNA" or molecule_type == "DNA":
            outmol.append("@<TRIPOS>HEADTAIL\nPA 1\nO3' 1\n@<TRIPOS>RESIDUECONNECT\n1 PA O3' 0 0 0 0\n")
        tmpedit = outmol[2].split()
        tmpedit[1] = int(tmpedit[1])-num_bonds_removed
        outmol[3] = " ".join([f"{x:>5}" for x in tmpedit])+"\n"
    with open(newmol2,"w") as of:
        for line in outmol:
            of.write(line)
    C(f"mv {mol2} old_mol2")
    C(f"mv {newmol2} {mol2}")
    return

def ImproveFRCMOD(frcmodfile,leaplog):
    with open(frcmodfile) as f:
        frcmodlines = f.readlines()
    with open(frcmodfile,"w") as f:
        for line in frcmodlines:
            f.write(line)
            if line.strip() == "BOND":
                f.write(GetMissingBonds(leaplog))
            if line.strip() == "ANGLE":
                f.write(GetMissingAngles(leaplog))
            if line.strip() == "DIHE":
                f.write(GetMissingDihedrals(leaplog))
    Self_Improve(leaplog)