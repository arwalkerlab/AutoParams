from ..utilities import *

def GetMissingParams(moltype,resname,mol2,leaplog,frcmod="",prmtop="test.prmtop",inpcrd="test.inpcrd",pdb="param.pdb"):
    # Write the tleap.in input file with the appropriate BASE forcefield for the moltype (RNA,DNA,Protein,Carbohydrate,etc.)
    with open("tleap.in","w") as f:
        f.write(LEAPRC_DICT[moltype])
        f.write("source leaprc.water.tip3p\n")
        f.write(f"{resname} = loadmol2 {mol2}\n")
        if all([frcmod != "",G(frcmod)]):
            f.write(f"loadamberparams {frcmod}\n")
        f.write(f"test = loadpdb {pdb}\n")
        f.write("addions test Na+ 0\n")
        f.write("addions test Cl- 0\n")
        f.write("solvatebox test TIP3PBOX 12.0\n")
        f.write(f"saveamberparm test {prmtop} {inpcrd}\n")
        f.write("run\n")
        f.write("quit\n")

    # Run tleap to obtain any missing parameters.
    S.call(f"tleap -f tleap.in > {leaplog}",shell=True)
    if all([G(prmtop),G(inpcrd)]):
        return False
    # Parse the tleap logfile to get the collection of missing parameters.
    missing_params={"MASSES":[],
                    "BONDS":[],
                    "ANGLES":[],
                    "DIHEDRALS":[],
                    "TORSIONS":[],
                    "NONBONDED":[]}
    for line in open(leaplog).readlines():
        if "No sp2 improper torsion term for" in line:
            missing_params["TORSIONS"].append(line.split()[-1])
            print("Found missing torsion term:",line.split()[-1])

        elif "Could not find bond parameter" in line:
            [a1,a2] = line.split(":")[-1].split("-")
            a1 = a1.strip()
            a2 = a2.strip()
            [a1,a2] = sorted([a1,a2])            
            missing_params["BONDS"].append(f"{a1:<2}-{a2:<2}")

        elif "Could not find angle parameter" in line:
            [a1,a2,a3] = line.split(":")[-1].split("-")
            a1 = a1.strip()
            a2 = a2.strip()
            a3 = a3.strip()
            [a1,a3] = sorted([a1,a3])            
            missing_params["ANGLES"].append(f"{a1:<2}-{a2:<2}-{a3:<2}")

        elif "** No torsion terms for atom types:" in line:
            [a1,a2,a3,a4] = line.replace("** No torsion terms for atom types:","").split("-")
            a1 = a1.strip()
            a2 = a2.strip()
            a3 = a3.strip()
            a4 = a4.strip()
            if [a3,a2] == sorted([a2,a3]):
                [a1,a2,a3,a4] = [a4,a3,a2,a1]
            missing_params["DIHEDRALS"].append(f"{a1:<2}-{a2:<2}-{a3:<2}-{a4:<2}")
        
        elif "could not find vdW (or other) parameters for type" in line:
            a1 = line.strip().split()[-1].replace("(","").replace(")","")
            missing_params["MASSES"].append(f"{a1:<2}")
            missing_params["NONBONDED"].append(f"{a1:<2}")
    
    # Reduce the lists to unique sets of terms - no need to define the same parameter multiple times.
    for key in missing_params.keys():
        missing_params[key] = list(set(missing_params[key]))

    return missing_params
