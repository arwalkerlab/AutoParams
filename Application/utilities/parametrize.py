from .mol2 import *
from .frcmod import *
from .utilities import *
# from .OLD_frcmod import HybridizeParams
from .OLD_mol2 import *

leaprc_dict = {"DNA":"source leaprc.DNA.OL15\n",
               "RNA":"source leaprc.RNA.OL3\n",
               "Protein":"source leaprc.protein.ff19SB\n",
               "Carbohydrate":"source leaprc.GLYCAM_06j-1\n"}
def CreateMol2File(pdbfile,mol2file,respcharges):
    ### Generate the mol2 file.
    GenerateMol2File(pdbfile,mol2file)
    ### Update the mol2 file with resp charges.
    writeMol2(respcharges, mol2file)

def GetMissingParams(moltype,resname,mol2,leaplog):
    # Write the tleap.in input file with the appropriate BASE forcefield for the moltype (RNA,DNA,Protein,Carbohydrate,etc.)
    with open("tleap.in","w") as f:
        f.write(leaprc_dict[moltype])
        f.write(f"{resname} = loadmol2 {mol2}\n")
        f.write(f"test = loadpdb param.pdb\n")
        f.write("check test\n")
        f.write("quit\n")

    # Run tleap to obtain any missing parameters.
    S.call(f"tleap -f tleap.in > {leaplog}",shell=True)

    # Parse the tleap logfile to get the collection of missing parameters.
    missing_params={"BONDS":[],
                    "ANGLES":[],
                    "DIHEDRALS":[],
                    "TORSIONS":[]}
    for i,line in enumerate(open(leaplog).readlines()):
        if "No sp2 improper torsion term for" in line:
            missing_params["TORSIONS"].append(line.split()[-1])

        elif "Could not find bond parameter for atom types:" in line:
            [a1,a2] = line.split(":")[-1].split("-")
            a1 = a1.strip()
            a2 = a2.strip()
            [a1,a2] = sorted([a1,a2])            
            missing_params["BONDS"].append(f"{a1:<2}-{a2:<2}")

        elif "Could not find angle parameter for atom types:" in line:
            [a1,a2,a3] = line.split(":")[-1].split("-")
            a1 = a1.strip()
            a2 = a2.strip()
            a3 = a3.strip()
            [a1,a3] = sorted([a1,a3])            
            missing_params["ANGLES"].append(f"{a1:<2}-{a2:<2}-{a3:<2}")

        elif "** No torsion terms for" in line:
            [a1,a2,a3,a4] = line.split(":")[-1].split("-")
            if [a3,a2] == sorted([a2,a3]):
                [a1,a2,a3,a4] = [a4,a3,a2,a1]
            missing_params["DIHEDRALS"].append(f"{a1:<2}-{a2:<2}-{a3:<2}-{a4:<2}")
    
    # Reduce the lists to unique sets of terms - no need to define the same parameter multiple times.
    for key in missing_params.keys():
        missing_params[key] = list(set(missing_params[key]))

    return missing_params

def GenerateParameters(file_list,respcharges,moltype,resname):
    CreateMol2File(file_list["Working PDB"],file_list["MOL2"],respcharges)

    ### Run tleap with just the mol2 to figure out what parameters are missing.
    missing_params = GetMissingParams(moltype,resname,file_list["MOL2"],file_list["LeapLog"])
    
    ### Generate FRCMOD from mol2
    GenerateFRCMODFile(file_list["MOL2"],file_list["FRCMOD"],missing_params)

    ### Modify/Update Mol2?
    # HybridizeParams(file_list["FRCMOD"],file_list["MOL2"],moltype)

    return