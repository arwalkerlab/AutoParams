from ..utilities import *
from ..parameter_lib import *

def GenerateFRCMODFile(base_ff,frcmod,missing_params):
    if not missing_params:
        S.call(f"touch {frcmod}",shell=True)
        return
    with open(frcmod,"w") as f:
        frcmod_name = frcmod.split("/")[-1].replace(".frcmod","")
        f.write(f"{frcmod_name}, use with {base_ff.replace('source','')}\n")
        f.write("MASS\n")
        for missing_mass in missing_params["MASSES"]:
            if missing_mass in ALL_KNOWN_PARAMS["MASSES"].keys():
                f.write(ALL_KNOWN_PARAMS["MASSES"][missing_mass])
            else:
                f.write(NewMass(missing_mass))
        f.write("\n")
        f.write("BOND\n")
        for missing_bond in missing_params["BONDS"]:
            if missing_bond in ALL_KNOWN_PARAMS["BONDS"].keys():
                f.write(ALL_KNOWN_PARAMS["BONDS"][missing_bond])
            else:
                f.write(NewBond(missing_bond))
                pass
        f.write("\n")
        f.write("ANGLE\n")
        for missing_angle in missing_params["ANGLES"]:
            if missing_angle in ALL_KNOWN_PARAMS["ANGLES"].keys():
                f.write(ALL_KNOWN_PARAMS["ANGLES"][missing_angle])
            else:
                f.write(NewAngle(missing_angle))
                pass
        f.write("\n")
        f.write("DIHE\n")
        for missing_dihe in missing_params["DIHEDRALS"]:
            if missing_dihe in ALL_KNOWN_PARAMS["DIHEDRALS"].keys():
                f.write(ALL_KNOWN_PARAMS["DIHEDRALS"][missing_dihe])
            else:
                f.write(NewDihedral(missing_dihe))
                pass
        f.write("\n")
        f.write("IMPROPER\n")
        for missing_torsion in missing_params["TORSIONS"]:
            if missing_torsion in ALL_KNOWN_PARAMS["TORSIONS"].keys():
                f.write(ALL_KNOWN_PARAMS["TORSIONS"][missing_torsion])
            else:
                # deal with unknown improper torsion here.
                pass
        f.write("\n")
        f.write("NONBON\n")
        for missing_nonbonded in missing_params["NONBONDED"]:
            if missing_nonbonded in ALL_KNOWN_PARAMS["NONBONDED"].keys():
                f.write(ALL_KNOWN_PARAMS["NONBONDED"][missing_nonbonded])
            else:
                f.write(NewNonBonded(missing_nonbonded))
        f.write("\n\n")
    return