from ..utilities import *
from ..parameter_lib import ALL_KNOWN_PARAMS

def GenerateFRCMODFile(base_ff,frcmod,missing_params):
    with open(frcmod,"w") as f:
        frcmod_name = frcmod.split("/")[-1].replace(".frcmod","")
        f.write(f"{frcmod_name}, built upon {base_ff}\n")
        f.write("BOND\n")
        for missing_bond in missing_params["BONDS"]:
            if missing_bond in ALL_KNOWN_PARAMS["BONDS"].keys():
                f.write(ALL_KNOWN_PARAMS["BONDS"][missing_bond])
            else:
                # deal with unknown bond here.
                pass
        f.write("\n")
        f.write("ANGLE\n")
        for missing_angle in missing_params["ANGLES"]:
            if missing_angle in ALL_KNOWN_PARAMS["ANGLES"].keys():
                f.write(ALL_KNOWN_PARAMS["ANGLES"][missing_angle])
            else:
                # deal with unknown angle here.
                pass
        f.write("\n")
        f.write("DIHE\n")
        for missing_dihe in missing_params["DIHEDRALS"]:
            if missing_dihe in ALL_KNOWN_PARAMS["DIHEDRALS"].keys():
                f.write(ALL_KNOWN_PARAMS["DIHEDRALS"][missing_dihe])
            else:
                # deal with unknown dihedral here.
                pass
        f.write("\n")
        f.write("IMPROPER\n")
        for missing_torsion in missing_params["TORSIONS"]:
            if missing_torsion in ALL_KNOWN_PARAMS["TORSIONS"].keys():
                f.write(ALL_KNOWN_PARAMS["TORSIONS"][missing_torsion])
            else:
                # deal with unknown dihedral here.
                pass
        f.write("\n\n")