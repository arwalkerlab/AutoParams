# from .future_mol2 import *
from .mol2 import *
from .frcmod import *
from .utilities import *
from .testing import *


def GenerateParameters(file_list,respcharges,moltype,resname):
    MAKE_MOL2_FILE(file_list["Working PDB"],file_list["MOL2"],respcharges)

    ### Run tleap with just the mol2 to figure out what parameters are missing.
    missing_params = GetMissingParams(moltype,resname,file_list["MOL2"],file_list["LeapLog"])
    if not missing_params:
        return
    print("Currently missing the following parameters:")
    print(missing_params)
    ### Generate FRCMOD from mol2
    GenerateFRCMODFile(file_list["MOL2"],file_list["FRCMOD"],missing_params)

    ### Modify/Update Mol2?
    # HybridizeParams(file_list["FRCMOD"],file_list["MOL2"],moltype)

    return