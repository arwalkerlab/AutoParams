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
        # If no parameters are missing, just skip the frcmod bit.
        return
    
    ### Generate FRCMOD from missing parameters.
    GenerateFRCMODFile(LEAPRC_DICT[moltype],file_list["FRCMOD"],missing_params)
    return