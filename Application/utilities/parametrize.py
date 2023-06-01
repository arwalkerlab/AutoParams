# from .future_mol2 import *
from .mol2 import *
from .frcmod import *
from .utilities import *
from .testing import *

def GenerateParameters(file_list,respcharges,moltype,resname,connections=[]):
    MAKE_MOL2_FILE(file_list["Working PDB"],file_list["MOL2"],respcharges,connections)
    missing_params = GetMissingParams(moltype,resname,file_list["MOL2"],file_list["LeapLog"],connections=connections)
    GenerateFRCMODFile(LEAPRC_DICT[moltype],file_list["FRCMOD"],missing_params)
    orig_missing_params = missing_params
    while missing_params:
        for key,val in orig_missing_params.items():
            for xx in val:
                if xx not in missing_params[key]:
                    missing_params[key].append(xx)
        orig_missing_params = missing_params
        GenerateFRCMODFile(LEAPRC_DICT[moltype],file_list["FRCMOD"],missing_params)
        ### Run tleap with just the mol2 to figure out what parameters are missing.
        missing_params = GetMissingParams(moltype,resname,file_list["MOL2"],file_list["LeapLog"],frcmod=file_list["FRCMOD"],connections=connections)
        if orig_missing_params == missing_params:
            return
    return