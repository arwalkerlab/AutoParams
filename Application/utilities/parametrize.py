from .mol2 import *
from .frcmod import *
from .utilities import *

def GenerateParameters(file_list,respcharges,moltype):

    ### Generate the mol2 file.
    GenerateMol2File(file_list["Working PDB"],file_list["MOL2"])
    
    ### Update the mol2 file with resp charges.
    writeMol2(respcharges, file_list["MOL2"])

    ### Generate FRCMOD from mol2
    GenerateFRCMODFile(file_list["MOL2"],file_list["FRCMOD"])

    ### Modify/Update Mol2?
    HybridizeParams(file_list["FRCMOD"],file_list["MOL2"],moltype)

    return