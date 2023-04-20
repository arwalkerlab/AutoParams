import glob
from rdkit import Chem
from rdkit.Chem import Draw
from openbabel import openbabel
from .defaults import *
from .utilities import *

def PDBtoChemDraw(pdb,output):
    smi_file = pdb.replace(".pdb",".smi")
    if not glob.glob(smi_file):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "smi")
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, pdb)
        obConversion.WriteFile(mol, smi_file)
    smiles = Chem.MolToSmiles(Chem.MolFromSmiles(open(smi_file).read().split()[0]))
    Draw.MolToFile(Chem.MolFromSmiles(smiles), output, size=(1268,720))
    # if not CheckSMILESinDB(smiles):
    #     with open(SMILES_DB,"a") as f:
    #         f.write(smiles+"\n")
    return smiles
