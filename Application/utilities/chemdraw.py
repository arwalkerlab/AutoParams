import glob
from rdkit import Chem
from rdkit.Chem import Draw
from openbabel import openbabel
from .defaults import *

def PDBtoSMI(pdb):
    smi_file = pdb.replace(".pdb",".smi")
    if not glob.glob(smi_file):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "smi")
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, pdb)
        obConversion.WriteFile(mol, smi_file)
    smiles = open(smi_file).read()
    return smiles.split()[0]

def CanonicalSmiles(pdb):
    smiles = PDBtoSMI(pdb)
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))

def Generate2DImage(smiles,image):
    Draw.MolToFile(Chem.MolFromSmiles(smiles), image, size=(1268,720))
    return 

def PDBtoChemDraw(pdb,output):
    try:
        Generate2DImage(PDBtoSMI(pdb), output)
        with open(SMILES_DB,"a") as f:
            f.write(CanonicalSmiles(pdb)+"\n")
    except:
        pass
    return CanonicalSmiles(pdb)
