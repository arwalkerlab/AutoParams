import glob
from rdkit import Chem
from rdkit.Chem import Draw
from openbabel import openbabel

def PDBtoSMI(pdb):
    smi_file = pdb.replace(".pdb",".smi")
    if not glob.glob(smi_file):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "smi")
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, pdb)
        obConversion.WriteFile(mol, smi_file)
    f = open(smi_file)
    smiles = f.read()
    f.close()
    return smiles.split()[0]

def Generate2DImage(smiles,image):
    Draw.MolToFile(Chem.MolFromSmiles(smiles), image, size=(1268,720))
    return

def PDBtoChemDraw(pdb,output):
    try:
        Generate2DImage(PDBtoSMI(pdb), output)
    except:
        pass
    return
