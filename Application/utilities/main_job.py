from .defaults import *
from .utilities import *
from .pdbclean import *
from .optimize import *
from .chemdraw import *

class MainJob():
    def __init__(self,jobid):
        ### Internal variables.
        self._job_id = jobid
        self._job_folder = UPLOAD_FOLDER + str(self._job_id)+ "/"
        self.file_list = {}
        self._opt_complete = False
        self._resp_complete = False
        self._param_complete = False
        self._charge = 0
        self._restype = "DNA"
        self._resname = "UNK"
        self._canon_smiles = ""
        self._multiplicity = 1
        
        ### Initialization Functions
        os.makedirs(self._job_folder,exist_ok=True)
        self.file_list["JobLog"] = os.path.join(self._job_folder, f"{self._job_id}.html")

    def LogJobMessage(self,message):
        with open(self.file_list["JobLog"],"a") as f:
            f.write(message)
            f.write("\n<br>\n")

    def UploadPDBFile(self,pdbfile):
        filename = secure_filename(pdbfile.filename)
        pdbsavepath = os.path.join(self._job_folder, filename)
        pdbfile.save(pdbsavepath)
        self.file_list["Original PDB"] = pdbsavepath
        self.file_list["FRCMOD"] = self.file_list["Original PDB"].replace(".pdb",".frcmod")
        self.file_list["MOL2"] = self.file_list["Original PDB"].replace(".pdb",".mol2")
        self.file_list["PRMTOP"] = self.file_list["Original PDB"].replace(".pdb",".prmtop")
        self.file_list["inpcrd"] = self.file_list["Original PDB"].replace(".pdb",".inpcrd")
        self.file_list["ChemDraw"] = self.file_list["Original PDB"].replace(".pdb",".png")
        self.file_list["Working PDB"] = self.file_list["Original PDB"]
        self.LogJobMessage("PDB file uploaded")

    def CheckPDBQuality(self):
        self._canon_smiles = PDBtoChemDraw(self.file_list["Original PDB"],self.file_list["ChemDraw"])
        self.LogJobMessage("ChemDraw figure generated")
        S.call(f'cp {self.file_list["ChemDraw"]} {STATIC_DIR}{self.file_list["ChemDraw"].split("/")[-1]}',shell=True)
        self.LogJobMessage("Copying ChemDraw figure to output.")
        SingleResidue(self.file_list["Working PDB"])
        self.file_list["Original PDB"] = self.file_list["Original PDB"]+".ORIG"
        self.LogJobMessage("PDB file checked for single molecule/residue.")
        self._resname = GetResName(self.file_list["Working PDB"])
        self.LogJobMessage(f"PDB residue name found: {self._resname}")
        try:
            self.LogJobMessage("Loading PDB into parmed")
            tmp = parmed.load_file(self.file_list["Working PDB"])
            self.LogJobMessage("PDB file cleaned")        
            return True
        except:
            self.LogJobMessage("PDB file cleaning failed.  Unable to use PDB.")
            return False
        
    def Optimize(self):
        if not self._opt_complete:
            OptimizePDB(self.file_list["Working PDB"],charge=self._charge,mult=self._multiplicity)
            self.LogJobMessage("PDB optimized")
        else:
            self.LogJobMessage("PDB optimization skipped.")
        return True

    def RESPCharges(self):
        self.LogJobMessage("In RESPCharges()")
        return True

    def Parametrize(self):
        self.LogJobMessage("In Parametrize()")
        return True
    
    def TestParams(self):
        self.LogJobMessage("In TestParams()")
        return True
    
    def AddResultsToDB(self):
        self.LogJobMessage("In AddResultsToDB()")
        return True