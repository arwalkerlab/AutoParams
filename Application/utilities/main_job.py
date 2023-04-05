from .defaults import *
from .utilities import *
from .pdbclean import *
from .optimize import *
from .chemdraw import *
from .respfitting import *


class MainJob():
    def __init__(self,jobid):
        ### Internal variables.
        self._job_id = jobid
        self._job_folder = UPLOAD_FOLDER + str(self._job_id)+ "/"
        self._DB_folder = DATABASE_DIR + str(self._job_id)+ "/"
        self.file_list = {}
        self._opt_complete = False
        self._resp_complete = False
        self._param_complete = False
        self._charge = 0
        self._restype = "DNA"
        self._resname = "UNK"
        self._canon_smiles = ""
        self._multiplicity = 1
        self._resp_charges = []
        
        ### Initialization Functions
        os.makedirs(self._job_folder,exist_ok=True)
        
        self.file_list["JobLog"] = os.path.join(self._job_folder, f"{self._job_id}.html")
        self.file_list["static_temp_png"] = f"{self._job_id}.png"

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
        self.file_list["INPCRD"] = self.file_list["Original PDB"].replace(".pdb",".inpcrd")
        self.file_list["ChemDraw"] = self.file_list["Original PDB"].replace(".pdb",".png")
        self.file_list["Working PDB"] = self.file_list["Original PDB"]
        self.file_list["SMILES"] = self.file_list["Original PDB"].replace(".pdb",".png").replace(".pdb",".smi")
        self.LogJobMessage("PDB file uploaded")
        self._canon_smiles = PDBtoChemDraw(self.file_list["Original PDB"],self.file_list["ChemDraw"])
        self.LogJobMessage("ChemDraw figure generated")
        S.call(f'cp {self.file_list["ChemDraw"]} {STATIC_DIR}{self.file_list["static_temp_png"]}',shell=True)
        self.LogJobMessage("Copying ChemDraw figure to output.")


    def CheckPDBQuality(self):
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
        self._resp_charges = GetRESPCharges(self.file_list["Working PDB"], self._charge, self._multiplicity, self._job_folder)
        if not self._resp_charges:
            return False
        print(self._resp_charges)
        return True

    def Parametrize(self):
        self.LogJobMessage("In Parametrize()")
        return True
    
    def TestParams(self):
        self.LogJobMessage("In TestParams()")
        return True
    
    def AddResultsToDB(self):
        self.LogJobMessage("In AddResultsToDB()")
        os.makedirs(self._DB_folder,exist_ok=True)
        S.call(f"cp {self.file_list['FRCMOD']} {self._DB_folder}{self._job_id[:10]}.frcmod",shell=True)
        S.call(f"cp {self.file_list['MOL2']} {self._DB_folder}{self._job_id[:10]}.mol2",shell=True)
        S.call(f"cp {self.file_list['PRMTOP']} {self._DB_folder}{self._job_id[:10]}.prmtop",shell=True)
        S.call(f"cp {self.file_list['INPCRD']} {self._DB_folder}{self._job_id[:10]}.inpcrd",shell=True)
        S.call(f"cp {self.file_list['ChemDraw']} {self._DB_folder}{self._job_id[:10]}.png",shell=True)
        S.call(f"echo '{self._canon_smiles}' >> {self._DB_folder}{self._job_id[:10]}.smi",shell=True)
        S.call(f"cp {self.file_list['Working PDB']} {self._DB_folder}{self._job_id[:10]}_cleaned.pdb",shell=True)
        S.call(f"cp {self.file_list['JobLog']} {TEMPLATES_DIR}logfiles/",shell=True)
        RefreshDB()
        return True