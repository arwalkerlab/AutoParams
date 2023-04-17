from .defaults import *
from .utilities import *
from .pdbclean import *
from .optimize import *
from .chemdraw import *
from .respfitting import *
from .parametrize import *
from .testing import *


class MainJob():
    def __init__(self,jobid):
        ### Internal variables.
        self._job_id = jobid
        self._job_folder = UPLOAD_FOLDER + str(self._job_id)+ "/"
        self._params_dir = self._job_folder+"Params/"
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
        self.LogJobMessage("<hr>")

    def LogJobMessage(self,message):
        with open(self.file_list["JobLog"],"a") as f:
            f.write(message)
            if message != "<hr>":
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
        self.file_list["LeapLog"] = os.path.join(self._job_folder, "leap.log")
        self.LogJobMessage("PDB file uploaded")
        self._canon_smiles = PDBtoChemDraw(self.file_list["Original PDB"],self.file_list["ChemDraw"])
        self.LogJobMessage("ChemDraw figure generated")
        S.call(f'cp {self.file_list["ChemDraw"]} {STATIC_DIR}{self.file_list["static_temp_png"]}',shell=True)
        self.LogJobMessage("Copying ChemDraw figure to output.")
        self.LogJobMessage("<hr>")


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
            self.LogJobMessage("<hr>")
            return True
        except:
            self.LogJobMessage("PDB file cleaning failed.  Unable to use PDB.")
            self.LogJobMessage("<hr>")
            return False
        
    def Optimize(self):
        if not self._opt_complete:
            OptimizePDB(self.file_list["Working PDB"],charge=self._charge,mult=self._multiplicity)
            self.LogJobMessage("PDB optimized")
            self.LogJobMessage("<hr>")
        else:
            self.LogJobMessage("PDB optimization skipped.")
            self.LogJobMessage("<hr>")
        return True

    def RESPCharges(self):
        self.LogJobMessage("Generating RESP Charges using PsiRESP.")
        self._resp_charges = GetRESPCharges(self.file_list["Working PDB"], self._charge, self._multiplicity, self._job_folder)
        if not self._resp_charges:
            self.LogJobMessage("Unable to generate RESP charges.")
            self.LogJobMessage("<hr>")
            return False
        self.LogJobMessage("RESP charges successfully generated.")
        self.LogJobMessage("<hr>")
        return True

    def Parametrize(self):
        self.LogJobMessage("Generating custom parameters.")
        os.makedirs(self._params_dir,exist_ok=True)
        S.call(f"cp {self.file_list['Working PDB']} {self._params_dir}/param.pdb",shell=True)
        os.chdir(self._params_dir)

        GenerateParameters(self.file_list,self._resp_charges,self._restype,self._resname)

        os.chdir(MAIN_DIR)
        if all([G(self.file_list["FRCMOD"]),G(self.file_list['MOL2'])]):
            self.LogJobMessage("Custom FRCMOD and MOL2 files generated.")
            return True
        self.LogJobMessage("Parameter generation failure.")
        return False
    
    def TestParams(self,rerun=True):
        self.LogJobMessage("Testing parameters for generation of MD inputs.")
        curr_miss_params = GetMissingParams(self._restype,
                            self._resname,
                            self.file_list["MOL2"],
                            self.file_list["LeapLog"],
                            frcmod = self.file_list["FRCMOD"],
                            prmtop=self.file_list["PRMTOP"],
                            inpcrd=self.file_list["INPCRD"],
                            pdb=self.file_list["Working PDB"])
        
        if not curr_miss_params:
            self.LogJobMessage("PRMTOP and INPCRD files successfully generated.")
            return True
        if rerun:
            GenerateParameters(self.file_list,self._resp_charges,self._restype,self._resname)
            self.TestParams(rerun=False)
        print(curr_miss_params)
        self.LogJobMessage("Unable to generate PRMTOP and INPCRD files.")
        return False
    
    def AddResultsToDB(self):
        self.LogJobMessage("In AddResultsToDB()")
        os.makedirs(self._DB_folder,exist_ok=True)
        MaybeCopy(self.file_list['FRCMOD'],f"{self._DB_folder}{self._job_id[:10]}.frcmod")
        MaybeCopy(self.file_list['MOL2'],f"{self._DB_folder}{self._job_id[:10]}.mol2")
        MaybeCopy(self.file_list['PRMTOP'],f"{self._DB_folder}{self._job_id[:10]}.prmtop")
        MaybeCopy(self.file_list['INPCRD'],f"{self._DB_folder}{self._job_id[:10]}.inpcrd")
        MaybeCopy(self.file_list['ChemDraw'],f"{self._DB_folder}{self._job_id[:10]}.png")
        MaybeCopy(self.file_list['Working PDB'],f"{self._DB_folder}{self._job_id[:10]}_cleaned.pdb")
        MaybeCopy(self.file_list['JobLog'],f"{TEMPLATES_DIR}logfiles/")
        S.call(f"echo '{self._canon_smiles}' >> {self._DB_folder}{self._job_id[:10]}.smi",shell=True)
        RefreshDB()
        self.file_list['FRCMOD'] = self.file_list['FRCMOD'].replace(UPLOAD_FOLDER,"")
        self.file_list['MOL2'] = self.file_list['MOL2'].replace(UPLOAD_FOLDER,"")
        return True