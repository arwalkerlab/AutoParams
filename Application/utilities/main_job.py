from .defaults import *
from .utilities import *
from .pdbclean import *

### Optimizer ###
# from .optimize import *
# from ._OPT.Psi4 import OptimizePDB
from ._OPT.TeraChem import OptimizePDB

### RESP Fitting ###
# from .respfitting import *
# from ._RESP.Psi4 import GetRESPCharges
from ._RESP.TeraChem import GetRESPCharges

from .chemdraw import *
from .parametrize import *
from .testing import *

class MainJob():
    def __init__(self,jobid,charge,multiplicity,restype,opt_bool,db_override_bool,level_of_theory,basis_set,connections):
        ### Internal variables.
        self._job_id = jobid
        self._charge = charge
        self._restype = restype
        self._multiplicity = multiplicity
        self._resname = "UNK"
        self._canon_smiles = ""
        self._resp_charges = []
        self._level_of_theory = level_of_theory
        self._basis_set = basis_set
        self._connections = connections

        ### Folder Heirarchy
        self._job_folder = UPLOAD_FOLDER + str(self._job_id)+ "/"
        self._DB_folder = DATABASE_DIR + str(self._job_id)+ "/"
        self._params_dir = self._job_folder+"Params/"
        self.file_list = {}

        ### Job Flags
        self._override_db_check = db_override_bool
        self._opt_complete = opt_bool
        self._resp_complete = False
        self._param_complete = False
        
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
        ### Check if the file has too many atoms/molecules
        if PDBisTooBig(self.file_list["Original PDB"]):
            self.LogJobMessage("Over 500 atoms found in provided structure.  This server is designed for single molecules or non-standard polymer subunits.")
            return False
        self._canon_smiles = PDBtoChemDraw(self.file_list["Original PDB"],self.file_list["ChemDraw"])
        S.call(f'cp {self.file_list["ChemDraw"]} {STATIC_DIR}{self.file_list["static_temp_png"]}',shell=True)
        if self._override_db_check:
            self.LogJobMessage("Overriding database check.")
            return True
        RefreshDB()
        if CheckSMILESinDB(self._canon_smiles):
            self.LogJobMessage("Parameters already exist in database.")
            MaybeCopy(self.file_list['JobLog'],f"{TEMPLATES_DIR}logfiles/")
            return False
        return True

    def CheckPDBQuality(self):
        SingleResidue(self.file_list["Working PDB"])
        self.file_list["Original PDB"] = self.file_list["Original PDB"]+".ORIG"
        self._resname = GetResName(self.file_list["Working PDB"])
        try:
            tmp = parmed.load_file(self.file_list["Working PDB"])
            del tmp
            return True
        except:
            self.LogJobMessage("PDB file cleaning failed.  Unable to use PDB.")
            self.LogJobMessage("<hr>")
            MaybeCopy(self.file_list['JobLog'],f"{TEMPLATES_DIR}logfiles/")
            return False
        
    def Optimize(self):
        if self._opt_complete:
            self.LogJobMessage("PDB optimization skipped.")
            self.LogJobMessage("<hr>")
            MaybeCopy(self.file_list['JobLog'],f"{TEMPLATES_DIR}logfiles/")
        else:
            OptimizePDB(self.file_list["Working PDB"],charge=self._charge,mult=self._multiplicity,method=f"{self._level_of_theory}/{self._basis_set}")
            self.LogJobMessage("PDB optimized")
            self.LogJobMessage("<hr>")
            MaybeCopy(self.file_list['JobLog'],f"{TEMPLATES_DIR}logfiles/")
            self._opt_complete = True
        return True

    def RESPCharges(self):
        self._resp_charges = GetRESPCharges(self.file_list["Working PDB"], self._charge, self._multiplicity, self._job_folder,level_of_theory=self._level_of_theory,basis_set=self._basis_set)
        if not self._resp_charges:
            self.LogJobMessage("Unable to generate RESP charges.")
            self.LogJobMessage("<hr>")
            MaybeCopy(self.file_list['JobLog'],f"{TEMPLATES_DIR}logfiles/")
            return False
        return True

    def Parametrize(self):
        os.makedirs(self._params_dir,exist_ok=True)
        S.call(f"cp {self.file_list['Working PDB']} {self._params_dir}/param.pdb",shell=True)
        os.chdir(self._params_dir)
        GenerateParameters(self.file_list,self._resp_charges,self._restype,self._resname,self._connections)
        os.chdir(MAIN_DIR)
        if all([G(self.file_list["FRCMOD"]),G(self.file_list['MOL2'])]):
            S.call(f"rm -rf {self._params_dir}",shell=True)
            return True
        self.LogJobMessage("Parameter generation failure.")
        MaybeCopy(self.file_list['JobLog'],f"{TEMPLATES_DIR}logfiles/")
        return False
    
    def TestParams(self):
        curr_miss_params = GetMissingParams(self._restype,
                            self._resname,
                            self.file_list["MOL2"],
                            self.file_list["LeapLog"],
                            frcmod = self.file_list["FRCMOD"],
                            prmtop=self.file_list["PRMTOP"],
                            inpcrd=self.file_list["INPCRD"],
                            pdb=self.file_list["Working PDB"],
                            connections=self._connections)
        
        if not curr_miss_params:
            return True
        self.LogJobMessage("Unable to generate PRMTOP and INPCRD files.")
        MaybeCopy(self.file_list['JobLog'],f"{TEMPLATES_DIR}logfiles/")
        return False
    
    def AddResultsToDB(self):
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
        self.file_list['PRMTOP'] = self.file_list['PRMTOP'].replace(UPLOAD_FOLDER,"")
        self.file_list['INPCRD'] = self.file_list['INPCRD'].replace(UPLOAD_FOLDER,"")
        return True