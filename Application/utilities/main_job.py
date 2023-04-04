from .defaults import *
from .utilities import *
from .pdbclean import *
from .optimize import *
from .chemdraw import PDBtoChemDraw

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
        self._job_logging = ""

        ### Initialization Functions
        os.makedirs(self._job_folder,exist_ok=True)

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
        self._job_logging += "PDB file uploaded\n"

    def CheckPDBQuality(self):
        SingleResidue(self.file_list["Original PDB"])
        try:
            tmp = parmed.load_file(self.file_list["Original PDB"])
            self.file_list["Working PDB"] = self.file_list["Original PDB"]
            self.file_list["Original PDB"] = self.file_list["Original PDB"]+".ORIG"
            self._job_logging += "PDB file cleaned\n"
            PDBtoChemDraw(self.file_list["Working PDB"],self.file_list["ChemDraw"])
            return True
        except:
            self._job_logging += "PDB file cleaning failed.  Unable to use PDB.\n"
            return False
        
    def Optimize(self):
        # make backup copy of original pdb to ORIG_{pdbfilename}
        if not self._opt_complete:
            TCOpt(self._job_folder,self.file_list["Working PDB"],self._charge)
            OptimToMainPDB(self.file_list["Working PDB"])
            self._job_logging += "PDB optimized.\n"
        else:
            print("Skipping Optimization.")

    def RESPCharges(self):
        print("In RESPCharges()")

    def Parametrize(self):
        print("In Parametrize()")

    def TestParams(self):
        print("In TestParams()")

    def AddResultsToDB(self):
        print("In AddResultsToDB()")
