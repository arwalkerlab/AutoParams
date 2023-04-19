### STATIC VARIABLES ONLY ###

# Location variables
MAIN_DIR = '/app/'
UPLOAD_FOLDER = '/app/uploads/'
STATIC_DIR = '/app/static/'
DATABASE_DIR = '/app/database/'
TEMPLATES_DIR = '/app/templates/'

# Extension variables
# ALLOWED_EXTENSIONS = {'txt', 'pdf', 'png', 'jpg', 'jpeg', 'pdb'}

# Session/Job variables.
CURRENT_JOBS = {}

# Database management variables.
SMILES_DB = DATABASE_DIR+"smilesdb.txt"
DATABASE_MOLS_SEEN = []

# PSI4 variables.
AVAILABLE_PSI4_MEMORY="4096 MB"


# TLEAP Variables.
LEAPRC_DICT = {"DNA":"source leaprc.DNA.OL15\n",
               "RNA":"source leaprc.RNA.OL3\n",
               "Protein":"source leaprc.protein.ff19SB\n",
               "Carbohydrate":"source leaprc.GLYCAM_06j-1\n"}