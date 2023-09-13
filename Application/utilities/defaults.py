### STATIC VARIABLES ONLY ###
import os
import sys

# Location variables
parent = os.path.dirname
MAIN_DIR = parent(parent(__file__))
BIN_DIR = os.path.join(parent(MAIN_DIR),"bin/")
sys.path.append(BIN_DIR)
UPLOAD_FOLDER = MAIN_DIR+'/uploads/'
STATIC_DIR = MAIN_DIR+'/static/'
DATABASE_DIR = MAIN_DIR+'/database/'
TEMPLATES_DIR = MAIN_DIR+'/templates/'
LOGFILES_DIR = TEMPLATES_DIR+'logfiles/'
os.makedirs(DATABASE_DIR,exist_ok=True)
os.makedirs(UPLOAD_FOLDER,exist_ok=True)
os.makedirs(LOGFILES_DIR,exist_ok=True)

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
               "Protein":"source leaprc.protein.ff14SB\n",
               "Carbohydrate":"source leaprc.GLYCAM_06j-1\n"}

# TeraChem Variables.
USE_TERACHEM = False