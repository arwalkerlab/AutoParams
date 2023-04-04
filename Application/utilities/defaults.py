MAIN_DIR = '/app/'
UPLOAD_FOLDER = '/app/uploads/'
STATIC_DIR = '/app/static/'
DATABASE_DIR = '/app/database/'
TEMPLATES_DIR = '/app/templates/'
ALLOWED_EXTENSIONS = {'txt', 'pdf', 'png', 'jpg', 'jpeg', 'pdb'}
CURRENT_JOBS = {}

SMILES_DB = DATABASE_DIR+"smilesdb.txt"
DATABASE_MOLS_SEEN = []

AVAILABLE_PSI4_MEMORY="2048 MB"