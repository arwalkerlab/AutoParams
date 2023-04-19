from .library_management import *
from .newparams import *

# Clean and sort the library files
CleanBondsLibraryFile()
CleanAnglesLibraryFile()
CleanDihedralsLibraryFile()
CleanTorsionsLibraryFile()
# Generate the known set of parameters
ALL_KNOWN_PARAMS = GenerateParameterDictionaries()