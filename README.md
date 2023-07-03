# WalkerLabs Automated Parametrizer

[![DOI](https://zenodo.org/badge/623450123.svg)](https://zenodo.org/badge/latestdoi/623450123)

This repository contains the necessary files to build and deploy a Docker container which provides a web-based generator for forcefield parameters to be used in classical molecular dynamics.

The current (public) implementation uses entirely open source engines on the back end.

- Psi4 for geometry optimization
- PsiRESP for RESP fitting
- AmberTools for generation of `.mol2` and `.frcmod` files.

The code base is designed, however, to allow users to modularly exchange certain engines.  
For instance, the geometry optimization routine currently calls `OptimizePDB(pdbfile,charge=0,mult=1,method="scf/6-31g**")` to produce an optimized structure of the given PDB file.  Users may create their own version of this function that uses any geometry optimization package they wish, including their own, so long as the function call remains the same.

#### Geometry Optimization
As the current implementation does not have job-queuing or job-recall options yet, users may need to leave their browser open for the duration of the process.  Depending on deployment server hardware, molecule size and complexity, and method (level of theory + basis set), geometry optimizations may take a **considerable** amount of time.  As such, we strongly recommend users begin with structures that are at or near a geometric minimum, or simply elect not to check the "Optimize Structure" checkbox during job submission.

Because RESP fitting and other parameters are at least somewhat dependent on geometry, the latter option is not recommended unless the submitted structure is already well-optimized.

#### RESP Fitting
The calculation of RESP charges in the current implementation uses [PsiRESP]<https://github.com/lilyminium/psiresp>, which is (currently) freely available and open-source.

#### Mol2 Generation
The construction of the `.mol2` file is performed by `antechamber` in the AmberTools suite.  However, to reduce reliance on external tools that may change in the future, work is ongoing to create a self-contained `.mol2` generator that can achieve the same results.



#### TO DO LIST
working on adding new parameters to library if unknown
