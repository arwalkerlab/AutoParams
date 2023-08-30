# WalkerLabs Automated Parametrizer

[![DOI](https://zenodo.org/badge/623450123.svg)](https://zenodo.org/badge/latestdoi/623450123)

This repository contains the necessary files to build and deploy a web-based generator for forcefield parameters to be used in classical molecular dynamics.

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

#### Polymer Connections
Many molecule types include connections to adjacent residues in a larger polymer, such as with amino acids in proteins or nucleotides in nucleic acid sequences.
The webserver defaults to "NONE" for these connections, however users may specify atom names in their submitted PDB which correspond to connections with the previous residue (Head Atom) or following residue (Tail Atom) in sequence. Please note that this requires all atoms to have unique names before submission.

## Installation Instructions
After downloading this repository to your location of choice, ensure all python packages listed in the `environment.yml` file are installed.  
Users who prefer to use `conda` installations can create a new environment from this file with 
```
conda create -f environment.yml
```
which will be named `autoparams`

## Usage - Batch Mode
To use batch mode (command line only), ensure that `bin/AutoParamsBatchMode.py` is accessible in your `$PATH`.

## Usage - WebServer Mode
To deploy the browser-based UI version for access on your local machine only:
```
cd Application/
flask run --port=8080
```

To deploy for external access, such as on a small academic network:
```
cd Application/
flask run --host=0.0.0.0 --port=8080
```

## Usage - Docker Container
For rapid deployment and ease of setup, we have also built several Docker images which are available on [DockerHub](https://hub.docker.com/repository/docker/markahix/auto-params/general).
Current builds are available for Python versions 3.8, 3.9, and 3.10.
With Docker installed, the user may retrieve and deploy the image as a container on their local machine.
```
docker run --name autoparams -p 8080:5310 markahix/auto-params:psi4-python3.10
```
In a web browser, the container may be accessed at `localhost:8080` or `127.0.0.1:8080`

#### TO DO LIST - Planned Features
- Set up conda installation method
- Additional module options for other QM packages.
- Add automatic monomer subunit capping.
