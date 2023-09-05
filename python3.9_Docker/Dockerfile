# syntax=docker/dockerfile:1
# This Dockerfile builds the WalkerLabs Automated Parameter Builder Webserver
# QM Optimization:      Psi4
# RESP Charges:         PsiRESP
# Parameter Generation: AmberTools - antechamber + prmchk2
# 2D Structure Drawing: RDKit, OpenBabel
# Parameter Testing:    AmberTools - tleap
# Current Version:      2023.08.26

FROM continuumio/miniconda3:latest

# Installation of necessary python packages for the Psi4/AmberTools Flask Application
# uncomment ARG CACHEBUST=1 if the environment.yml file has been modified.
# ARG CACHEBUST=1
RUN conda create --name autoparams
RUN conda install -n autoparams -y python=3.9
RUN conda install -n autoparams -y -c psi4 psi4
RUN conda install -n autoparams -y -c conda-forge psiresp-base
RUN conda install -n autoparams -y -c conda-forge rdkit
RUN conda install -n autoparams -y -c conda-forge openbabel
RUN conda install -n autoparams -y -c conda-forge flask
RUN conda install -n autoparams -y -c conda-forge ambertools
RUN conda install -n autoparams -y -c conda-forge parmed

# Copy Flask Application data into /app/ folder.
ARG CACHEBUST=1
COPY Application/ /app/
COPY PsiRESPJob.py /bin/PsiRESPJob
RUN chmod +x /bin/PsiRESPJob

# Move into /app/ where Flask Application is now located.
WORKDIR /app/

# Publish ports
EXPOSE 5390

# On container initialization, run the flask app, exposed
ENTRYPOINT ["conda", "run", "-n", "autoparams", "flask", "run","--host=0.0.0.0","--port=5390"]
# CMD ["conda", "run", "-n", "autoparams", "flask", "run","--host=0.0.0.0","--port=5380"]
# RUN echo "conda activate autoparams" >> /root/.bashrc
# RUN echo "flask run --host=0.0.0.0 --port=5380" >> /root/.bashrc 

# RUN echo "exit" >> /root/.bashrc
