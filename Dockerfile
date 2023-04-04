# syntax=docker/dockerfile:1
# This Dockerfile builds the WalkerLabs Automated Parameter Builder Webserver
# QM Optimization:      Psi4
# RESP Charges:         PsiRESP
# Parameter Generation: AmberTools - antechamber + prmchk2
# 2D Structure Drawing: RDKit, OpenBabel
# Parameter Testing:    AmberTools - tleap
# Current Version:      2023.03.30

FROM continuumio/miniconda
# Installation of necessary python packages for the Psi4/AmberTools Flask Application
COPY environment.yml .
RUN conda env create -f environment.yml
RUN conda install -y -c psi4 libecpint --name psi4flask
# RUN conda create --name psi4flask python=3.7.11
# RUN conda install -y -c conda-forge rdkit==2020.09.1.0 --name psi4flask
# RUN conda install -y -c conda-forge openbabel==3.1.1 --name psi4flask
# RUN conda install -y -c conda-forge flask==2.1.3 --name psi4flask
# RUN conda install -y -c conda-forge psiresp-base=0.4.2 --name psi4flask
# RUN conda install -y -c conda-forge werkzeug==2.0.3 --name psi4flask
# RUN conda install -y -c conda-forge numpy==1.21.2 --name psi4flask
# RUN conda install -y -c psi4 psi4=1.5 --name psi4flask
# RUN conda install -y -c conda-forge ambertools --name psi4flask
# RUN conda install -y -c conda-forge parmed --name psi4flask
# RUN echo "conda activate psi4flask" >> /etc/profile
# RUN echo "conda activate psi4flask" >> /root/.bashrc
# RUN echo "alias flask=/opt/conda/envs/psi4flask/bin/flask" >> /etc/profile
# RUN echo "PATH=$PATH:/opt/conda/envs/psi4flask/" >> /etc/profile
# COPY startapp.sh .
# RUN chmod +x ./startapp.sh

# Copy Flask Application data into /app/ folder.
ARG CACHEBUST=1
COPY Application/ /app/

# Make RUN commands use the new environment:
# SHELL ["conda", "run", "-n", "psi4flask", "/bin/bash", "-c"]

#COPY testpsi4.py .
#COPY test.pdb .
# ENTRYPOINT [ "conda", "run", "-n", "psi4flask", "python", "-u", "testpsi4.py"]
# ENTRYPOINT [ "conda", "run", "-n", "psi4flask", "/bin/bash", "startapp.sh"]
# WORKDIR /app/
# ENTRYPOINT ["conda", "run", "-n", "psi4flask", "flask", "run","--host=0.0.0.0"]
# CMD [ "sh", "./startapp.sh" ]
# FROM alpine:latest
# FROM conda/miniconda3
# CMD ["ls","-lrth"]

