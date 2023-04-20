# syntax=docker/dockerfile:1
# This Dockerfile builds the WalkerLabs Automated Parameter Builder Webserver
# QM Optimization:      Psi4
# RESP Charges:         PsiRESP
# Parameter Generation: AmberTools - antechamber + prmchk2
# 2D Structure Drawing: RDKit, OpenBabel
# Parameter Testing:    AmberTools - tleap
# Current Version:      2023.04.04

FROM continuumio/miniconda

# Installation of necessary python packages for the Psi4/AmberTools Flask Application
# uncomment ARG CACHEBUST=1 if the environment.yml file has been modified.
ARG CACHEBUST=1
COPY environment.yml .
RUN conda env create -f environment.yml
RUN conda install -y -c psi4 libecpint --name psi4flask

# Copy Flask Application data into /app/ folder. Included CACHEBUST to force rebuild.
ARG CACHEBUST=1
COPY Application/ /app/
COPY PsiRESPJob.py /bin/PsiRESPJob
RUN chmod +x /bin/PsiRESPJob

# Move into /app/ where Flask Application is now located.
WORKDIR /app/

# Publish ports
EXPOSE 5005

# On container initialization, run the flask app, exposed
# ENTRYPOINT ["conda", "run", "-n", "psi4flask", "flask", "run","--host=0.0.0.0","--port=5005"]
CMD ["conda", "run", "-n", "psi4flask", "flask", "run","--host=0.0.0.0","--port=5005"]
# RUN echo "conda activate psi4flask" >> /root/.bashrc
# RUN echo "flask run --host=0.0.0.0 --port=5005" >> /root/.bashrc 
# 
# RUN echo "exit" >> /root/.bashrc
