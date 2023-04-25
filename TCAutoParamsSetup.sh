#!/bin/bash
conda env create -f environment.yml

cat << EOF >> RunAutoParamsSite.sh
conda activate tc-autoparams-flask
source /opt/intel/oneapi/setvars.sh
module load TeraChem/1.9
cd Application/
flask run --host=0.0.0.0 --port=5005
EOF