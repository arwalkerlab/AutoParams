#!/bin/bash

## Check for conda installation
if [ -z $(which conda) ]
then
  ## Install miniconda
  wget https://repo.anaconda.com/miniconda/Miniconda3-py38_23.5.2-0-Linux-x86_64.sh
  chmod +x Miniconda3-py38_23.5.2-0-Linux-x86_64.sh
  bash ./Miniconda3-py38_23.5.2-0-Linux-x86_64.sh -b -f -p $HOME/miniconda
  eval "$($HOME/miniconda/bin/conda shell.${SHELL%%*/} hook)"
  conda init
  conda config --set auto_activate_base false
fi

## Install conda environment if it doesn't already exist
if { conda env list | grep 'autoparams'; } >/dev/null 2>&1
then
  echo "conda environment 'autoparams' already exists."
else
  conda env create -f environment.yml
  conda install -y -c psi4 libecpint --name autoparams
fi

## Add run function to .bashrc if "-b" or "--bashrc" found in arguments.
for var in $@
do
  if [ $var == "-b" ] || [ $var == "--bashrc" ]
  then
    echo "autoparams() {" >> $HOME/.bashrc
    echo "cd $PWD/Application" >> $HOME/.bashrc
    echo "flask run --host=0.0.0.0 --port=5005" >> $HOME/.bashrc
    echo "}" >> $HOME/.bashrc
    echo "Added 'autoparams' function to your bashrc!"
  fi 
done

## Run AutoParams Server if "--run" or "-r" found in arguments.
for var in $@
do
  if [ $var == "-r" ] || [ $var == "--run" ]
  then
    export AUTOPARAMS_DIR=$PWD
    cd $AUTOPARAMS_DIR/Application/
    echo "Starting Autoparams server!  Open your web browser"
    flask run --host=0.0.0.0 --port=5005
  fi
done
