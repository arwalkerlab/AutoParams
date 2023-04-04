#!/bin/bash

########### Variables ###########
# Name of the image (instead of the hash of numbers and letters)
IMAGENAME=flask-auto-params

# Name for the resulting container to be generated.
CONTAINERNAME=AutoParametrizer

# Port assignment from inside the docker container to the host system. 
HOST_PORT=5000
DOCKER_PORT=5005
DOCKER_PORT_MAPPING="-p $HOST_PORT:$DOCKER_PORT"

# Restart container if it experiences a failure (useful for error handling)
#DOCKER_RUN_FLAGS="--restart on-failure"
DOCKER_RUN_FLAGS=""

######## Docker Commands ########
# Build the image from the Dockerfile
docker build --tag $IMAGENAME . 

docker run -it $IMAGENAME
# Run the image in a container using the variables above.
# docker run --name $CONTAINERNAME $DOCKER_PORT_MAPPING $DOCKER_RUN_FLAGS $IMAGENAME

# Purge the container
docker rm $(docker ps -aq)

# Purge the image
docker image rm $(docker images -aq)
