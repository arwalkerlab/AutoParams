#!/bin/bash
########### Variables ###########
# Name of the image (instead of the hash of numbers and letters)
IMAGENAME=markahix/auto-params:psi4-python3.10
# Name for the resulting container to be generated.
CONTAINERNAME=AutoParametrizer
# Port assignment from inside the docker container to the host system. 
HOST_PORT=8088
DOCKER_PORT=5310
DOCKER_PORT_MAPPING="-p $HOST_PORT:$DOCKER_PORT"

# External Directory mounted inside Container
# Allows read/write to external database to prevent loss at container termination.
HOST_MOUNT_DIR='./JunkDir'
DOCKER_MNT_DIR=/app/uploads 
# DOCKER_MNT_DIR=/app/database
DOCKER_MOUNT_COMMAND="--mount src=$HOST_MOUNT_DIR,target=$DOCKER_MNT_DIR,type=bind"

# Restart container if it experiences a failure (useful for error handling)
#DOCKER_RUN_FLAGS="--restart on-failure"
DOCKER_RUN_FLAGS=""

######## Docker Commands ########
# Build the image from the Dockerfile
docker build --tag $IMAGENAME . 


###### Interactive Container ######
# docker run -it $DOCKER_MOUNT_COMMAND $IMAGENAME

###### Stand-alone Container ######
docker run --network=host --name $CONTAINERNAME $DOCKER_MOUNT_COMMAND $DOCKER_PORT_MAPPING $DOCKER_RUN_FLAGS $IMAGENAME

#### Cleanup on container exit ####
# Purge the container
#docker rm -f $(docker ps -aq)

# Purge the image
#docker image rm -f $(docker images -aq)
