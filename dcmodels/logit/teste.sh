#!/bin/sh
docker run --cap-add=SYS_PTRACE --security-opt seccomp=unconfined --rm --name teste --volumes-from postgresdb --net host -it --entrypoint /bin/bash gcr.io/sainf-abear/dcmodels-main
