#!/bin/sh
docker run --cap-add=SYS_PTRACE --security-opt seccomp=unconfined --rm --name debug-gentbls --volumes-from postgresdb --net host -it --entrypoint /bin/bash gcr.io/sainf-abear/gentbls-main -c "gdb main -tui"
