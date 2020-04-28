#!/bin/sh
docker run --cap-add=SYS_PTRACE --security-opt seccomp=unconfined --rm --name debug --volumes-from postgresdb --net host -it --entrypoint /bin/bash gcr.io/mycloud/mixedLogit-main -c "gdb main -tui"
