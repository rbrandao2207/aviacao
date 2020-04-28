#!/bin/sh
docker run --cap-add=SYS_PTRACE --security-opt seccomp=unconfined --rm --name debug-blp --volumes-from postgresdb -v "$PWD"/results:/app/results --net host -it --entrypoint /bin/bash sainfibre/blp-main -c "gdb main"
