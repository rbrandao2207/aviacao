#!/bin/sh
docker run --cap-add=SYS_PTRACE --security-opt seccomp=unconfined --rm --name debug-stbls --volumes-from postgresdb -v "$PWD"/results:/app/results --net host -it --entrypoint /bin/bash sainfibre/stbls-main -c "gdb main"
