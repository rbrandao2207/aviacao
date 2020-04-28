#!/bin/sh
docker run --rm --name run_bash --volumes-from postgresdb -v "$PWD"/results:/app/results --net host -it --entrypoint /bin/bash sainfibre/atpi-main
