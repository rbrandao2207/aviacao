#!/bin/sh
docker run --rm --name run_fillmats --volumes-from postgresdb -v "$PWD"/results:/app/results --net host sainfibre/atpi-main:latest fillmats &
