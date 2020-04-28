#!/bin/sh
docker run --rm --name run_idxs --volumes-from postgresdb -v "$PWD"/results:/app/results --net host sainfibre/atpi-main:latest idxs &
