#!/bin/sh
docker run --rm --name run_gentbls --volumes-from postgresdb --net host rbrandao2207/gentbls-main mktsinfo &
