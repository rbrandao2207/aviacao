#!/bin/sh
docker run --rm --name run_gentbls --volumes-from postgresdb --net host gcr.io/sainf-abear/gentbls-main &
