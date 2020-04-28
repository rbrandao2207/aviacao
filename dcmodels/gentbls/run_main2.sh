#!/bin/sh
docker run --rm --name run_gentbls2 --volumes-from postgresdb --net host gcr.io/sainf-abear/gentbls-main &
