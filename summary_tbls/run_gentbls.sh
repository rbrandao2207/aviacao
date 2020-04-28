#!/bin/sh
docker run --rm --name run_gentbls --volumes-from postgresdb -v "$PWD"/results:/app/results --net host sainfibre/stbls-main:latest gentbls &
