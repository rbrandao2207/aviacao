#!/bin/sh
docker run --rm --name run_estimation --volumes-from postgresdb -v "$PWD"/results:/app/results --net host sainfibre/blp-main:latest estimation &
