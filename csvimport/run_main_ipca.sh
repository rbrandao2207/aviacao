#!/bin/sh
docker run --rm --name run_gentbls --volumes-from postgresdb --net host sainfibre/gentbls-main ipca &
