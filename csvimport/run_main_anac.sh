#!/bin/sh
docker run --rm --volumes-from postgresdb --net host --name run_import sainfibre/gentbls-main anac &
