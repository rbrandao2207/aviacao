#!/bin/sh
docker run --rm --volumes-from postgresdb --net host --name run_import rbrandao2207/gentbls-main anac &
