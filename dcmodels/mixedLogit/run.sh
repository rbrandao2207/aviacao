#!/bin/sh
docker run --rm --name run -v "$PWD"/results:/app/results --net host gcr.io/mycloud/mixedLogit-main &
