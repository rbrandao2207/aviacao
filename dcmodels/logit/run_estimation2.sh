#!/bin/sh
docker run --rm --name run_estimation2 -v "$PWD"/results:/app/results --net host gcr.io/sainf-abear/dcmodels-main estimation &
