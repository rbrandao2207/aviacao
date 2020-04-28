#!/bin/sh
docker run --rm -it --entrypoint make -v "$PWD"/src:/app gcr.io/sainf-abear/dcmodels-main
