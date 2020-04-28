#!/bin/sh
docker run --rm -it --entrypoint make -v "$PWD"/src:/app gcr.io/mycloud/mixedLogit-main
