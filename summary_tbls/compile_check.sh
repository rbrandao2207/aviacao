#!/bin/sh
docker run --rm -it --entrypoint make -v "$PWD"/src:/app sainfibre/stbls-main
