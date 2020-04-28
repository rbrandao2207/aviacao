#!/bin/sh
docker run --rm --volumes-from postgresdb --net host gcr.io/sainf-abear/gentbls-main csvimport
