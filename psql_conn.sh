#!/bin/bash
docker run -it --rm --net host --name psql mpostgres psql -h localhost -U postgres -d aviacao
