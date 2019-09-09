#!/bin/bash

set -e

docker run -it -v $PWD:/in skipping_build make -f /in/Makefile $@
