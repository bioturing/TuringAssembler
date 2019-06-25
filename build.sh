#!/bin/bash

set -e

tar -zxvf KMC-lib.tar.gz
cd KMC
make kmc_lib

cd ..
make
