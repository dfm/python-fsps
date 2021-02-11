#!/bin/sh -l
set -e

git clone https://github.com/cconroy20/fsps $SPS_HOME
cp .github/workflows/Makefile $SPS_HOME/src/Makefile
make -C $SPS_HOME/src all
