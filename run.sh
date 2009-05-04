#!/bin/bash
pushd nrn
nrngui test_ss_soma.hoc
popd
pushd py
python test_ss_soma.py
popd
python plot.py

