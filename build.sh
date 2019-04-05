#!/bin/bash

version=$1

# rm dist/*whl
python setup.py bdist_wheel
pip install --prefix=/ei/software/testing/bgrrl/${version}/x86_64 -U dist/bgrrl-${version}-*.whl
