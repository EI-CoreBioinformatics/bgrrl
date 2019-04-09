#!/bin/bash

version=$1

# rm dist/*whl
python setup.py bdist_wheel
pip install --prefix=/ei/software/testing/bgrrl/${version}/x86_64 -U dist/bgrrl-${version}-*.whl

cp /ei/software/testing/bin/bgrrl-0.6.1 /ei/software/testing/bin/bgrrl-${version}
