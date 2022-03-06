#!/bin/bash
rm -r build
rm -r dist
rm -r spkmeans.egg-info
rm spkmeans.cpython-39-x86_64-linux-gnu.so

python3 setup.py build_ext --inplace
python3 setup.py build
python3 setup.py install > installation_info

