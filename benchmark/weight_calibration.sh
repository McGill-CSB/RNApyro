#!/bin/bash
#python RNAPyroProfile.py -d test -p profile.txt -a 0.5 -m 20 -no_profile -b 5
cd ../src
python RNAPyroProfile.py -d "../benchmark/RF00001_40.in" -a 0.5 -m 15 -no_profile -s_gc 0.5 100 -gc_sec_struct -gc_data weights.dat 
