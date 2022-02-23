#!/bin/bash
g++ -g -o GEM_fitter_v3.o GEM_fitter_v3.cpp $(root-config --cflags --glibs --ldflags) -lMinuit
./GEM_fitter_v3.o
