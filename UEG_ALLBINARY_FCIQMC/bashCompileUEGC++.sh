#!/bin/bash
g++ -c planeWaves.C -o planeWaves.o -O3 -std=c++1y
g++ -c UEGHamiltonian.C -o UEGHamiltonian.o -O3 -std=c++1y
g++ -c UEG_MAIN_binarytest.C  -o FCIQMC_MAIN.o -O3 -std=c++1y
g++ planeWaves.o UEGHamiltonian.o FCIQMC_MAIN.o -o FCIQMC_MAINexec -O3 -std=c++1y
exit 0
