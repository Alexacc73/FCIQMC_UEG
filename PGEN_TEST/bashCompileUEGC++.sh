#!/bin/bash
g++ -c planeWaves.C -o planeWaves.o -O3 -std=c++11
g++ -c UEGHamiltonian.C -o UEGHamiltonian.o -O3 -std=c++11
g++ -c excitor_test.C  -o excitor_test.o -O3 -std=c++11
g++ planeWaves.o UEGHamiltonian.o excitor_test.o -o pGenTestExec -O3 -std=c++11
exit 0
