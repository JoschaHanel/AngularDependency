#!/bin/bash

rm WomDSimulation
g++ WomDSimulation.cpp `root-config --libs --cflags` -o WomDSimulation

if [ -e WomDSimulation ]
	then
	echo "Compilation successful."
else
	echo "Compilation error. Check settings and program for errors and try again."
	exit
fi

./WomDSimulation