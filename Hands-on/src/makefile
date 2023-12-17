.PHONY: all
all: main

main: Shape.hpp HyperRectangle.cpp HyperSphere.cpp Montecarlo.cpp main.cpp
	mpicc -I../muparser-2.3.4/include main.cpp -o main -L../muparser-2.3.4 -lmuparser -lstdc++ -lm -fopenmp -Wall

.PHONY: clear
clear:
	rm -f main
