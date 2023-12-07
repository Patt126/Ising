#include "Shape.hpp"
#include <fstream>
#include <cmath>
#include <random>
#include <string>
#include <vector>
#include <mpi.h>

#ifndef HYPERSPHERE_HPP
#define HYPERSPHERE_HPP

using namespace std;

class HyperSphere : public Shape {
    struct Edges {
        double x,y; 
    }; 

private:
        double radius = 0.0; 
        vector<double> center;
        int n = 0;
        double sum; 

public:
    HyperSphere(string inputFile);

    double
        getRadius();

#endif