#ifndef HYPERSPHERE_HPP
#define HYPERSPHERE_HPP

#include "Shape.hpp"
#include <fstream>
#include <cmath>
#include <random>
#include <string>
#include <vector>
#include <mpi.h>

class HyperSphere : public Shape {
public:
    HyperSphere(std::string inputFile);

    double getRadius() const { return radius; }
    std::vector<double> getCenter() const { return center; } 

private:
    double radius = 0.0;
    std::vector<double> center;
    std::vector<Edges> bounds;
    double sum;

    void calculateNorma();
};

#endif
