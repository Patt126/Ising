#ifndef HYPERRECTANGLE_HPP
#define HYPERRECTANGLE_HPP

#include "Shape.hpp"
#include <vector>
#include <fstream>
#include <random>
#include <string>
#include <cstdlib>
#include <mpi.h>

class HyperRectangle : public Shape {
    struct Edges {
        double x;
        double y;
    };

public:
    HyperRectangle(std::string inputFile);

#endif