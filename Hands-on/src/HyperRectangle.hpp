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
public:
    HyperRectangle(std::string filename);
};

#endif
