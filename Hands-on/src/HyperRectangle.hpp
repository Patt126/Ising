#ifndef HYPERRECTANGLE_HPP
#define HYPERRECTANGLE_HPP

#include <mpi.h>

#include <cstdlib>
#include <fstream>
#include <random>
#include <vector>
#include <string>

#include "Shape.hpp"


class HyperRectangle : public Shape {
public:
    explicit HyperRectangle(std::string filename) {
        std::ifstream inFile;
        inFile.open(filename, std::ios_base::in);

        inFile >> n;
        inFile >> function;
        bounds.reserve(n);
        bounds.resize(n);

        int i;

        for (i = 0; i < n; i++) {
            inFile >> bounds[i].x >> bounds[i].y;
            if (bounds[i].x == bounds[i].y) exit(-1);
            if (bounds[i].x > bounds[i].y)
                std::swap(bounds[i].x, bounds[i].y);
        }

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        engine.seed(rank);

        inFile.close();

        if (rank == 0) calculateNorm();
    }

    std::vector<double> generateVector() override {
        std::vector<double> point;
        point.reserve(n);
        point.resize(n);

        for (int j = 0; j < n; ++j) {
            std::uniform_real_distribution<double> distribution(bounds.[j].x, bounds[j].y);
            point.push_back(distribution(engine));
        }

        return point;
    }

private:
    std::vector<Edges> bounds;
    std::string function;

    void calculateNorm() {
        double volTotal = norm;

        for (int i = 0; i < n; i++) {
            volTotal *= std::abs(bounds[i].x - bounds[i].y);
        }
        norm = volTotal;
    }
};

#endif
