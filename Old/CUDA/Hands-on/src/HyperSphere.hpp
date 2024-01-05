#ifndef HYPERSPHERE_HPP
#define HYPERSPHERE_HPP

#include <mpi.h>

#include <cmath>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "Shape.hpp"

class HyperSphere : public Shape {
 public:
  explicit HyperSphere(std::string filename) {
    std::ifstream inFile;
    inFile.open(filename, std::ios_base::in);

    double radius;

    inFile >> n;
    inFile >> function;
    inFile >> radius;

    center.reserve(n);
    center.resize(n);
    bounds.reserve(n);
    bounds.resize(n);

    if (radius == 0.0) exit(-1);

    for (int i = 0; i < n; ++i) {
      inFile >> center[i];
      bounds[i].x = center[i] - radius;
      bounds[i].y = center[i] + radius;
    }

    int r;  // rank
    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    engine.seed(r);

    inFile.close();

    if (r == 0) calculateNorm();
  }
  std::vector<double> generateVector() override {
    std::vector<double> point;
    point.reserve(n);
    point.resize(n);

    // Creare il generatore con distribuzione in [0, 1]
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    std::vector<double> sample(n);

    // Generare il vettore con numeri casuali nel range [0, 1]
    for (int j = 0; j < n; ++j) {
        sample[j] = distribution(engine);
    }

    // Effettuare la trasformazione lineare per scalare nel range desiderato
    for (int j = 0; j < n; ++j) {
      point[j] = bounds[j].x + (bounds[j].y - bounds[j].x) * sample[j];
    }

    // Calcolare sum
    double sum = 0.0;
    #pragma omp parallel for reduction(+ : sum)
    for (int j = 0; j < n; ++j) {
      sum += (point[j] - center[j]) * (point[j] - center[j]);
    }

    if (sum <= radius * radius) {
      return std::vector<double>();
    }

    return point;
  }

 private:
  double radius = 0.0;
  std::vector<double> center;
  std::vector<Edges> bounds;
  double sum;

  void calculateNorm() {
    double volTotale = 1.;

    volTotal *= std::pow(radius, n);
    volTotal *= std::pow(M_PI, (n / 2.));
    volTotal /= std::tgamma((n / 2.) + 1);

    norm = volTotal;
  }
};

#endif
