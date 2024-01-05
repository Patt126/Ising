#ifndef ABSTRACT_MONTE_CARLO_SIMULATION_H
#define ABSTRACT_MONTE_CARLO_SIMULATION_H

#include <array>

template <std::size_t N>
class AbstractMonteCarloSimulation {
public:
    virtual void simulatePhaseTransition() = 0;
private:
    virtual void createRandVect(std::array<int, N>& randVect) = 0;
    virtual void flip(std::array<int, N>& lattice, std::array<float, 2>& prob, int site, int& M, int& E) = 0;
    virtual void simulateStep(std::array<float, 2> prob, std::array<int, N>& lattice, std::array<int, N> randomVect, int& M, int& E) = 0;
    virtual float mexact(float T){
        return std::pow((1.0 - std::pow(std::sinh(2 * J / T), -4)), 1.0 / 8.0);
    }

    virtual ~AbstractMonteCarloSimulation() = default;
};

#endif // ABSTRACT_MONTE_CARLO_SIMULATION_H
