#ifndef MY_PROJECT_MONTE_CARLO_SIMULATION_H
#define MY_PROJECT_MONTE_CARLO_SIMULATION_H

#include <array>
#include "AbstractMonteCarloSimulation.h"
#include "SquareLattice.h"

namespace MyProject {

    template <std::size_t N>
    class MonteCarloSimulation : public AbstractMonteCarloSimulation<N> {
    public:
        void simulatePhaseTransition();

    private:
        void createRandVect(std::array<int, N>& randVect);
        float mexact(float T);
        void translateMatrix(std::array<int, N>& inputMatrix);
        void flip(std::array<int, N>& lattice, std::array<float, 2>& prob, int site, int& M, int& E);
        void simulateStep(std::array<float, 2> prob, std::array<int, N>& lattice, std::array<int, N> randomVect, int& M, int& E);

        SquareLattice<N> lattice;  // Use SquareLattice as a private member
    };

} // namespace MyProject

#include "MonteCarloSimulation.tpp"  // Include the source file

#endif // MY_PROJECT_MONTE_CARLO_SIMULATION_H

