// MonteCarloSimulation.tpp

#include "MonteCarloSimulation.h"
#include <cstdlib>

namespace MyProject {

    template <std::size_t N>
    void MonteCarloSimulation<N>::simulatePhaseTransition() {
        int E = 0;
        int M = 0;
        float m = 0;
        float mExact = 1;
        float T = 0.1;

        float error = 0;
        const double tolerance = 0.001; // Tolerance of a 0.1% not aligned spin
        int step = 0;
        std::array<float, 2> prob;
        std::array<int, N> random;

        while (T < 0.2) {
            prob[0] = std::exp(-4 * 1.0 / T);
            prob[1] = std::exp(-8 * 1.0 / T);
            step = 0;
            mExact = mexact(T);
            m = static_cast<float>(M) / N;
            error = std::abs(std::abs(m) - mExact);
            while (error > tolerance) {
                m = static_cast<float>(M) / N;
                error = std::abs(std::abs(m) - mExact);
                createRandVect(random);
                simulateStep(prob, lattice.getLattice(), random, M, E);
                step++;
            }

            T += 0.1;
        }
    }

    template <std::size_t N>
    void MonteCarloSimulation<N>::createRandVect(std::array<int, N>& randVect) {
        for (int j = 0; j < N; j++) {
            randVect[j] = rand() % N;
        }
    }

    template <std::size_t N>
    float MonteCarloSimulation<N>::mexact(float T) {
        // Implement the logic to calculate Mexact for a given temperature T
        // You can use the formula provided in your original code
        return 0.0;  // Placeholder, replace with the actual calculation
    }

    template <std::size_t N>
    void MonteCarloSimulation<N>::translateMatrix(std::array<int, N>& inputMatrix) {
        std::unique_ptr<int[]> localCopy(new int[N]);
        std::memcpy(localCopy.get(), inputMatrix.data(), N * sizeof(int));

        for (int i = 0; i < N; ++i) {
            // Check if the new indices are within bounds
            if ((i + L) < N) {
                if ((i + 1) % L != 0)
                    inputMatrix[i + L + 1] = localCopy[i];
                else
                    inputMatrix[i + 1] = localCopy[i];
            } else if ((i + 1) % L != 0) {
                inputMatrix[i + L + 1 - N] = localCopy[i];
            } else {
                inputMatrix[0] = localCopy[i];
            }
        }
        // The memory managed by std::unique_ptr is automatically deallocated when it goes out of scope.
    }

    template <std::size_t N>
    void MonteCarloSimulation<N>::flip(std::array<int, N>& lattice, std::array<float, 2>& prob, int site, int& M, int& E) {
        int sum = 0;

        if (site < L) {
            sum += lattice[site + L * (L - 1)];
        } else {
            sum += lattice[site - L];
        }
        if (site % L == 0) {
            sum += lattice[site + (L - 1)];
        } else {
            sum += lattice[site - 1];
        }

        if (site >= L * (L - 1)) {
            sum += lattice[site - L * (L - 1)];
        } else {
            sum += lattice[site + L];
        }
        if ((site + 1) % L == 0) {
            sum += lattice[site - (L - 1)];
        } else {
            sum += lattice[site + 1];
        }
        int delta = 2 * sum * lattice[site];
        if (delta <= 0) {
            lattice[site] = -lattice[site];
        } else if (delta == 4) {
            float rnd = (rand() % 10000) / 1e4;
            if (rnd < prob[0]) {
                lattice[site] = -lattice[site];
            } else {
                return;
            }
        } else if (delta == 8) {
            float rnd = (rand() % 10000) / 1e4;
            if (rnd < prob[1]) {
                lattice[site] = -lattice[site];
            } else {
                return;
            }
        }
        M += 2 * lattice[site];
    }

    template <std::size_t N>
    void MonteCarloSimulation<N>::simulateStep(std::array<float, 2> prob, std::array<int, N>& lattice, std::array<int, N> randomVect, int& M, int& E) {
        for (unsigned long int i = 0; i < (N); i++) {
            int n = randomVect[i];
            if (n != -1) {
                flip(lattice, prob, n, M, E);
            }
        }
    }

} // namespace MyProject
