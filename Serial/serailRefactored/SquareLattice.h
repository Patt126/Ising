// SquareLattice.h

#ifndef SQUARELATTICE_H
#define SQUARELATTICE_H

#include "AbstractLattice.h"
#include <iostream>
#include <array>

template <int N>
class SquareLattice : public AbstractLattice<N> {
public:
    SquareLattice(float interactionStrength, int latticeSize);

    void printLattice() const override;
    float evaluateEnergy() const override;
    void initialize() override;

private:
    const float J;
    const int L;
    std::array<int, N> lattice;
};

#include "SquareLattice.tpp"  // Include template implementation

#endif // SQUARELATTICE_H
