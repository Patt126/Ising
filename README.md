# MonteCarlo23-24-Marini-Patricelli-Pagani
# Ising Model Simulation in C++

This project is an implementation of the Ising model simulation using C++. It aims to model magnetic dipole moments of atomic "spins" on a lattice that are either in one of two states (+1 or -1). The spins are arranged in a grid, and the model can simulate phase transitions. Parallel computing techniques with OpenMP are used to accelerate the simulation.

## Features

- 2D Ising Model simulation with periodic boundary conditions.
- Parallel computation using OpenMP to leverage multicore processors.
- Configurable simulation parameters such as lattice size, temperature range, and number of iterations.
- Observation of physical properties like total energy and magnetization as a function of temperature.

## Prerequisites

To build and run the simulation, you need:

- C++ Compiler with C++17 or later support 
- OpenMP installed and enabled in your compiler

