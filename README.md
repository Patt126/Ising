
# Ising Model Simulation in C++ with OpenMP and CUDA

This project is an implementation of the Ising model simulation using C++. It models the magnetic dipole moments of atomic "spins" on a lattice that are either in one of two states (+1 or -1). The spins are arranged in a grid, and the model can simulate phase transitions for ferromagnetic materials. 

## Features

- 2D Ising Model simulation with periodic boundary conditions.
- Hybrid parallel computation utilizing both OpenMP (for multicore CPUs) and CUDA (for NVIDIA GPUs).
- Configurable simulation parameters such as lattice size, temperature range, and number of iterations.
- Observation of physical properties like total energy and magnetization as a function of temperature.
- CUDA-based random number generation and lattice initialization for high-performance execution on GPUs.

## Prerequisites

To build and run the CPU-based simulation, you need:

- C++ Compiler with C++17 or later support.
- OpenMP installed and enabled in your compiler.

To build and run the GPU-accelerated simulation, you need:

- CUDA Toolkit (version 10.0 or later is recommended).
- NVIDIA GPU with Compute Capability 6.1 or higher.
- An NVIDIA GPU driver and CUDA-capable compiler (nvcc) installed and configured.

## Building the Simulation

To compile the CPU-based version of the simulation with OpenMP support, use the following command:


```bash
g++ -std=c++17 -fopenmp ising_simulation.cpp -o ising_simulation
```

For the GPU-accelerated version with CUDA, use the nvcc compiler provided by the CUDA Toolkit:

```bash
nvcc -std=c++17  ising_simulation.cu -o ising_simulation
```

**Running the Simulation**

To run the simulation, simply execute the generated binary:

```bash
./ising_simulation
```

The program's output will display the simulation results, including the energy and magnetization of the system as it evolves over time and temperature changes.

