## Benchmarking Strategies for Parallel 2D Ising Model Simulation     
This project implements the Ising model simulation using C++ and offers parallel computation using MPI (for distributed memory systems), OpenMP (for multicore CPUs), and CUDA (for NVIDIA GPUs). It models the magnetic dipole moments of atomic "spins" on a lattice, simulating phase transitions in ferromagnetic materials.

## Features
2D Ising Model simulation with periodic boundary conditions.

- Serial computation with C++ with the Wolff algorithm.
- (Hybrid) parallel computation utilizing MPI, OpenMP, for the Simulation via the  Replica Exchange method, where each replica is simulated via parallelized Swendsen-Wang Algorithm with OPneMP, and each replica is simulated as a separate MPI process.
- Configurable parameters such as lattice size, temperature range, and number of iterations.
- Analysis of physical properties like energy and magnetization per site.

## Prerequisites

- C++ Compiler with C++17 or later support.
- OpenMP for multicore CPU parallelization.
- MPI implementation (e.g., MPICH or OpenMPI) for distributed system parallelization.
- CUDA Toolkit and NVIDIA GPU for GPU-accelerated execution. The Code has been developed and tested with CUDA compute capability higher than 6.1 and CUDA version higher than 11.

## Building the Simulation

**For GPU-accelerated Parallelization with CUDA:**
```bash
nvcc -std=c++17 metropolis_GPU.cu -o metropolis_GPU
```

## Running the Simulation

**For CUDA Version:**
```bash
./metropolis_GPU
```




## Project Structure Explanation

The structure of this CUDA program is specifically designed for development on cloud-based platforms like Google Colab and Kaggle. These platforms were chosen for their suitability in handling CUDA programming, especially for those without access to local CUDA hardware. The streamlined code structure, avoiding traditional `.cuh` and `.cu` files and `CMake`, is more practical for these environments, focusing on ease of use and efficiency in building and executing the program.




