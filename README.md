## Benchmarking Strategies for Parallel 2D Ising Model Simulation     
This project implements the Ising model simulation using C++ and offers parallel computation using MPI (for distributed memory systems), OpenMP (for multicore CPUs), and CUDA (for NVIDIA GPUs). It models the magnetic dipole moments of atomic "spins" on a lattice, simulating phase transitions in ferromagnetic materials.

## Features

- 2D Ising Model simulation with periodic boundary conditions.
- Serial computation with C++ with two differend algorithms, Metropolis-Hastings and Wolff.
- (Hybrid) parallel computation utilizing MPI, OpenMP, and CUDA.
- Configurable parameters such as lattice size, temperature range, and number of iterations.
- Analysis of physical properties like energy and magnetization per site.

## Prerequisites

- C++ Compiler with C++17 or later support.
- OpenMP for multicore CPU parallelization.
- MPI implementation (e.g., MPICH or OpenMPI) for distributed system parallelization.
- CUDA Toolkit and NVIDIA GPU for GPU-accelerated execution. The Code has been developed and tested with CUDA compute capability higher than 6.1 and CUDA version higher than 11.

## Building the Simulation

**For CPU-based Parallelization with OpenMP:**
```bash
g++ -std=c++17 -fopenmp ising_simulation.cpp -o ising_simulation
```

**For GPU-accelerated Parallelization with CUDA:**
```bash
nvcc -std=c++17 ising_simulation.cu -o ising_simulation
```

**For Distributed System Parallelization with MPI:**
```bash
mpicxx -std=c++17 -fopenmp ising_simulation_mpi.cpp -o ising_simulation_mpi
```
or
```bash
mpic++ -std=c++17 -fopenmp ising_simulation_mpi.cpp -o ising_simulation_mpi
```

## Running the Simulation

**For OpenMP and CUDA Versions:**
```bash
./ising_simulation
```

**For the MPI Version:**
```bash
mpirun -np <number_of_processes> ./ising_simulation_mpi
```
Replace `<number_of_processes>` with the desired number of MPI processes.



## Project Structure Explanation

The structure of this CUDA program is specifically designed for development on cloud-based platforms like Google Colab and Kaggle. These platforms were chosen for their suitability in handling CUDA programming, especially for those without access to local CUDA hardware. The streamlined code structure, avoiding traditional `.cuh` and `.cu` files and `CMake`, is more practical for these environments, focusing on ease of use and efficiency in building and executing the program.




