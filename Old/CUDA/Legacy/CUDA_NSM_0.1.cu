#include <iostream>
#include <curand.h>
#include <curand_kernel.h>

const int ITER = 5e2;
const int L = 1024;                    // Lattice size (LxL)
const int N = L * L;
const float J = 1.0f;
const float beta = 0.5f;        // Critical temperature inverse for the 2D Ising model
const int THREADS_PER_BLOCK = 16;     // Assuming the number of threads per block side
const int BLOCKS_PER_GRID = L / THREADS_PER_BLOCK;

// CUDA kernel to setup the random state
__global__ void setupStates(curandState* states, unsigned long seed) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int idy = threadIdx.y + blockIdx.y * blockDim.y;
    int id = idy * L + idx;
    if (idx < L && idy < L) {
        curand_init(seed, id, 0, &states[id]);
    }
}
__global__ void calculateEnergy(int* lattice, float* energy, float J) {
    int global_idx = blockIdx.x * blockDim.x + threadIdx.x;
    int global_idy = blockIdx.y * blockDim.y + threadIdx.y;
    int id = global_idy * L + global_idx;

    if (global_idx < L && global_idy < L) {
        int siteSpin = lattice[id];
        int rightSpin = (global_idx + 1 < L) ? lattice[id + 1] : lattice[global_idy * L];
        int downSpin = (global_idy + 1 < L) ? lattice[id + L] : lattice[global_idx];
        atomicAdd(energy, -J * siteSpin * (rightSpin + downSpin));
    }
}


// Ising model Metropolis algorithm utilizing shared memory
__global__ void metroIsing(int* lattice, curandState* states, float beta, int parity) {
    __shared__ int sharedLattice[THREADS_PER_BLOCK + 2][THREADS_PER_BLOCK + 2];

    int global_idx = blockIdx.x * blockDim.x + threadIdx.x;
    int global_idy = blockIdx.y * blockDim.y + threadIdx.y;
    int local_idx = threadIdx.x + 1;
    int local_idy = threadIdx.y + 1;
    int id = global_idy * L + global_idx;

    if (global_idx < L && global_idy < L) {
        // Load block into shared memory including halo
        sharedLattice[local_idy][local_idx] = lattice[id];

        // Load top halo
        if (threadIdx.y == 0) {
            sharedLattice[0][local_idx] = lattice[id - L];
        }

        // Load bottom halo
        if (threadIdx.y == blockDim.y - 1 || global_idy == L - 1) {
            sharedLattice[THREADS_PER_BLOCK + 1][local_idx] = lattice[id + L];
        }

        // Load left halo
        if (threadIdx.x == 0) {
            sharedLattice[local_idy][0] = lattice[id - 1];
        }

        // Load right halo
        if (threadIdx.x == blockDim.x - 1 || global_idx == L - 1) {
            sharedLattice[local_idy][THREADS_PER_BLOCK + 1] = lattice[id + 1];
        }

        __syncthreads();

        // Apply Metropolis algorithm only to the inner part if parity matches
        if ((global_idx + global_idy) % 2 == parity && global_idx < L && global_idy < L) {
            int siteSpin = sharedLattice[local_idy][local_idx];
            int spinSum = sharedLattice[local_idy + 1][local_idx] + sharedLattice[local_idy - 1][local_idx]
                + sharedLattice[local_idy][local_idx + 1] + sharedLattice[local_idy][local_idx - 1];
            int deltaE = 2 * siteSpin * spinSum;

            // Metropolis acceptance criteria
            if (deltaE <= 0 || curand_uniform(&states[id]) < expf(-beta * deltaE)) {
                lattice[id] = -siteSpin;
            }
        }
    }
}

// Main function
int main() {
    int* d_lattice;
    curandState* d_states;
    cudaMalloc(&d_lattice, N * sizeof(int));
    cudaMalloc(&d_states, N * sizeof(curandState));
    float* d_energy;
    cudaMalloc(&d_energy, sizeof(float));

    dim3 blocks(BLOCKS_PER_GRID, BLOCKS_PER_GRID);
    dim3 threads(THREADS_PER_BLOCK, THREADS_PER_BLOCK);

    // Initialization of lattice to all ones
    int* h_lattice = new int[N];
    for (int i = 0; i < N; ++i) {
        h_lattice[i] = 1;  // Set all spins to 1 (up)
    }
    cudaMemcpy(d_lattice, h_lattice, N * sizeof(int), cudaMemcpyHostToDevice);

    // Setup RNG states
    setupStates << <blocks, threads >> > (d_states, time(0));

    // Create CUDA events for timing
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // Run Monte Carlo simulation for each t from 2 to 0.2
    for (float t = 0.1f; t <= 3.0f; t += 0.1f) {
        float beta = 1.0f / t; // Compute beta as the inverse of temperature t

        cudaEventRecord(start);

        // Reset lattice to all ones for each new temperature
        //cudaMemcpy(d_lattice, h_lattice, N * sizeof(int), cudaMemcpyHostToDevice);

        for (int iter = 0; iter < ITER; ++iter) {
            metroIsing << <blocks, threads >> > (d_lattice, d_states, beta, iter % 2);
        }

        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, stop);

        // Retrieve final lattice state
        cudaMemcpy(h_lattice, d_lattice, N * sizeof(int), cudaMemcpyDeviceToHost);

        // Compute magnetization
        int mag = 0;
        for (int i = 0; i < N; ++i) {
            mag += h_lattice[i];
        }

        float energy = 0.0f;
        cudaMemcpy(d_energy, &energy, sizeof(float), cudaMemcpyHostToDevice);

        calculateEnergy << <blocks, threads >> > (d_lattice, d_energy, J);
        cudaMemcpy(&energy, d_energy, sizeof(float), cudaMemcpyDeviceToHost);

        // Average energy per site
        float energyPerSite = energy / N;

        std::cout << "Temperature: " << t << ", Final magnetization per site: " << static_cast<float>(mag) / N << ", Energy per site: " << energyPerSite << ", Time: " << milliseconds << " ms" << std::endl;
    }


    // Destroy CUDA events
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    delete[] h_lattice;
    cudaFree(d_lattice);
    cudaFree(d_states);

    return 0;
}
