#include <cuda.h>
#include <curand_kernel.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>

#define L 100
#define N (L*L)
#define J 1.00
#define IT 5e8 // Number of iterations, should be divisible by 2 for even updates
#define NTHREADS 256 // Number of GPU threads

__device__ int get_index(int row, int col) {
    return (row * L + col) % N;
}

// The lattice uses boolean values, true for spin up (equivalent to 1) and false for spin down (equivalent to -1)
__device__ int delta_energy(bool* lattice, int r, int c) {
    int sum = lattice[get_index((r - 1 + L) % L, c)]
        + lattice[get_index((r + 1) % L, c)]
        + lattice[get_index(r, (c - 1 + L) % L)]
        + lattice[get_index(r, (c + 1) % L)];
    sum = 2 * sum - 4; // Convert sum from [0, 4] to [-4, 4] to match the original spin values
    int spin = lattice[get_index(r, c)] ? 1 : -1; // Convert bool to equivalent spin value
    return 2 * spin * sum;
}

__global__ void flip_spins(bool* lattice, float* prob, float* energy, int* M, curandState* states, bool update_black) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;

    int r = idx / L;
    int c = idx % L;
    bool is_black = ((r + c) % 2 == 0);

    if (is_black == update_black) {
        int delta = delta_energy(lattice, r, c);
        float rnd = curand_uniform(&states[idx]);

        if (delta <= 0 || (delta == 4 && rnd < prob[0]) || (delta == 8 && rnd < prob[1])) {
            lattice[get_index(r, c)] = !lattice[get_index(r, c)];
            atomicAdd(energy, delta * J);
            int spin_change = lattice[get_index(r, c)] ? 2 : -2; // Convert bool to equivalent spin change
            atomicAdd(M, spin_change);
        }
    }
}


__global__ void setup_rand_kernel(curandState* state, unsigned long seed) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        curand_init(seed, idx, 0, &state[idx]);
    }
}

__global__ void initialize_lattice_kernel(bool* lattice, float* energy, int* M, curandState* states) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        float randVal = curand_uniform(&states[idx]);
        lattice[idx] = (randVal < 0.5f) ? true : false;

        // Calculate magnetization
        atomicAdd(M, lattice[idx] ? 1 : -1);
    }
}
int main() {
    bool* dev_lattice;
    cudaMalloc((void**)&dev_lattice, N * sizeof(bool));

    curandState* dev_states;
    cudaMalloc((void**)&dev_states, N * sizeof(curandState));

    dim3 blocksPerGrid((N + NTHREADS - 1) / NTHREADS, 1, 1);
    dim3 threadsPerBlock(NTHREADS, 1, 1);

    unsigned long seed = static_cast<unsigned long>(time(nullptr));
    setup_rand_kernel << <blocksPerGrid, threadsPerBlock >> > (dev_states, seed);

    float* dev_energy;
    int* dev_M;
    cudaMalloc((void**)&dev_energy, sizeof(float));
    cudaMalloc((void**)&dev_M, sizeof(int));

    float energy = 0.0f;
    int M = 0;
    cudaMemcpy(dev_energy, &energy, sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_M, &M, sizeof(int), cudaMemcpyHostToDevice);

    initialize_lattice_kernel << <blocksPerGrid, threadsPerBlock >> > (dev_lattice, dev_energy, dev_M, dev_states);

    float* dev_probabilities;
    cudaMalloc((void**)&dev_probabilities, 2 * sizeof(float));

    for (float T = 0.2f; T <= 3.0f; T += 0.1f) {
        clock_t start_time = clock();

        float prob[2] = { exp(-4 * J / T), exp(-8 * J / T) };
        cudaMemcpy(dev_probabilities, prob, 2 * sizeof(float), cudaMemcpyHostToDevice);

        // Ensure an even number of iterations for a complete Monte Carlo sweep.
        for (unsigned long i = 0; i < IT / N; i += 2) {
            flip_spins << <blocksPerGrid, threadsPerBlock >> > (dev_lattice, dev_probabilities, dev_energy, dev_M, dev_states, true);
            flip_spins << <blocksPerGrid, threadsPerBlock >> > (dev_lattice, dev_probabilities, dev_energy, dev_M, dev_states, false);
        }

        cudaDeviceSynchronize();

        clock_t end_time = clock();

        double elapsed_secs = static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC;

        cudaMemcpy(&energy, dev_energy, sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(&M, dev_M, sizeof(int), cudaMemcpyDeviceToHost);

        std::cout << "Temperature: " << T << std::endl;
        std::cout << "Final Energy: " << energy / N << std::endl;
        std::cout << "Final Magnetization: " << static_cast<float>(M) / N << std::endl;
        std::cout << "Simulation time (seconds): " << elapsed_secs << std::endl << std::endl;
    }

    cudaFree(dev_lattice);
    cudaFree(dev_states);
    cudaFree(dev_energy);
    cudaFree(dev_M);
    cudaFree(dev_probabilities);

    return 0;
}
