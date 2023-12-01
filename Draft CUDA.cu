#include <cuda.h>
#include <curand_kernel.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>

#define L 300
#define N (L*L)
#define J 1.00
#define IT 5e8 // Number of iterations
#define NTHREADS 256 // Number of GPU threads (This is a starting point and should be tuned for the specific GPU hardware)

__device__ int get_index(int row, int col) {
    return (row * L + col) % N;
}

__device__ int delta_energy(int* lattice, int r, int c) {
    // Compute the change in energy for flipping the spin at (r, c).
    int sum = lattice[get_index((r-1+L)%L, c)] 
            + lattice[get_index((r+1)%L, c)] 
            + lattice[get_index(r, (c-1+L)%L)] 
            + lattice[get_index(r, (c+1)%L)];
    return 2 * lattice[get_index(r, c)] * sum;
}

__global__ void flip_spins(int* lattice, float* prob, float* energy, int* M, curandState* states) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;

    int r = idx / L;
    int c = idx % L;
    int delta = delta_energy(lattice, r, c);
    float rnd = curand_uniform(&states[idx]);

    if (delta <= 0 || (delta == 4 && rnd < prob[0]) || (delta == 8 && rnd < prob[1])) {
        lattice[get_index(r, c)] *= -1;
        atomicAdd(energy, 2 * delta * J);
        atomicAdd(M, 2 * lattice[get_index(r, c)]);
    }
}

__global__ void setup_rand_kernel(curandState *state, unsigned long seed) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    curand_init(seed, idx, 0, &state[idx]);
}

void initialize_lattice(int* lattice, float& energy, int& M, curandState* states) {
    // Randomly initialize lattice on GPU using CUDA and cuRAND.
    // Energy and M are updated on the host side based on the initial lattice state.
    // Implement this kernel to update lattice with random values and calculate corresponding energy and M.
}

__global__ void initialize_lattice_kernel(int* lattice, float* energy, int* M, curandState* states) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;

    // Generate a random state
    curandState localState = states[idx];
    float randVal = curand_uniform(&localState);
    states[idx] = localState; // Update the state for the next usage

    // Initialize with a random spin value
    lattice[idx] = randVal >= 0.5f ? 1 : -1;

    // Calculate contribution to energy and magnetization from the current spin
    atomicAdd(M, lattice[idx]);

    __syncthreads();

    if (threadIdx.x == 0) {
        // Only one thread calculates the energy to avoid race conditions
        float local_energy = 0.0f;
        for (int i = 0; i < N; i++) {
            int r = i / L;
            int c = i % L;
            local_energy -= J * lattice[i] * (
                lattice[get_index((r+1)%L, c)] +
                lattice[get_index(r, (c+1)%L)]
            );
        }
        atomicAdd(energy, local_energy);
    }
}
int main() {
    // Allocate memory for lattice on the GPU.
    int* dev_lattice;
    cudaMalloc((void**)&dev_lattice, N * sizeof(int));

    // Allocate memory for cuRAND states on the GPU.
    curandState* dev_states;
    cudaMalloc((void**)&dev_states, N * sizeof(curandState));

    // Setup CUDA grid and block dimensions.
    dim3 blocksPerGrid((N + NTHREADS - 1) / NTHREADS, 1, 1);
    dim3 threadsPerBlock(NTHREADS, 1, 1);

    // Initialize cuRAND states.
    setup_rand_kernel<<<blocksPerGrid, threadsPerBlock>>>(dev_states, static_cast<unsigned long>(time(nullptr)));
    cudaDeviceSynchronize(); // Synchronize to make sure cuRAND states are initialized.

    // Allocate memory on the GPU for energy and magnetization.
    float* dev_energy;
    int* dev_M;
    cudaMalloc((void**)&dev_energy, sizeof(float));
    cudaMalloc((void**)&dev_M, sizeof(int));

    // Set initial values on the host.
    float energy = 0.0f;
    int M = 0;

    // Copy initial values of energy and magnetization to the GPU.
    cudaMemcpy(dev_energy, &energy, sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_M, &M, sizeof(int), cudaMemcpyHostToDevice);

    // Initialize the lattice on the GPU.
    initialize_lattice_kernel<<<blocksPerGrid, threadsPerBlock>>>(dev_lattice, dev_energy, dev_M, dev_states);
    cudaDeviceSynchronize(); // Synchronize to ensure that initialization is complete.

    // Copy the updated energy and magnetization back to the host.
    cudaMemcpy(&energy, dev_energy, sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&M, dev_M, sizeof(int), cudaMemcpyDeviceToHost);

    // Allocate memory for probabilities on the GPU.
    float* dev_probabilities;
    cudaMalloc((void**)&dev_probabilities, 2 * sizeof(float));

    // Main simulation loop over temperature.
    float T = 2.0f; // Example temperature
    float prob[2] = {exp(-4 * J / T), exp(-8 * J / T)};
    cudaMemcpy(dev_probabilities, prob, 2 * sizeof(float), cudaMemcpyHostToDevice);

    // Simulate using the flip_spins kernel.
    flip_spins<<<blocksPerGrid, threadsPerBlock>>>(dev_lattice, dev_probabilities, dev_energy, dev_M, dev_states);
    cudaDeviceSynchronize(); // Synchronize to ensure the simulation is complete.

    // Copy the final results back to the host.
    cudaMemcpy(&energy, dev_energy, sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&M, dev_M, sizeof(int), cudaMemcpyDeviceToHost);

    // Output the results.
    std::cout << "Final Energy: " << energy << std::endl;
    std::cout << "Final Magnetization: " << ((float)M)/N << std::endl;

    // Cleanup resources.
    cudaFree(dev_lattice);
    cudaFree(dev_probabilities);
    cudaFree(dev_energy);
    cudaFree(dev_M);
    cudaFree(dev_states);

    return 0;
}
