#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <unordered_set>
#include <queue>
#include <map>
#include <numeric>
#include <omp.h>
#include <mpi.h>
#include <random>

#define L 100
#define N (L*L)
#define J 1.00
#define IT 31
#define NUM_REPLICAS 10
#define SWAP_INTERVAL 5 // Interval for attempting replica exchange

//It must be tested more extensiveley

void print_lattice(std::vector <int>& lattice) {

    int i;
    for (i = 0; i < N; i++) {
        if (i % L == 0) std::cout << std::endl;
        if (lattice[i] == -1) {
            std::cout << "o" << " ";
        }
        else {
            std::cout << "x" << " ";


        }
    }
    std::cout << std::endl;
}

void create_rand_vect(std::vector<int>& rand_vect_0) {
    int i;
    for (i = 0; i < IT; i++) {
        rand_vect_0.push_back(rand() % N);
    }
}


//if needed function to correctly evaluate energy
float evaluate(std::vector<int>& lattice) {
    int sum = 0;
    for (int i = 0; i < N; i++) {
        //energy contibute
        if (i >= L) { //NO FIRST ROW
            sum += lattice[i - L] * lattice[i] * 2;
        }
        if (i % L != 0) { //NO FIRST COLUMN
            sum += lattice[i - 1] * lattice[i] * 2;
        }

        if (i >= L * (L - 1)) { //LAST ROW
            sum += lattice[i - L * (L - 1)] * lattice[i] * 2; //times 2 two count also contribute where i = 0
        }
        if ((i + 1) % L == 0) { //LAST COLUMN
            sum += lattice[i - (L - 1)] * lattice[i] * 2;
        }
    }
    return -J * sum;
}


void initialize_lattice(std::vector<int>& lattice, float& energy, int& M) {
    int k;
    int size = N;
    int sum = 0;
    for (int i = 0; i < size; i++) {
        k = rand() % 2;
        if (k == 0) {
            lattice[i] = -1;
            M -= 1;
        }
        else {
            lattice[i] = 1;
            M += 1;
        }
        //energy contibute
        if (i >= L) { //NO FIRST ROW
            sum += lattice[i - L] * lattice[i] * 2;
        }
        if (i % L != 0) { //NO FIRST COLUMN
            sum += lattice[i - 1] * lattice[i] * 2;
        }

        if (i >= L * (L - 1)) { //LAST ROW
            sum += lattice[i - L * (L - 1)] * lattice[i] * 2; //times 2 two count also contribute where i = 0
        }
        if ((i + 1) % L == 0) { //LAST COLUMN
            sum += lattice[i - (L - 1)] * lattice[i] * 2;
        }

    }
    energy += -J * sum;
}





float calculate_magnetization_per_site(const std::vector<int>& lattice) {
    int total_spin = std::accumulate(lattice.begin(), lattice.end(), 0);
    return static_cast<float>(total_spin) / N;
}

#include <omp.h>

int find_set(int x, std::vector<int>& parent) {
    while (x != parent[x]) {
        parent[x] = parent[parent[x]];  // Path compression
        x = parent[x];
    }
    return x;
}


void union_sets(int x, int y, std::vector<int>& parent, std::vector<int>& rank) {
    x = find_set(x, parent);
    y = find_set(y, parent);
    if (x != y) {
        if (rank[x] < rank[y]) {
            std::swap(x, y);
        }
        parent[y] = x;
        if (rank[x] == rank[y]) {
            rank[x]++;
        }
    }
}

void swendsen_wang_update(std::vector<int>& lattice, float T, std::vector<int>& parent, std::vector<int>& rank) {
    float P = 1 - exp(-2 * J / T);
    unsigned int seed = 0; // Declare the variable "seed"


    // Reset parent and rank for each site
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        parent[i] = i;
        rank[i] = 0;
    }

    // Form clusters
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        int right = (i % L == L - 1) ? i - (L - 1) : i + 1;
        int down = (i + L) % N;

        if (lattice[i] == lattice[right] && static_cast<float>(rand_r(&seed)) / RAND_MAX < P) {
            union_sets(i, right, parent, rank);
        }
        if (lattice[i] == lattice[down] && static_cast<float>(rand_r(&seed)) / RAND_MAX < P) {
            union_sets(i, down, parent, rank);
        }
    }

    // Flip clusters
    std::vector<int> flip_decision(N);
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        flip_decision[i] = rand_r(&seed) % 2;
    }

    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        if (flip_decision[find_set(i, parent)]) {
            lattice[i] *= -1;
        }
    }
}


void attempt_replica_exchange(std::vector<int>& lattice, float& T, int world_rank, int world_size, MPI_Comm comm) {
    int partner;
    if (world_rank % 2 == 0) {
        partner = world_rank + 1;
    } else {
        partner = world_rank - 1;
    }

    // Boundary condition for the replicas
    if (partner < 0 || partner >= world_size) return;

    float my_energy = evaluate(lattice);
    float partner_energy;
    float my_T = T, partner_T;

    // Exchange temperatures with the partner
    MPI_Sendrecv(&my_T, 1, MPI_FLOAT, partner, 0, &partner_T, 1, MPI_FLOAT, partner, 0, comm, MPI_STATUS_IGNORE);

    // Only exchange energies and attempt swapping if temperatures are different
    if (my_T != partner_T) {
        MPI_Sendrecv(&my_energy, 1, MPI_FLOAT, partner, 1, &partner_energy, 1, MPI_FLOAT, partner, 1, comm, MPI_STATUS_IGNORE);

        // Metropolis criterion
        float delta = (1.0/my_T - 1.0/partner_T) * (partner_energy - my_energy);
        if (delta < 0 || exp(-delta) > static_cast<float>(rand()) / RAND_MAX) {
            // Accept the swap
            std::vector<int> partner_lattice(N);
            MPI_Sendrecv(lattice.data(), N, MPI_INT, partner, 2, partner_lattice.data(), N, MPI_INT, partner, 2, comm, MPI_STATUS_IGNORE);
            lattice.swap(partner_lattice);
        }
    }
}


int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    srand(static_cast<unsigned>(time(nullptr)) + world_rank);

    std::vector<int> lattice(N);
    
    // Assign a unique temperature to each replica
    float T = 0.4 * (world_rank + 1); // This will give different temperatures to different replicas

    float energy = 0.0;
    int M = 0;

    initialize_lattice(lattice, energy, M);

    std::vector<int> parent(N), rank(N, 0);
    std::iota(parent.begin(), parent.end(), 0); // Initialize parent vector

    for (int iter = 0; iter < IT; ++iter) {
        swendsen_wang_update(lattice, T, parent, rank);

        // Attempt to exchange configuration with neighboring replica every 5 iterations
        if (iter % 5 == 0 && world_size > 1) {
            attempt_replica_exchange(lattice, T, world_rank, world_size, MPI_COMM_WORLD);
        }

        // Calculate and print magnetization per site
        float magnetization = calculate_magnetization_per_site(lattice);
        // Print for the first replica or every 30 iterations
        if (iter % 30 == 0) {
            std::cout << "Replica " << world_rank << ", Iteration " << iter << ", Temperature " << T 
                      << ", Magnetization per site: " << magnetization << std::endl;
        }

        // Additional measurements or operations can be added here
    }

    MPI_Finalize();
    return 0;
}
