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


#define L 2500
#define N (L*L)
#define J 1.00
#define IT 30

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

    // Reset parent and rank for each site
    std::fill(parent.begin(), parent.end(), 0);
    std::fill(rank.begin(), rank.end(), 0);
    std::iota(parent.begin(), parent.end(), 0); // Set each site as its own parent initially

    // Form clusters
    for (int i = 0; i < N; ++i) {
        int right = (i % L == L - 1) ? i - (L - 1) : i + 1;
        int down = (i + L) % N;

        if (lattice[i] == lattice[right] && static_cast<float>(rand()) / RAND_MAX < P) {
            union_sets(i, right, parent, rank);
        }
        if (lattice[i] == lattice[down] && static_cast<float>(rand()) / RAND_MAX < P) {
            union_sets(i, down, parent, rank);
        }
    }

    // Flip clusters
    std::vector<int> flip_decision(N);
    for (int i = 0; i < N; ++i) {
        flip_decision[i] = rand() % 2;
    }

    for (int i = 0; i < N; ++i) {
        if (flip_decision[find_set(i, parent)]) {
            lattice[i] *= -1;
        }
    }
}

float simulate(float T, std::vector<int>& lattice, float& energy, int& M, std::vector<int>& rand_vect) {
    using namespace std;
    vector<float> energy_vec(1, energy);
    vector<float> m(1, (float)M / N);
    vector<int> t_axis(1, 0);

    std::vector<int> parent(N), rank(N, 0);
    for (int i = 1; i < IT; ++i) {
        // Pass the additional arguments to the function
        swendsen_wang_update(lattice, T, parent, rank);
        if (i % N == 0) {
            m.push_back((float)M / N);
            energy_vec.push_back(energy);
            t_axis.push_back(i / N);
        }
    }
    return m.back();
}

float calculate_magnetization_per_site(const std::vector<int>& lattice) {
    int total_spin = std::accumulate(lattice.begin(), lattice.end(), 0);
    return static_cast<float>(total_spin) / N;
}

int main() {
    using namespace std;
    srand(static_cast<unsigned>(time(nullptr)));
    vector<int> lattice(N);
    float energy = 0;
    int M = 0;
    initialize_lattice(lattice, energy, M);
    float T = 0.1;
    vector<float> results(1, 1);
    vector<int> rand_vect;
    create_rand_vect(rand_vect);
    //print_lattice(lattice);
    auto start = chrono::high_resolution_clock::now();

    while (T <= 0.2) {
        results.push_back(simulate(T, lattice, energy, M, rand_vect));
        T += 0.1;
        //print_lattice(lattice);
        cout << "Magnetization: " << abs(results.back()) << endl;
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "Time taken for simulation: " << elapsed.count() << " seconds" << endl;
   // cout << "Maximum number of threads used: " << omp_get_max_threads() << endl;

    // Calculate and print the magnetization per site
    float magnetization_per_site = calculate_magnetization_per_site(lattice);
    cout << "Magnetization per site: " << magnetization_per_site << endl;

    return 0;
}
