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

#define L 1000
#define N (L*L)
#define J 1.00
#define IT 2*1e9

void initialize_lattice(std::vector<int>& lattice, float& energy, int& M);
void print_lattice(const std::vector<int>& lattice);

int find_set(int x, std::vector<int>& parent) {
    if (x != parent[x]) {
        parent[x] = find_set(parent[x], parent);
    }
    return parent[x];
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

void swendsen_wang_update(std::vector<int>& lattice, float T) {
    float P = 1 - exp(-2 * J / T);
    std::vector<int> parent(N), rank(N, 0);
    for (int i = 0; i < N; ++i) {
        parent[i] = i;
    }

    for (int i = 0; i < N; ++i) {
        int right = i + 1;
        if (i % L == L - 1) right -= L;
        int down = (i + L) % N;

        if (lattice[i] == lattice[right] && ((float) rand() / RAND_MAX) < P) {
            union_sets(i, right, parent, rank);
        }
        if (lattice[i] == lattice[down] && ((float) rand() / RAND_MAX) < P) {
            union_sets(i, down, parent, rank);
        }
    }

    for (int i = 0; i < N; ++i) {
        if (rand() % 2 == 0) {
            lattice[find_set(i, parent)] *= -1;
        }
    }
}

float simulate(float T, std::vector<int>& lattice, float& energy, int& M, std::vector<int>& rand_vect) {
    using namespace std;
    vector<float> energy_vec(1, energy);
    vector<float> m(1, (float)M / N);
    vector<int> t_axis(1, 0);

    for (int i = 1; i < IT; ++i) {
        swendsen_wang_update(lattice, T);
        if (i % N == 0) {
            m.push_back((float)M / N);
            energy_vec.push_back(energy);
            t_axis.push_back(i / N);
        }
    }
    return m.back();
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
    print_lattice(lattice);
    auto start = chrono::high_resolution_clock::now();
    while (T <= 0.2) {
        results.push_back(simulate(T, lattice, energy, M, rand_vect));
        T += 0.1;
        print_lattice(lattice);
        cout << abs(results.back()) << endl;
    }
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "Time taken for simulation: " << elapsed.count() << " seconds" << endl;
    return 0;
}
