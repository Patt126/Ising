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

#define L 100
#define N (L*L)
#define J 1.00
#define IT 2*1e3

//It must be tested more extensiveley. Chronometers time doesn0t seemtobe accurate of the time taken by the whole simulation.

void print_lattice(std::vector <int> & lattice) {

    int i;
    for (i = 0; i < N; i++) {
        if(i%L == 0) std::cout<<std::endl;
        if (lattice[i] == -1) {
            std::cout << "o" << " ";
        } else {
            std::cout << "x" << " ";


        }
    }
    std::cout<<std::endl;
}

void create_rand_vect(std::vector<int>& rand_vect_0) {
    int i;
    for (i = 0; i < IT; i++) {
        rand_vect_0.push_back(rand() % N );
    }
}


//if needed function to correctly evaluate energy
float evaluate(std::vector<int>& lattice) {
    int sum=0;
    for (int i=0;i<N;i++){
        //energy contibute
        if (i >= L) { //NO FIRST ROW
            sum += lattice[i - L]*lattice[i]*2;
        }
        if (i%L != 0) { //NO FIRST COLUMN
            sum += lattice[i-1]*lattice[i]*2;
        }

        if(i >= L*(L - 1)) { //LAST ROW
            sum += lattice[i-L*(L-1)]*lattice[i]*2; //times 2 two count also contribute where i = 0
        }
        if((i + 1) % L  == 0) { //LAST COLUMN
            sum += lattice[i-(L-1)]*lattice[i]*2;
        }
        }
    return -J*sum;
    }

void flip(std::vector <int> & lattice, std::vector<float>& prob,float& energy, int& M, int n) {

	int sum = 0;

	if (n < L) {
		sum += lattice[n+L*(L-1)];
	}
	else {
		sum += lattice[n-L];
	}
	if (n % L == 0) {
		sum += lattice[n + (L - 1)];
	}
	else {
		sum += lattice[n - 1];
	}

	if (n >= L*(L - 1)) {
		sum += lattice[n - L*(L-1)];
	}
	else {
		sum += lattice[n + L];
	}
	if ((n+1) % L == 0) {
		sum += lattice[n - (L-1)];
	}
	else {
		sum += lattice[n + 1];
	}
    int delta = 2*sum*lattice[n];
	if (delta <= 0) {
        lattice[n] = -lattice[n];
	}
	else if (delta == 4) {
        float rnd = (rand() % 10000)/1e4;
		if (rnd < prob[0] ){
            lattice[n] = -lattice[n];
		}
        else{
            return;
        }
	}
	else if (delta==8){
        float rnd = (rand() % 10000)/1e4;
		if (rnd < prob[1]) {
            lattice[n] = -lattice[n];
		}
        else{
            return;
        }
	}

	energy += 2*delta*J;
    M += 2*lattice[n];

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
                sum += lattice[i - L]*lattice[i]*2;
            }
            if (i%L != 0) { //NO FIRST COLUMN
                sum += lattice[i-1]*lattice[i]*2;
            }

            if(i >= L*(L - 1)) { //LAST ROW
                sum += lattice[i-L*(L-1)]*lattice[i]*2; //times 2 two count also contribute where i = 0
            }
            if((i + 1) % L  == 0) { //LAST COLUMN
                sum += lattice[i-(L-1)]*lattice[i]*2;
            }

		}
    energy+= -J*sum;
	}






int write_file(std::vector<float>& energy_vec, std::vector<float>& m, std::vector<int>& t){
    std::ofstream myfile ("data.txt");
    if (myfile.is_open())
    {
        for(int i = 0; i < t.size(); i ++){
            myfile << energy_vec[i] << " "<<abs(m[i])<<' '<<t[i]<<'\n';
        }
        myfile.close();
    }
    else std::cout << "Unable to open file";
    return 0;
}


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
