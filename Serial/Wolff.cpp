#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include<fstream>
#include<ctime>
#include <chrono>



#define L 1000
#define N (L*L)
#define J 1.00
#define IT 2*1e9 //number of iterations

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



void add_to_cluster(std::vector<int>& lattice, std::unordered_set<int>& cluster, std::queue<int>& spin_queue, float P, int i) {
    if (cluster.find(i) == cluster.end()) {
        cluster.insert(i);
        spin_queue.push(i);
        std::vector<int> neighbors = {
            (i + L) % N,
            (i - L + N) % N,
            (i + 1) % L + (i / L) * L,
            (i - 1 + L) % L + (i / L) * L
        };
        for (int neighbor : neighbors) {
            if (lattice[neighbor] == lattice[i] && (rand() / (float)RAND_MAX) < P) {
                spin_queue.push(neighbor);
            }
        }
    }
}

void wolff_update(std::vector<int>& lattice, float T) {
    float P = 1 - exp(-2 * J / T);
    int seed = rand() % N;
    std::unordered_set<int> cluster;
    std::queue<int> spin_queue;
    spin_queue.push(seed);
    while (!spin_queue.empty()) {
        int current = spin_queue.front();
        spin_queue.pop();
        add_to_cluster(lattice, cluster, spin_queue, P, current);
    }
    for (int i : cluster) {
        lattice[i] *= -1;
    }
}

float simulate(float T,std::vector <int> & lattice, float& energy, int& M,std::vector<int>& rand_vect) {
    using namespace std;
    vector<float> prob(2);
    prob[0] = exp(-4 * J/T);
    prob[1] = exp(-8 * J/T);
    vector<float> energy_vec(1);
    energy_vec[0] =  energy;
    vector<float> m(1);
    m[0] = (float)M/N;
    vector<int> t_axis(1);
    t_axis[0] = 0;
    int n;
    for (int i = 1; i < IT; i++) {
        wolff_update(lattice, T);
        if (i % N == 0) {
            m.push_back((float) M / N);
            energy_vec.push_back(energy);
            t_axis.push_back(i / N);
        }
    }
    return m.back();
}

void create_rand_vect(std::vector<int>& rand_vect_0) {
    int i;
    for (i = 0; i < IT; i++) {
        rand_vect_0.push_back(rand() % N );
    }
}




int main() {
    using namespace std;
    unsigned seed = time(0);
    srand(seed);
    vector<int> lattice(N);
    float energy = 0;
    int M =0;
    initialize_lattice(lattice,energy,M);
    float T = 0.1;
    vector<float> results (1);
    results[0] = 1;
    std::vector<int> rand_vect;
    create_rand_vect(rand_vect);
    print_lattice(lattice);
    auto start = std::chrono::high_resolution_clock::now();
    while(T<=0.2){
        results.push_back(simulate(T,lattice,energy,M,rand_vect));
        T += 0.1;
        print_lattice(lattice);
        cout << abs(results.back()) << endl;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << "Time taken for simulation : "<< ": " << elapsed.count() << " seconds" << endl;
    return 0;
}


