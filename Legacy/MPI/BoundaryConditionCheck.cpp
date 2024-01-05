#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include<fstream>
#include<ctime>
#include <chrono>
#include <random>
#include <mpi.h>
#include <omp.h>


#define L 30
#define N (L*L)
#define A 50//Block side lenght
#define J 1.00
#define IT 6*1e6//number of iterations
#define CHUNKSIZE 1e6


void print_lattice(std::vector<int>& lattice) {

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

int atomicflip(std::vector <int> & lattice, std::vector<float>& prob, int site) {
    int sum = 0;
    int spin;
    if (site < L) {
        spin = lattice[site+L*(L-1)];
        sum += lattice[site+L*(L-1)];
    }
    else {
        spin = lattice[site-L];
        sum += lattice[site-L];
    }
    if (site % L == 0) {
        spin = lattice[site+(L-1)];
        sum += lattice[site + (L - 1)];
    }
    else {
        spin = lattice[site-1];
        sum += lattice[site - 1];
    }

    if (site >= L*(L - 1)) {
        spin = lattice[site-L*(L-1)];
        sum += lattice[site - L*(L-1)];
    }
    else {
        spin = lattice[site+L];
        sum += lattice[site + L];
    }
    if ((site+1) % L == 0) {
        spin = lattice[site-(L-1)];
        sum += lattice[site - (L-1)];
    }
    else {
        spin = lattice[site+1];
        sum += lattice[site + 1];
    }
    int delta = 2*sum*lattice[site];
    if (delta <= 0) {
#pragma omp atomic write
        lattice[site] = -lattice[site];
    }

    else if (delta == 4) {
        float rnd = (rand() % 10000)/1e4;
        if (rnd < prob[0] ){
#pragma omp atomic write
            lattice[site] = -lattice[site];
        }
        else{
            return 0;
        }
    }
    else if (delta==8){
        float rnd = (rand() % 10000)/1e4;
        if (rnd < prob[1]) {
#pragma omp atomic write
            lattice[site] = -lattice[site];
        }
        else{
            return 0;
        }
    }
    return 2*lattice[site];

}

int flip(std::vector <int> & lattice, std::vector<float>& prob, int site) {

    int sum = 0;

    if (site < L) {
        sum += lattice[site+L*(L-1)];
    }
    else {
        sum += lattice[site-L];
    }
    if (site % L == 0) {
        sum += lattice[site + (L - 1)];
    }
    else {
        sum += lattice[site - 1];
    }

    if (site >= L*(L - 1)) {
        sum += lattice[site - L*(L-1)];
    }
    else {
        sum += lattice[site + L];
    }
    if ((site+1) % L == 0) {
        sum += lattice[site - (L-1)];
    }
    else {
        sum += lattice[site + 1];
    }
    int delta = 2*sum*lattice[site];
    if (delta <= 0) {
        lattice[site] = -lattice[site];
    }

    else if (delta == 4) {
        float rnd = (rand() % 10000)/1e4;
        if (rnd < prob[0] ){
            lattice[site] = -lattice[site];
        }
        else{
            return 0;
        }
    }
    else if (delta==8){
        float rnd = (rand() % 10000)/1e4;
        if (rnd < prob[1]) {
            lattice[site] = -lattice[site];
        }
        else{
            return 0;
        }
    }
    return 2*lattice[site];

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



float simulate(float T, std::vector<int>& lattice, std::vector<int> num_vect, const int NUMBLOCKS, std::vector<bool> boundary, int rank, int size) {
    std::vector<float> prob(2);
    prob[0] = exp(-4 * J / T);
    prob[1] = exp(-8 * J / T);

    // Calculate the number of iterations per process
    int local_iterations = IT / size;
    int start_iter = rank * local_iterations;
    int end_iter = (rank + 1) * local_iterations;

    int M_local = 0;  // Local magnetization

    for (int iter = start_iter; iter < end_iter; ++iter) {
        int n = num_vect[iter];
        if (n % A == 0 || (n + 1) % A == 0 || n % (A * L) == 0 || (n + L) % (A * L) == 0) {
            M_local += atomicflip(lattice, prob, n);
        } else {
            M_local += flip(lattice, prob, n);
        }
    }

    // Aggregate local magnetizations to get the total magnetization
    int globalM;
    MPI_Allreduce(&M_local, &globalM, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return globalM;
}




void create_rand_vect(std::vector<int> &rand_vect,const int NUMBLOCKS, std::vector<int> & tStart,std::vector<bool>& boundary) {
    int i;

    for (unsigned long int i = 0; i < IT/CHUNKSIZE; i++) {
        //first term is first site in our block
        // shift = column*NUMTHREAD + row*NUMBLOCKLINE*NUMTHREAD due to current htread work
        for (int j = 0;j<CHUNKSIZE;j++) {
            int r = rand() % A, c = rand() % A;
            rand_vect.push_back(tStart[i%NUMBLOCKS] + r * L + c);
            if(r==0 || r==(A-1)||c==0||c==(A-1)){
                boundary.push_back(true);
            }
            else{
                boundary.push_back(false);
            }
        }
    }
}


//return the number of block in a line = num  blocks in a col
const int setBlockSize(int dimSideBlock,std::vector<int>& tStart) {
    if (L % dimSideBlock == 0) {
        const int numCellRow = L / dimSideBlock; // = numero di celle per colonna
        const int NUMBLOCKS = numCellRow*numCellRow;
        for(int i=0;i<  NUMBLOCKS;i++){
            tStart.push_back(floor(i/numCellRow)*A*L + i%numCellRow*A); //thread at which each block start
        }
        return NUMBLOCKS;
    }
    else {
        std::cout << "La dimensione del reticolo non è multiplo della dimensione del blocco";
    }
    return 0;
}


// Check if a site is at the boundary
bool isBoundarySite(int site, int latticeSize, int start, int end) {
    // Check if the site is at the start or end of the section
    if (site < latticeSize + start || site >= end - latticeSize) return true;
    // Check if the site is at the left or right edge of the lattice
    if (site % latticeSize == 0 || (site + 1) % latticeSize == 0) return true;
    return false;
}
/*

// Function to check boundary conditions
bool checkBoundaryConditions(const std::vector<int>& lattice, int rank, int size) {
    int sectionSize = N / size;
    int start = rank * sectionSize;
    int end = (rank + 1) * sectionSize;

    for (int i = start; i < end; ++i) {
        if (isBoundarySite(i, L, start, end)) {
            int sum = 0;
            // Calculate sum of neighboring spins
            sum += (i >= L) ? lattice[i - L] : lattice[N - L + i];  // Top neighbor
            sum += (i < N - L) ? lattice[i + L] : lattice[i + L - N];  // Bottom neighbor
            sum += (i % L != 0) ? lattice[i - 1] : lattice[i + L - 1];  // Left neighbor
            sum += ((i + 1) % L != 0) ? lattice[i + 1] : lattice[i - L + 1];  // Right neighbor

            // Check if the spin of the site 'i' is updated correctly
            int majoritySpin = (sum >= 0) ? 1 : -1;
            if (lattice[i] != majoritySpin) {
                return false;
            }
        }
    }
    return true;
} */

void setupSimpleLattice(std::vector<int>& lattice) {
    for (int i = 0; i < N; ++i) {
        lattice[i] = (i % 2 == 0) ? 1 : -1;  // Simple checkerboard pattern
    }
}


bool checkBoundaryConditions(const std::vector<int>& lattice, int rank, int size) {
    bool result = true;
    int sectionSize = N / size;
    int start = rank * sectionSize;
    int end = (rank + 1) * sectionSize;

    for (int i = start; i < end; ++i) {
        if (isBoundarySite(i, L, start, end)) {
            int sum = 0;
            // Calculate sum of neighboring spins (with wrap-around logic)
            sum += (i >= L) ? lattice[i - L] : lattice[N - L + i];
            sum += (i < N - L) ? lattice[i + L] : lattice[i + L - N];
            sum += (i % L != 0) ? lattice[i - 1] : lattice[i + L - 1];
            sum += ((i + 1) % L != 0) ? lattice[i + 1] : lattice[i - L + 1];

            int majoritySpin = (sum >= 0) ? 1 : -1;
            if (lattice[i] != majoritySpin) {
                result = false;
                break;
            }
        }
    }
    return result;
}

// Main function to run the unit tests
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<int> lattice(N, 1);  // Initialize lattice with all spins up

    // Perform test for boundary conditions
    bool boundaryConditionTest = checkBoundaryConditions(lattice, rank, size);

    if (rank == 0) {  // Only the root process prints the result
        if (boundaryConditionTest) {
            std::cout << "Boundary condition test passed." << std::endl;
        } else {
            std::cout << "Boundary condition test failed." << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}
