#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include<fstream>
#include<ctime>
#include <chrono>
#include <omp.h>
#include <random>



constexpr int L = 100  ;
constexpr int N = (L*L);
constexpr int J = 1.00;




void print_lattice(std::array<int,N>& lattice) {

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
float evaluate(std::array<int,N>& lattice) {
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

void flip(std::array <int,N> & lattice, std::array<float,2>& prob, int site, int& M, int& E ) {

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
            return ;
        }
    }
    else if (delta==8){
        float rnd = (rand() % 10000)/1e4;
        if (rnd < prob[1]) {
            lattice[site] = -lattice[site];
        }
        else{
            return ;
        }
    }
    M += 2*lattice[site];

}

void initialize_lattice(std::array<int,N>& lattice, int& energy, int& M) {
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



void simulate(std::array<float,2> prob, std::array <int,N> & lattice, std::array<int,N> randomVect, int& M ,int& E ) {

    
        for (unsigned long int i = 0; i < (N); i++) {
                //first term is first site in our block
                // shift = column*NUMTHREAD + row*NUMBLOCKLINE*NUMTHREAD due to current htread work
                //int r = num_vect[i], c = block_vect[i];
                //int flipingSite = tStart[omp_get_thread_num()] + r * L  + c;
                //if (boundary[i]) {
                int n = randomVect[i];
                if (n != -1) {
                    flip(lattice, prob, n , M, E);
                }
            }
        

    
}



void create_rand_vect(std::array<int,N> &rand_vect) {
    //  N/NUMBLOCKS is the number of iteration for each task
    for (int j = 0;j<N;j++) {
        rand_vect[j] = rand()%N;
        }
}    




//Evaluate exact equilibrium value of per site magnetization from Onsager's formula for 2D case
float Mexact(float T){
    return std::pow((1.0 - std::pow(std::sinh(2 * J / T), -4)), 1.0 / 8.0);
    }
     

void translateMatrix( std::array<int,N>& inputMatrix) {
//#pragma omp parallel for num_threads(16) schedule(static, (int)CHUNKSIZE)
        int* localCopy = new int[N];
        std::memcpy(localCopy, inputMatrix.data(), N * sizeof(int));
        for (int i = 0; i < N; ++i) {
            // Check if the new indices are within bounds
            if ((i + L) < N) { //we are not moving last row
                if ((i + 1) % L != 0) //We are not in last column
                    inputMatrix[i + L + 1] = localCopy[i];
                else { //if we are in last column but not last row put element at the beginning of new row
                    inputMatrix[i + 1] = localCopy[i];
                }
            } else if ((i + 1) % L != 0) //We are not in last column of last row
            {
                // If out of bounds, put the element at the beginning
                inputMatrix[i + L + 1 - N] = localCopy[i]; //+1 to translate col -N to move to first row
            } else {//edge
                inputMatrix[0] = localCopy[i];
            }
        }
        delete[] localCopy;

    }


int main() {
    unsigned seed = time(0);
    srand(seed);
    std::array<int,N> lattice;
 

    int E = 0;
    int M = 0;
    float m = 0;
    float mExact = 1;
    initialize_lattice(lattice, E, M);
    print_lattice(lattice);
    float T = 0.1;
  

    float error = 0;

    const double tollerance = 0.001; // tollerance of a 0.1% not aligned spin
    int step = 0;
    auto start = std::chrono::high_resolution_clock::now();  // Start timing before simulation
    std::array<float,2> prob;
    std::array<int,N> random;

            
    while (T < 0.2) {
        
        prob[0] = exp(-4 * J / T);
        prob[1] = exp(-8 * J / T);
        step = 0;
        mExact = Mexact(T);
        m = static_cast<float>(M) / N;
        error = abs(abs(m)-mExact);
        while(error>tollerance ){
            m = static_cast<float>(M) / N;
            error = abs(abs(m)-mExact);
            create_rand_vect(random);
            simulate(prob,lattice, random, M, E);
            step++;
            }  
        
        T+=0.1;
        }

    auto end = std::chrono::high_resolution_clock::now();  // End timing after simulation
    std::chrono::duration<double> elapsed = end - start;  // Calculate elapsed time
    print_lattice(lattice);
    std::cout << elapsed.count() << std::endl;
    std::cout<<step*N;



    return 0;

}

