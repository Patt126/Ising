#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include<fstream>
#include<ctime>
#include <chrono>
#include <omp.h>
#include <random>



#define L 100
#define N (L*L)
#define A 25//Block side lenght
#define J 1.00
#define IT 6*1e7//number of iterations

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

int flip(std::vector <int> & lattice, std::vector<float>& prob,float& energy, int site) {

    int sum = 0;
    int ID = omp_get_thread_num();
    //std::cout<<"ID: "<<" site: "<<site<<std::endl;

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



float simulate(float T,std::vector <int> & lattice, float& energy, int& M,std::vector<int> row_vect, std::vector<int> col_vect,const int NUMTHREAD, const int NUMBLOCKLINE) {

    using namespace std;
    vector<float> prob(2);
    prob[0] = exp(-4 * J / T);
    prob[1] = exp(-8 * J / T);


    vector<int> t_axis(1);
    t_axis[0] = 0;
#pragma omp parallel num_threads (NUMTHREAD)
    {
#pragma omp single
        {
            //define the number of bloch team iteration, in each iteration a flip in each block is attempted
            // even iteration is for black blocks, odd is for white one
            for (unsigned long int i = 0; i < IT / (NUMTHREAD * NUMTHREAD); i++) {
                bool color = (i % 2 == 0);

                //this two for spawns a group of thread that must be synchronized each group iteration
                for (int row = 0; row < NUMBLOCKLINE; row += 1) {
                    int column;
                    if (color) { column = row % 2; } //If black we have that even row start with black
                    else { column = 1 - row % 2; } //If white we have that odd row start with white
                    //cout<<" Column: " <<column<<endl;
                    while (column < NUMBLOCKLINE) {
                        // (find the row)*row_lwnght + BlockColumn + columnINblock
                        //first term is first site in our block
                #pragma omp task
                        {
                            // shift = column*NUMTHREAD + row*NUMBLOCKLINE*NUMTHREAD due to current htread work
                            int iteration = i * NUMTHREAD + (row * NUMBLOCKLINE + column) *
                                                            NUMTHREAD; //count the number of total iteration
                            for (int j = 0; j < NUMTHREAD; j++) {
                                int r = row_vect[iteration + j], c = col_vect[iteration + j];
                                int flipingSite = (row * A + r) * L + column * A + c;

                                flip(lattice, prob, energy, flipingSite); //true è nero

                            }
                        }
                        //cout<<"COLOR: "<<color<<" ITERATION: "<<iteration;

#pragma omp taskwait


                    }
                }
            }
        }
    }
/*
     #pragma omp parallel for num_threads (numThread) reduction(+ : M)
     {
         for (int i = 0; i < IT-1; i+=2){
             int M_local = 0;
             M += flip(lattice, prob, energy, rand_vect[i], true); //true è nero
             M += flip(lattice, prob, energy, rand_vect[i+1], false); //true è nero

         }
     }*/


    //float m = (float)M/N;
return 0 ;
}

void create_rand_vect(std::vector<int> &row_vect,std::vector<int> &col_vect) {
    int i;
    for (i = 0; i < IT; i++) {
        row_vect.push_back(rand() % A);//ROW
        col_vect.push_back(rand() % A);//COL

    }
}

//return the number of block in a line = num  blocks in a col
const int setBlockSize(int dimSideBlock) {
    if (L % dimSideBlock == 0) {
        const int numCellRow = L / dimSideBlock; // = numero di celle per colonna
        if (numCellRow % 2 == 0) {
            return numCellRow;
        } else {
            std::cout << "Con il valore inserito non posso definire una scacchiera";
        }
    } else {
        std::cout << "La dimensione del reticolo non è multiplo della dimensione del blocco";
    }
    return 0;
}


int main() {
    using namespace std;
    unsigned seed = time(0);
    srand(seed);
    vector<int> lattice(N);
    const int NUMBLOCKLINE = setBlockSize(A);
    const int NUMTHREAD = NUMBLOCKLINE*NUMBLOCKLINE/2;
    cout << NUMTHREAD <<" "<<NUMBLOCKLINE<< endl;
    float energy = 0;
    int M = 0;
    initialize_lattice(lattice, energy, M);
    print_lattice(lattice);
    float T = 0.1;
    vector<float> results(1);
    results[0] = 1;

    vector<int> row_vect;
    vector<int> col_vect;
    create_rand_vect(row_vect,col_vect);

    auto start = std::chrono::high_resolution_clock::now();  // Start timing before simulation
    while (T <= 0.2) {

        results.push_back(simulate(T, lattice, energy, M, row_vect, col_vect, NUMTHREAD, NUMBLOCKLINE));

        print_lattice(lattice);
        T += 0.1;


        cout << abs(results.back()) << endl;
    }
    auto end = std::chrono::high_resolution_clock::now();  // End timing after simulation
    std::chrono::duration<double> elapsed = end - start;  // Calculate elapsed time
    cout << elapsed.count() << endl;


    return 0;

}

