#include <math.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include<fstream>



#define L 100
#define N (L*L)
#define J 1.00
#define IT 6*1e7 //number of iterations

void print_lattice(std::vector < std::vector<int> >& matrix) {

    int i, j;
    for (i = 0; i < L; i++) {
        std::cout << std::endl;
        for (j = 0; j < L; j++) {
            if (matrix[i][j] == -1) {
                std::cout << "o" << " ";
            } else {
                std::cout << "x" << " ";

            }
        }
    }
}

//if needed function to correctly evaluate energy
float evaluate(std::vector < std::vector<int> >& matrix) {
    int sum=0;
    for (int i=0;i<L;i++){
        for (int j=0; j<L;j++){

            if (i != 0) {
                sum += matrix[i - 1][j]*matrix[i][j]*2;
            }
            if (j != 0) {
                sum += matrix[i][j-1]*matrix[i][j]*2;
            }

            if(i == L - 1) {
                sum += matrix[0][j]*matrix[i][j]*2; //times 2 two count also contribute where i = 0
            }
            if(j == L - 1) {
                sum += matrix[i][0]*matrix[i][j]*2;
            }
        }
        }
    return -J*sum;
    }

void flip(std::vector < std::vector<int> >& matrix, std::vector<float>& prob,float& energy, int& M) {

	int r, c;
	r = rand() % L;
	c = rand() % L;
	int sum = 0;

	if (r == 0) {
		sum += matrix[L - 1][c];
	}
	else {
		sum += matrix[r - 1][c];
	}
	if (c == 0) {
		sum += matrix[r][L - 1];
	}
	else {
		sum += matrix[r][c - 1];
	}

	if (r == L - 1) {
		sum += matrix[0][c];
	}
	else {
		sum += matrix[r + 1][c];
	}
	if (c == L - 1) {
		sum += matrix[r][0];
	}
	else {
		sum += matrix[r][c + 1];
	}
    int delta = 2*sum*matrix[r][c];
	if (delta <= 0) {
		matrix[r][c] = -matrix[r][c];
	}
	else if (delta == 4) {
        float rnd = (rand() % 10000)/1e4;
		if (rnd < prob[0] ){
			matrix[r][c] = -matrix[r][c];
		}
        else{
            return;
        }
	}
	else if (delta==8){
        float rnd = (rand() % 10000)/1e4;
		if (rnd < prob[1]) {
			matrix[r][c] = -matrix[r][c];
		}
        else{
            return;
        }
	}

	energy += 2*delta*J;
    M += 2*matrix[r][c];

}

void initialize_lattice(std::vector < std::vector<int> >& matrix, float& energy, int& M) {
	int k;
	int size = L;
    int sum = 0;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			k = rand() % 2;
			if (k == 0) {
				matrix[i][j] = -1;
                M -= 1;
			}
			else {
				matrix[i][j] = 1;
                M += 1;
			}
            //energy contibute
            if (i != 0) {
                sum += matrix[i - 1][j]*matrix[i][j]*2;
            }
            if (j != 0) {
                sum += matrix[i][j-1]*matrix[i][j]*2;
            }

            if(i == L - 1) {
                sum += matrix[0][j]*matrix[i][j]*2; //times 2 two count also contribute where i = 0
            }
            if(j == L - 1) {
                sum += matrix[i][0]*matrix[i][j]*2;
            }
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



float simulate(float T,std::vector < std::vector<int> >& lattice, float& energy, int& M) {

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
	for (int i = 1; i < IT; i++){
		flip(lattice, prob, energy, M);
        if (i % N == 0) {
            // aggiornare l'energia e la magnetizzazione ogni N step (definizione di tempo nel metropolis come da libro)
            m.push_back((float) M / N); // magnetizzazione media per sito
            energy_vec.push_back(energy);
            t_axis.push_back(i / N);

        }

	}
    return m.back();
}


int main() {
    using namespace std;
    vector<vector<int> > lattice(L, vector<int>(L));
    float energy = 0;
    int M =0;
    initialize_lattice(lattice,energy,M);
    float T = 0.1;
    vector<float> results (1);
    results[0] = 1;
    while(T<=3.5){
        results.push_back(simulate(T,lattice,energy,M));
        T += 0.1;
        cout<<abs(results.back())<<endl;
    }



    return 0;

}

