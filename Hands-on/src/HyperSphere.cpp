#include "HyperSphere.hpp"

class HyperSphere : public Shape{
public:
    explicit HyperSphere (std::string filename){
        std::ifstream inFile;
        inFile.open(filename, std::ios_base::in);

        inFile >> n;
        inFile >> function;
        inFile >> radius;

        center.reserve(n);
        center.resize(n);
        bounds.reserve(n);
        bounds.resize(n);

        if(radius == 0.0)
            exit(-1);

        for (int i = 0; i < nDimensions; ++i){
            inFile >> center.at(i);
            bounds.at(i).x = center.at(i) - radius;
            bounds.at(i).y = center.at(i) + radius;
        }

        int r; //rank
        MPI_Comm_rank(MPI_COMM_WORLD, &r);
        engine.seed(r);

        inFile.close();

        if(r == 0) 
            calculateNorma();
    }
    std::vector<double> generateVector() override{
        std::vector<double> point;
        point.reserve(n);
        point.resize(n);

        sum = 0.0;
       
        #pragma omp parallel for reduction(+:sum)

            for (int j = 0; j < n; ++j){
                std::uniform_real_distribution<double> distribution(bounds.at(j).x, bounds.at(j).y);
                point.at(j) = distribution(engine);
                sum += (point.at(j) - center.at(j)) * (point.at(j) - center.at(j));
            }

        if(sum <= radius * radius) {
            return std::vector<double>();
        }

        return point;
    }

    void calculateNorma(){
        double volTotale = 1.;

        volTotale *= std::pow(radius, n);
        volTotale *= std::pow(M_PI, (n / 2.));
        volTotale /= std::tgamma((n / 2.) + 1);

        norma = volTotale;
    }
}
