#include <iostream>
#include <vector>
#include <string>
#include <random>

#ifndef SHAPE_HPP
#define SHAPE_HPP

class Shape {
    public:

       Shape()= default;

       Shape(int nDim) : n(nDim);

       virtual std::vector<double> generateVector() = 0;
       virtual void calculateNorma() = 0;
       

       int getDimensions() const { return n; }
       double getNorma () const { return norma; }
       std::vector<double> getPoint () const { return point; }
       std::vector<Edges> bounds () const { return bounds; }
       const std::string getFunction() const { return function; }

    protected:
    int n{};
    double norma{};
    std::vector<double> point {};
    std::string function {};
    std::vector<Edges> bounds;
    std::mt19937 engine;
}

#endif 
