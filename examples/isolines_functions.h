#ifndef ISOLINES_FUNCTIONS_H
#define ISOLINES_FUNCTIONS_H
#include <cmath>
using namespace std;

double func1(double z,double x){
    return sqrt(z-x*x);
}
double func2(double z,double x){
    return -sqrt(z-x*x);
}
double rfunc1(double z, double x){
    return x*x + std::sqrt(z-pow(1-x,2))/10;
}
double rfunc2(double z, double x){
    return x*x - std::sqrt(z-pow(1-x,2))/10;
}
double func3p(double z, double x){
    return 1.0 + sqrt(z - pow(x-2.5,2));
}
double func3n(double z, double x){
    return 1.0 - sqrt(z - pow(x-2.5,2));
}
#endif // ISOLINES_FUNCTIONS_H
