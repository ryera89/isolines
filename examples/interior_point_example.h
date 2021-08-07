#ifndef INTERIOR_POINT_EXAMPLE_H
#define INTERIOR_POINT_EXAMPLE_H

#include <cmath>
#include "linalg/mathmatrix.h"

using std::exp;
using std::pow;
struct funcObj{
    double operator()(const VecDoub &x){
        double temp = x[0]*x[1]*x[2]*x[3]*x[4];

        return exp(temp) - 0.5*pow((pow(x[0],3) + pow(x[1],3) + 1),2);
    }
    void df(const VecDoub &x,VecDoub &g){
        g.resize(x.size());

        double temp = x[0]*x[1]*x[2]*x[3]*x[4];
        temp = exp(temp);

        double temp2 = (pow(x[0],3) + pow(x[1],3) + 1);

        double temp1 = x[1]*x[2]*x[3]*x[4];
        g[0] = temp1*temp - temp2*3*pow(x[0],2);

        temp1 = x[0]*x[2]*x[3]*x[4];
        g[1] = temp1*temp - temp2*3*pow(x[1],2);

        temp1 = x[0]*x[1]*x[3]*x[4];
        g[2] = temp1*temp;

        temp1 = x[0]*x[1]*x[2]*x[4];
        g[3] = temp1*temp;

        temp1 = x[0]*x[1]*x[2]*x[3];
        g[4] = temp1*temp;
    }
    void df2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());

        double temp = x[0]*x[1]*x[2]*x[3]*x[4];
        temp = exp(temp);

        double temp2 = (pow(x[0],3) + pow(x[1],3) + 1);
        double temp1 = x[1]*x[2]*x[3]*x[4];
        double temp3 = x[0]*x[2]*x[3]*x[4];

        H(0,0) = temp1*temp1*temp - 6*x[0]*temp2 - 9*pow(x[0],4);
        H(0,1) = temp3*temp1*temp - 9*x[0]*x[0]*x[1]*x[1];

        temp3 = x[0]*x[1]*x[3]*x[4];
        H(0,2) = temp3*temp1*temp;

        temp3 = x[0]*x[1]*x[2]*x[4];
        H(0,3) = temp3*temp1*temp;

        temp3 = x[0]*x[1]*x[2]*x[3];
        H(0,4) = temp3*temp1*temp;

        temp1 =x[0]*x[2]*x[3]*x[4];
        H(1,1) = temp1*temp1*temp - 6*x[1]*temp2 - 9*pow(x[1],4);

        temp3 = x[0]*x[1]*x[3]*x[4];
        H(1,2) = temp3*temp1*temp;

        temp3 = x[0]*x[1]*x[2]*x[4];
        H(1,3) = temp3*temp1*temp;

        temp3 = x[0]*x[1]*x[2]*x[3];
        H(1,4) = temp3*temp1*temp;

        temp1 =x[0]*x[1]*x[3]*x[4];
        H(2,2) = temp1*temp1*temp;

        temp3 = x[0]*x[1]*x[2]*x[4];
        H(2,3) = temp3*temp1*temp;

        temp3 = x[0]*x[1]*x[2]*x[3];
        H(2,4) = temp3*temp1*temp;

        temp1 =x[0]*x[1]*x[2]*x[4];
        H(3,3) = temp1*temp1*temp;

        temp3 = x[0]*x[1]*x[2]*x[3];
        H(3,4) = temp3*temp1*temp;

        temp1 =x[0]*x[1]*x[2]*x[3];
        H(4,4) = temp1*temp1*temp;

        for (size_t i = 0; i < H.rows(); ++i){
            for (size_t j = i+1; j < H.cols(); ++j){
                H(j,i) = H(i,j);
            }
        }
    }
};

//Restricciones
double eq_res1(const VecDoub &x){
    double res = 0.0;
    for (size_t i = 0;  i < x.size(); ++i){
        res+=x[i]*x[i];
    }
    res-=10;
    return res;
}
void eq_res1D1(const VecDoub &x,VecDoub &g){
    g.resize(x.size());

    g[0] = 2*x[0];
    g[1] = 2*x[1];
    g[2] = 2*x[2];
    g[3] = 2*x[3];
    g[4] = 2*x[4];
}
void eq_res1D2(const VecDoub &x,MatDoub &H){
    H.resize(x.size(),x.size());
    H.setToZero();

    for (size_t i = 0; i < x.size(); ++i){
        H(i,i) = 2;
    }

}
double eq_res2(const VecDoub &x){
    return x[1]*x[2] - 5*x[3]*x[4];
}
void eq_res2D1(const VecDoub &x,VecDoub &g){
    g.resize(x.size());

    g[0] = 0;
    g[1] = x[2];
    g[2] = x[1];
    g[3] = -5*x[4];
    g[4] = -5*x[3];
}
void eq_res2D2(const VecDoub &x,MatDoub &H){
    H.resize(x.size(),x.size());
    H.setToZero();

    H(1,2) = 1;
    H(3,4) = -5;

    H(2,1) = 1;
    H(4,3) = -5;
}
double eq_res3(const VecDoub &x){
    return pow(x[0],3) + pow(x[1],3) + 1;
}
void eq_res3D1(const VecDoub &x,VecDoub &g){
    g.resize(x.size());

    g[0] = 3*x[0]*x[0];
    g[1] = 3*x[1]*x[1];
    g[2] = 0;
    g[3] = 0;
    g[4] = 0;
}
void eq_res3D2(const VecDoub &x,MatDoub &H){
    H.resize(x.size(),x.size());
    H.setToZero();

    H(0,0) = 6*x[0];
    H(1,1) = 6*x[1];
}
#endif // INTERIOR_POINT_EXAMPLE_H
