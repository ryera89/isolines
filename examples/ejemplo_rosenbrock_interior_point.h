#ifndef EJEMPLO_ROSENBROCK_INTERIOR_POINT_H
#define EJEMPLO_ROSENBROCK_INTERIOR_POINT_H

#include "linalg/mathmatrix.h"

struct ejemplo_rosenbrock
{
    double operator()(const VecDoub &pi){
        return 100*pow(pi[1] - pi[0]*pi[0],2) + pow(1-pi[0],2);
    }
    void df(const VecDoub &x,VecDoub &g){
        g.resize(x.size());
        g[0] = -400*x[0]*(x[1] - x[0]*x[0]) - 2*(1-x[0]);
        g[1] = 200*(x[1] - x[0]*x[0]);
    }

    void df2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H(0,0) = -400*(x[1] - 3.0*x[0]*x[0]) + 2.0;
        H(0,1) = -400*x[0];
        H(1,0) = -400*x[0];
        H(1,1) = 200;
    }

    //Restricciones
    static double rfunc1(const VecDoub &x){
        return x[0] + 5;
    }
    static void rfunc1d1(const VecDoub &x,VecDoub &g){
        g.resize(x.size());

        g[0] = 1;
        g[1] = 0;
    }
    static void rfunc1d2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H.setToZero();
    }
    static double rfunc2(const VecDoub &x){
        return x[1] + 3;
    }
    static void rfunc2d1(const VecDoub &x,VecDoub &g){
        g.resize(x.size());

        g[0] = 0.0;
        g[1] = 1.0;
    }
    static void rfunc2d2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H.setToZero();
    }
    static double rfunc3(const VecDoub &x){
        return x[1] - x[0]*x[0] - 3*x[0];
    }
    static void rfunc3d1(const VecDoub &x,VecDoub &g){
        g.resize(x.size());

        g[0] = -2.0*x[0] - 3;
        g[1] = 1;
    }
    static void rfunc3d2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H.setToZero();
        H(0,0) = -2.0;
    }
};
#endif // EJEMPLO_ROSENBROCK_INTERIOR_POINT_H
