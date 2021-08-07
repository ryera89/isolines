#ifndef EJEMPLO_JUANMA_H
#define EJEMPLO_JUANMA_H

#include "linalg/mathmatrix.h"

struct ejemplo_juanma
{
    double operator()(const VecDoub &x){
        return 2.0*x[0]*x[0] + x[1]*x[1];
    }
    void df(const VecDoub &x,VecDoub &g){
        g.resize(x.size());

        g[0] = 4.0*x[0];
        g[1] = 2.0*x[1];
    }
    void df2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H(0,0) = 4.0;
        H(0,1) = 0.0;
        H(1,0) = 0.0;
        H(1,1) = 2.0;
    }

    //Restricciones
    static double rfunc1(const VecDoub &x){
        return -2.0*x[0] + x[1] - 2.0;
    }
    static void rfunc1d1(const VecDoub &x,VecDoub &g){
        g.resize(x.size());

        g[0] = -2;
        g[1] = 1;
    }
    static void rfunc1d2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H.setToZero();
    }
    static double rfunc2(const VecDoub &x){
        return -2.0*x[0] - x[1] + 4.0;
    }
    static void rfunc2d1(const VecDoub &x,VecDoub &g){
        g.resize(x.size());

        g[0] = -2;
        g[1] = -1;
    }
    static void rfunc2d2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H.setToZero();
    }
};
#endif // EJEMPLO_JUANMA_H
