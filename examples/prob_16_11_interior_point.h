#ifndef PROB_16_11_INTERIOR_POINT_H
#define PROB_16_11_INTERIOR_POINT_H

#include "linalg/mathmatrix.h"

struct problema_16_11{

    double operator()(const VecDoub &x){
        return x[0]*x[0] + 2.0*x[1]*x[1] - 2.0*x[0] - 6.0*x[1] - 2.0*x[0]*x[1];
    }
    void df(const VecDoub &x,VecDoub &g){
        g.resize(x.size());

        g[0] = 2.0*x[0] - 2 - 2.0*x[1];
        g[1] = 4.0*x[1] - 6 - 2.0*x[0];
    }
    void df2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H(0,0) = 2.0;
        H(0,1) = -2.0;
        H(1,0) = -2.0;
        H(1,1) = 4.0;
    }

    //Restricciones
    static double rfunc1(const VecDoub &x){
        return 1.0 - 0.5*x[0] - 0.5*x[1];
    }
    static void rfunc1d1(const VecDoub &x,VecDoub &g){
        g.resize(x.size());

        g[0] = -0.5;
        g[1] = -0.5;
    }
    static void rfunc1d2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H.setToZero();
    }
    static double rfunc2(const VecDoub &x){
        return 2.0  + x[0] - 2.0*x[1];
    }
    static void rfunc2d1(const VecDoub &x,VecDoub &g){
        g.resize(x.size());

        g[0] = 1.0;
        g[1] = -2.0;
    }
    static void rfunc2d2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H.setToZero();
    }
    static double rfunc3(const VecDoub &x){
        return x[0];
    }
    static void rfunc3d1(const VecDoub &x,VecDoub &g){
        g.resize(x.size());

        g[0] = 1.0;
        g[1] = 0;
    }
    static void rfunc3d2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H.setToZero();
    }
    static double rfunc4(const VecDoub &x){
        return x[1];
    }
    static void rfunc4d1(const VecDoub &x,VecDoub &g){
        g.resize(x.size());

        g[0] = 0.0;
        g[1] = 1.0;
    }
    static void rfunc4d2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H.setToZero();
    }
};



#endif // PROB_16_11_INTERIOR_POINT_H
