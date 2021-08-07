#ifndef EJEM_16_04_INTERIOR_POINT_H
#define EJEM_16_04_INTERIOR_POINT_H

#include "linalg/mathmatrix.h"
#include "cmath"

using namespace std;

struct ejemplo_16_04{
    double operator()(const VecDoub &x){
        return pow(x[0]-1,2) + pow(x[1]-2.5,2);
    }
    void df(const VecDoub &x,VecDoub &g){
        g.resize(x.size());

        g[0] = 2.0*(x[0]-1);
        g[1] = 2.0*(x[1]-2.5);
    }
    void df2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H.setToZero();

        H(0,0) = 2.0;
        H(1,1) = 2.0;
    }

    //Restricciones
    static double eq_res1(const VecDoub &x){
        return x[0] - 2.0*x[1] + 2.0;
    }
    static void eq_res1D1(const VecDoub &x,VecDoub &g){
        g.resize(x.size());
        g[0] = 1;
        g[1] = -2;
    }
    static void eq_res1D2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H.setToZero();
    }
    static double eq_res2(const VecDoub &x){
        return -x[0] - 2.0*x[1] + 6.0;
    }
    static void eq_res2D1(const VecDoub &x,VecDoub &g){
        g.resize(x.size());

        g[0] = -1;
        g[1] = -2;
    }
    static void eq_res2D2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H.setToZero();
    }
    static double eq_res3(const VecDoub &x){
        return -x[0] + 2.0*x[1] + 2.0;
    }
    static void eq_res3D1(const VecDoub &x,VecDoub &g){
        g.resize(x.size());

        g[0] = -1;
        g[1] = 2;
    }
    static void eq_res3D2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H.setToZero();
    }
    static double eq_res4(const VecDoub &x){
        return x[0];
    }
    static void eq_res4D1(const VecDoub &x,VecDoub &g){
        g.resize(x.size());

        g[0] = 1;
        g[1] = 0;
    }
    static void eq_res4D2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H.setToZero();
    }
    static double eq_res5(const VecDoub &x){
        return x[1];
    }
    static void eq_res5D1(const VecDoub &x,VecDoub &g){
        g.resize(x.size());

        g[0] = 0;
        g[1] = 1;
    }
    static void eq_res5D2(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H.setToZero();
    }
};



#endif // CEJEM_16_04_INTERIOR_POINT_H
