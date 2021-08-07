#ifndef FUNCTORS_H
#define FUNCTORS_H

#include <cmath>
#include <vector>
#include "linalg/mathmatrix.h"

using std::pow;
using std::vector;

struct RosenbrockFunc{
    int fun_ev;
    RosenbrockFunc():fun_ev(0){}
    double operator()(const VecDoub &pi){
        ++fun_ev;
        return 100*pow(pi[1] - pi[0]*pi[0],2) + pow(1-pi[0],2);
    }

    void df(const VecDoub &pi, VecDoub &dxi){
        dxi[0] = -400*pi[0]*(pi[1] - pi[0]*pi[0]) - 2*(1-pi[0]);
        dxi[1] = 200*(pi[1] - pi[0]*pi[0]);
    }

    VecDoub& d1f(const VecDoub &x,VecDoub &g){
        g.resize(x.size());
        g[0] = -400*x[0]*(x[1] - x[0]*x[0]) - 2*(1-x[0]);
        g[1] = 200*(x[1] - x[0]*x[0]);
        return g;
    }

    MatDoub& d2f(const VecDoub &x,MatDoub &H){
        H.resize(x.size(),x.size());
        H(0,0) = -400*(x[1] - 3.0*x[0]*x[0]) + 2.0;
        H(0,1) = -400*x[0];
        H(1,0) = -400*x[0];
        H(1,1) = 200;

        return H;
    }
};

struct NonConvFunc{
    double operator()(const VecDoub &pi){
        return 8*pi[0] + 12*pi[1] + pi[0]*pi[0] - 2.0*pi[1]*pi[1];
    }

    void df(const VecDoub &pi, VecDoub &dxi){
        dxi[0] = 8 + 2*pi[0];
        dxi[1] = 12 - 4*pi[1];
    }
};
#endif // FUNCTORS_H
