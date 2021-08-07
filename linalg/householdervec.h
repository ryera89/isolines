#ifndef HOUSEHOLDERVEC_H
#define HOUSEHOLDERVEC_H

#include "mathvector.h"
#include "mathmatrix.h"

/* Estructura para computar el vector de Householder que se
 * utiliza para conformar la matriz de reflexion de Householder
 * utilizada en algoritmos como la factorizacion QR, este algoritmo mapea a e1
 * osea P = I - bet*(v*vt)  Px = |x|*e1; Golub, G.H., and Van Loan, C.F. 1996,
 * Matrix Computations, 3rd ed. (Baltimore: Johns Hopkins University Press)
 * Algorithm 5.1.1*/
struct householdervec
{
    VecDoub v; //Vector de householder
    double bet; // bet = 2/v*v

    void operator()(const VecDoub &x){
        //TODO To avoid everflow invoque x = x/|x|;
        v = x;
        v[0] = 1.0;
        double sig = 0.0;

        for (size_t i = 1; i < x.size(); ++i) sig+=x[i]*x[i];

        if (sig == 0.0) bet = 0.0;
        else{
            double mu = sqrt(x[0]*x[0] + sig);

            if (x[0] <= 0.0) v[0] = x[0] - mu;
            else v[0] = -sig/(x[0] + mu);

            double aux = v[0]*v[0];
            bet = 2.0*aux/(sig + aux);
            aux = v[0];
            v/=aux;
        }
    }
    //mref es como si fuera un vector tiene dimensiones (m-j,1)
    void operator()(const MatDoub2DRef &mref){
        size_t n = mref.size();
        v.resize(n,1.0);
        double sig = 0.0;
        for (size_t j = 1; j < n; ++j){
            double xj = mref(j,0);
            v[j] = xj;
            sig += xj*xj;
        }
        //v[0] = 1.0;
        if (sig == 0.0) bet = 0.0;
        else{
            double x0 = mref(0,0);
            double mu = sqrt(x0*x0 + sig);
            if (x0 <= 0.0) v[0] = x0 - mu;
            else v[0] = -sig/(x0 + mu);

            double aux = v[0]*v[0];
            bet = 2.0*aux/(sig + aux);
            aux = v[0];
            v/=aux;
        }
    }
};


#endif // HOUSEHOLDERVEC_H
