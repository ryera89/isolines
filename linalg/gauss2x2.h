#ifndef GAUSS2X2_H
#define GAUSS2X2_H

#include "mathmatrix.h"

struct gauss2x2
{
    gauss2x2(){}

    void inverse(const MatDoub &m,MatDoub &inv_m){
         assert(m.rows() == 2 && m.cols() == 2);

         inv_m.resize(2,2);

         double det = m(0,0)*m(1,1) - m(0,1)*m(1,0);

         //TODO manejar si determinante igual a cero

         inv_m(0,0) = m(1,1)/det;
         inv_m(0,1) = -m(0,1)/det;
         inv_m(1,0) = -m(1,0)/det;
         inv_m(1,1) = m(0,0)/det;
    }
};


#endif // GAUSS2X2_H
