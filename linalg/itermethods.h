#ifndef ITERMETHODS_H
#define ITERMETHODS_H
#include "mathmatrix.h"

/*Algoritmo 5.2 CG adaptado para su aplicacion a sistemas reducidos en esta implementacion
 * no se usa precondicioner (ie W = I) J. Nocedal Numerical Optimization
 * problema min(x)[ 1/2*(x_t*Z_t*G*Z*x + x_t*c ] */
inline void cg_redsys(const MatDoub &Z, const MatDoub &G, const VecDoub &c,VecDoub &x, const double tol = 1.0e-8){
    assert(Z.rows() == G.cols());
    VecDoub vaux = Z*x;
    vaux = G*vaux;
    vaux = vaux*Z;

    VecDoub r(c);
    r+=vaux;

    VecDoub p(r.size());

    for(size_t i = 0; i < p.size(); ++i) p[i] = -r[i];

    double r_norm = r.norm();

    while (r_norm > tol){
        vaux = Z*p;
        vaux = G*vaux;
        vaux = vaux*Z;

        double a = r_norm*r_norm/(p*vaux);

        x += a*p;

        r += a*vaux;

        double r_norm_new = r.norm();

        double bet = r_norm_new*r_norm_new/(r_norm*r_norm);

        p*=bet;
        p-=r;

        r_norm = r_norm_new;
    }
}
inline void cg_redsys(const MatDoub2DRef &Z, const MatDoub &G, const VecDoub &c,VecDoub &x, const double tol = 1.0e-8){
    assert(Z.rows() == G.cols());
    VecDoub vaux = Z*x;
    vaux = G*vaux;
    vaux = vaux*Z;

    VecDoub r(c);
    r+=vaux;

    VecDoub p(r.size());

    for(size_t i = 0; i < p.size(); ++i) p[i] = -r[i];

    double r_norm = r.norm();

    while (r_norm > tol){
        vaux = Z*p;
        vaux = G*vaux;
        vaux = vaux*Z;

        double a = r_norm*r_norm/(p*vaux);

        x += a*p;

        r += a*vaux;

        double r_norm_new = r.norm();

        double bet = r_norm_new*r_norm_new/(r_norm*r_norm);

        p*=bet;
        p-=r;

        r_norm = r_norm_new;
    }
}
/* Algoritmo 5.2 CG  en esta implementacion
 * no se usa precondicioner (ie W = I) J. Nocedal Numerical Optimization
 * problema min(x)[ 1/2*(x_t*A*x + x_t*b ] */
inline void cg_linsys(const MatDoub &A, const VecDoub &b,VecDoub &x, const double tol = 1.0e-8){
    VecDoub vaux = A*x;

    VecDoub r(vaux - b);

    VecDoub p(r.size());

    for(size_t i = 0; i < p.size(); ++i) p[i] = -r[i];

    double r_norm = r.norm();

    while (r_norm > tol){
        vaux = A*p;

        double a = r_norm*r_norm/(p*vaux);

        x += a*p;

        r += a*vaux;

        double r_norm_new = r.norm();

        double bet = r_norm_new*r_norm_new/(r_norm*r_norm);

        p*=bet;
        p-=r;

        r_norm = r_norm_new;
    }
}
/* Algoritmo 16.1 CG precondicionado adaptado para su aplicacion a sistemas reducidos
 * en esta  J. Nocedal Numerical Optimization problema
 * min(x)[ 1/2*(x_t*Z_t*G*Z*x + x_t*c ] */
/*inline void cg_redsys_prec(const MatDoub &Z, const MatDoub &G, const VecDoub &c,VecDoub &x, const double tol = 1.0e-8){

}*/

#endif // ITERMETHODS_H
