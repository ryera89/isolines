#ifndef INTERPOLMIN_H
#define INTERPOLMIN_H
#include <cmath>

using std::pow;
using std::sqrt;

/* funcion para encontrar en minimizador de la funcion interpoladora cubica entre
 * los ptos "a0" y "a1" conociendo los valores de la funcion y de las derivadas en
 * estos puntos */
inline double min_cub_int(const double &a0, const double &a1, const double &f0, const double &f1, const double &df0, const double &df1){
    //"a0" y "a1" puntos entre los que encuentra el valor de a buscado
    //"f0" y "f1" valores de la funcion en "a0" y "a1" respectivamente
    //"df0" y "df1" valores de la derivada de la funcion en "a0" y "a1" respectivamente
    double h = a1 - a0;
    double F = f1 - f0;
    double G = (df1 - df0)*h;
    double c = G - 2.0*(F - df0*h);
    double den2 = pow(G - 3.0*c,2) - 12.0*c*df0*h;
    double den1 = G - 3.0*c;
    double den = den1 + sqrt(den2);
    double gam = - (2.0*df0*h)/den;
    return gam*h + a0;

}
/* funcion para encontrar en minimizador de la funcion interpoladora cuadratica entre
 * los ptos "a0" y "a1" conociendo los valores de la funcion en los puntos y de la derivada
 * en "a0" */
inline double min_quadratic_int(const double &a0, const double &a1, const double &f0, const double &f1, const double &df1){
    double h = a1 - a0;
    double F = f1 - f0;
    double b = df1*h - F;
    double gam = (1.0 - F/b)/2.0;
    return gam*h + a0;
}

/* funcion para encontrar en minimizador de la funcion interpoladora cuadratica entre
 * los ptos "a0" y "a1" conociendo los valores de la derivada en los puntos */
inline double secant(const double &a0, const double &a1, const double &df0, const double &df1){
    double den = df1-df0;
    if (den == 0.0) return (a0 + a1)*0.5;
    return (a0*df1 - a1*df0)/den;
}


#endif // INTERPOLMIN_H
