#ifndef CGFRPR_H
#define CGFRPR_H

#include "linesearch.h"

/*Clase donde se implementa el metodo de gradientes conjugados
 * para problemas no lineales el algo implementado es el Fletcher-Revees
 * con las modificaciones de Polak-Riviere conocido como CG FR-PR ver
 * J Nocedal - Numerical Optimization */
template<typename T>
class CG_FR_PR : public Linesearch<T>{
private:
    int iter; //Numero de iteracionoes que tomo resolver el problema
    double fret; //Valor de la funcion en el minimo
    const double tol; //Tolerancia
    VecDoub old_pi;  //Direccion anterior
    VecDoub old_gradf; //gradiente anterior

    //TODO: Cambiar todo esto con la nueva clase de vectores q hice
    double beta_calc(){
        size_t n = pi.size();

        double num_pr = 0.0; //Para Polack-Ribiere
        double num_fr = 0.0; //Para Fletcher-Reeves
        double den1 = 0.0;

        for (size_t j = 0; j < n; ++j){
            num_pr += pi[j]*(pi[j] - old_gradf[j]); //Polack-Ribiere
            num_fr += pi[j]*pi[j];   //Fletcher-Reeves
            den1 += old_gradf[j]*old_gradf[j];
        }
        double bet_pr = num_pr/den1;
        double bet_fr = num_fr/den1;

        double bet = bet_pr;

        /* Para asegurar que la nueva direccion sea de descenso,
         *  ver J. Nocedal, Numerical optimization pag 123 */
        if (bet_pr < -bet_fr) bet = -bet_fr;
        if (bet_pr > bet_fr) bet = bet_fr;

        return bet;
    }
public:
    //Los miembros de una funcion template no se heredan automaticamente
    using Linesearch<T>::xi;   //punto anterior
    using Linesearch<T>::pi;   //direccion de decrecimiento
    using Linesearch<T>::funcd; //funcion objetivo
    using Linesearch<T>::linminSWC; //metodo de busqueda lineal
    using Linesearch<T>::xmin; //parametro encontrado
    using Linesearch<T>::df0;
    using Linesearch<T>::dfmin; //valor de la derivada direccional en el min encontrado
    using Linesearch<T>::gft; //Gradiente en xt + apt
    using Linesearch<T>::gf0; //gradiente en xt

    CG_FR_PR(T &func, double toll = 3.0e-8):Linesearch<T>(func), tol(toll){}

    int IterNumber() const{return iter;}
    double FMin() const{return fret;}

    VecDoub minimize(const VecDoub &xxi){
        iter = 0;
        double a0 = 1.0; //valor inicial para la busqueda lineal
        //const double EPS = 1.0e-18; //Para rectificar el caso especial de converger a 0.0
        const double GTOL = 1.0e-8; //criterio de convergencia para el test de cero gradiente

        size_t n = xxi.size();
        xi = xxi;

        pi.resize(n);
        gft.resize(n);
        gf0.resize(n);

        double fx = funcd(xi);  //valor de la funcion objetivo en el pto xi
        funcd.df(xi,gf0);  //pi se almacena es gradiente de la funcion objetivo en el pto xi
        if (gf0.norm() < tol) {
            fret = fx;
            return xi;
        }
        //Primera direccion de descenso
        for (size_t i = 0; i < n; ++i) pi[i] = -gf0[i];

        for (int its = 1; true; ++its){
            iter = its;
            /* c2 = 0.1: 0 < c1 < c2 < 0.5 para garantizar convergencia
             * ver J. Nocedal, Numerical Optimization */
            fret = linminSWC(a0);
            //TODO estudiar criterios de parada
            //if (2.0*abs(fret - fx) <= tol*(abs(fret) + abs(fx) + EPS)) return xi; //en caso que el minimo sea 0.0 o no haya progreso

            //TODO: Estudiar eleccion del proximo a de prueba
            //a0 = 2.0*(fret - fx)/df0; //Para el valor inicial de prueba en la proxima iteracion
            //a0 = min(1.0,1.01*a0); //ajuste del valor inicial de prueba

            fx = fret;
            double test = 0.0; //Para test de convergencia cero gradiente
            //double den = max(fx,1.0);
            for (size_t j = 0; j < n; ++j){
                //double temp = abs(pi[j]*max(abs(xi[j]),1.0))/den;
                double temp = abs(gft[j]);
                if (temp > test) test = temp; //Norma infinito
            }
            //TODO: Estudiar test de gradiente cero
            if (test < GTOL*(1.0+abs(fx))) return xi; //test de gradiente cero
            double den = gf0.norm_2();
            double bet_pr = gft*(gft - gf0)/den;
            double bet_fr = gft.norm_2()/den;
            double bet = bet_pr;
            /* Para asegurar que la nueva direccion sea de descenso,
             * ver J. Nocedal, Numerical optimization pag 123 */
            if (bet_pr < -bet_fr) bet = -bet_fr;
            if (bet_pr > bet_fr) bet = bet_fr;
            //La nueva direccion de descenso
            pi = bet*pi - gft;
        }
    }
    VecDoub minimize(const VecDoub &xxi, ConvAnal &conv_analizer){
        iter = 0;
        double a0 = 1.0; //valor inicial para la busqueda lineal
        //const double EPS = 1.0e-18; //Para rectificar el caso especial de converger a 0.0
        const double GTOL = 1.0e-8; //criterio de convergencia para el test de cero gradiente

        size_t n = xxi.size();
        xi = xxi;

        pi.resize(n);
        gft.resize(n);
        gf0.resize(n);

        double fx = funcd(xi);  //valor de la funcion objetivo en el pto xi
        funcd.df(xi,gf0);  //pi se almacena es gradiente de la funcion objetivo en el pto xi

        conv_analizer.addIter(0.0,xi,fx); //Para el analisis de convergencia

        if (gf0.norm() < tol) {
            fret = fx;
            return xi;
        }
        //Primera direccion de descenso
        for (size_t i = 0; i < n; ++i) pi[i] = -gf0[i];

        for (int its = 1; true; ++its){
            iter = its;
            /* c2 = 0.1: 0 < c1 < c2 < 0.5 para garantizar convergencia
             * ver J. Nocedal, Numerical Optimization */
            fret = linminSWC(a0);
            conv_analizer.addIter(xmin,xi,fret); //Para el analisis de convergencia
            //TODO estudiar criterios de parada
            //if (2.0*abs(fret - fx) <= tol*(abs(fret) + abs(fx) + EPS)) return xi; //en caso que el minimo sea 0.0 o no haya progreso

            //TODO: Estudiar eleccion del proximo a de prueba
            //a0 = 2.0*(fret - fx)/df0; //Para el valor inicial de prueba en la proxima iteracion
            //a0 = min(1.0,1.01*a0); //ajuste del valor inicial de prueba

            fx = fret;
            double test = 0.0; //Para test de convergencia cero gradiente
            //double den = max(fx,1.0);
            for (size_t j = 0; j < n; ++j){
                //double temp = abs(pi[j]*max(abs(xi[j]),1.0))/den;
                double temp = abs(gft[j]);
                if (temp > test) test = temp; //Norma infinito
            }
            //TODO: Estudiar test de gradiente cero
            if (test < GTOL*(1.0+abs(fx))) return xi; //test de gradiente cero
            double den = gf0.norm_2();
            double bet_pr = gft*(gft - gf0)/den;
            double bet_fr = gft.norm_2()/den;
            double bet = bet_pr;
            /* Para asegurar que la nueva direccion sea de descenso,
             * ver J. Nocedal, Numerical optimization pag 123 */
            if (bet_pr < -bet_fr) bet = -bet_fr;
            if (bet_pr > bet_fr) bet = bet_fr;
            //La nueva direccion de descenso
            pi = bet*pi - gft;
        }
    }
};
#endif // CGFRPR_H
