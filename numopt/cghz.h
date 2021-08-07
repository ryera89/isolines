#ifndef CGHZ_H
#define CGHZ_H

#include "linesearch.h"

/*Clase donde se implementa el metodo de gradientes conjugados
 * para problemas no lineales ver J. Nocedal, Numerical Optimization
 * W. W. HAGER AND H. ZHANG, A new conjugate gradient method
 * with guaranteed descent and an efficient line search,
 * SIAM Journal on Optimization, 16 (2005), pp. 170–192. */
template<typename T>
class CG_HZ : public Linesearch<T>{
private:
    int iter; //Numero de iteracionoes que tomo resolver el problema
    double fret; //Valor de la funcion en el minimo
    const double tol; //Tolerancia
public:
    //Los miembros de una funcion template no se heredan automaticamente
    using Linesearch<T>::xi;   //punto anterior
    using Linesearch<T>::pi;   //direccion de decrecimiento
    using Linesearch<T>::funcd; //funcion objetivo
    using Linesearch<T>::linminMWC; //metodo de busqueda lineal
    using Linesearch<T>::xmin; //parametro encontrado
    using Linesearch<T>::df0;
    using Linesearch<T>::dfmin; //valor de la derivada direccional en el min encontrado
    using Linesearch<T>::gft; //Gradiente en xt + apt
    using Linesearch<T>::gf0; //gradiente en xt


    CG_HZ(T &func, double toll = 3.0e-8):Linesearch<T>(func), tol(toll){}

    int IterNumber() const{return iter;}
    double FMin() const{return fret;}

    VecDoub minimize(const VecDoub &xxi, const double ETA = 0.01){
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
            fret = linminMWC(a0);
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
            //TODO: MEJORAR ESTO
            if (test < GTOL*(1.0+abs(fx))) return xi; //test de gradiente cero

            VecDoub vec_aux = gft - gf0;
            double aux1 = vec_aux.norm_2();
            double aux2 = vec_aux*pi;

            double old_gradf_mod = gf0.norm();
            double old_pi_mod = pi.norm();

            double ETAK = -1.0/(old_pi_mod*min(ETA,old_gradf_mod));

            double aux3 = 2.0*aux1/aux2;

            double bet = (vec_aux - aux3*pi)*gft;

            bet/=aux2;

            bet = max(bet,ETAK);

            //Creo que old_pi no hace falta ya
            //La nueva direccion de descenso
            pi = bet*pi - gft;
        }
    }
    VecDoub minimize(const VecDoub &xxi, ConvAnal &conv_analizer, const double ETA =  0.01){
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
            fret = linminMWC(a0);
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

            //TODO: MEJORAR ESTO
            if (test < GTOL*(1.0+abs(fx))) return xi; //test de gradiente cero

            VecDoub vec_aux = gft - gf0;
            double aux1 = vec_aux.norm_2();
            double aux2 = vec_aux*pi;

            double old_gradf_mod = gf0.norm();
            double old_pi_mod = pi.norm();

            double ETAK = -1.0/(old_pi_mod*min(ETA,old_gradf_mod));

            double aux3 = 2.0*aux1/aux2;

            double bet = (vec_aux - aux3*pi)*gft;

            bet/=aux2;

            bet = max(bet,ETAK);

            //Creo que old_pi no hace falta ya
            //La nueva direccion de descenso
            pi = bet*pi - gft;
        }
    }
};
#endif // CGHZ_H
