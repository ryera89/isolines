#ifndef STEEPESTDESCENT_H
#define STEEPESTDESCENT_H

#include "linesearch.h"

/* En esta clase esta implementado el metodo de Steepest descent donde se toma como
 * direccion de decreciemiento de la funcion objetivo pi = -Grad(f) este es un metodo
 * que puede llegar a ser muy lento solo esta implementado por razones de didacticas no
 * es recomendable para aplicaciones reales */
template<typename T>
class SteepestDescent : public Linesearch<T>{
private:
    int iter; //Numero de iteracionoes que tomo resolver el problema
    double fret; //Valor de la funcion en el minimo
    const double tol; //Tolerancia

public:
    //Los miembros de una funcion template no se heredan automaticamente
    using Linesearch<T>::xi;   //punto anterior
    using Linesearch<T>::pi;   //direccion de decrecimiento
    using Linesearch<T>::funcd; //funcion objetivo
    using Linesearch<T>::linminSWC; //metodo de busqueda lineal
    using Linesearch<T>::xmin;  //parametro encontrado
    using Linesearch<T>::df0;
    using Linesearch<T>::dfmin; //valor de la derivada direccional en el min encontrado
    using Linesearch<T>::gft; //Gradiente en xt + apt
    using Linesearch<T>::gf0; //gradiente en xt

    SteepestDescent(T &func, double toll = 3.0e-8):Linesearch<T>(func), tol(toll){}

    int IterNumber() const{return iter;}
    double FMin() const{return fret;}

    VecDoub minimize(const VecDoub &xxi, ConvAnal &conv_analizer){
        iter = 0;
        double a0 = 1.0; //valor inicial para la busqueda lineal
        const double EPS = 1.0e-18; //Para rectificar el caso especial de converger a 0.0
        const double GTOL = 1.0e-8; //criterio de convergencia para el test de cero gradiente

        size_t n = xxi.size();
        xi = xxi;

        pi.resize(n);
        gft.resize(n);
        gf0.resize(n);

        double fx = funcd(xi);  //valor de la funcion objetivo en el pto xi
        funcd.df(xi,pi);  //pi se almacena es gradiente de la funcion objetivo en el pto xi

        conv_analizer.addIter(0.0,xi,fx); //Para el analisis de convergencia

        bool is_cero_grad = true;

        for(double &a : pi){
            if (a != 0.0) is_cero_grad = false;
            a = -a;  //pi = -grad(f)
        }
        if (is_cero_grad) return xi; // por si a la primera chocamos con el minimo, muy poco probable

        for (int its = 1; true; ++its){
            iter = its;
            fret = linminSWC(a0);

            conv_analizer.addIter(xmin,xi,fret); //Para el analisis de convergencia

            //TODO: Mejorar esto
            if (2.0*abs(fret - fx) <= tol*(abs(fret) + abs(fx) + EPS)) return xi; //en caso que el minimo sea 0.0 o no haya progreso

            //TODO: Estudiar proximo estimar un nuevo a de prueba
            //a0 = 2.0*(fret - fx)/df0;
            //a0 = min(1.0,1.01*a0);

            fx = fret;
            //FIXME: Ya yo tengo el gradiente en el nuevo punto no tengo que volver a calcularlo
            funcd.df(xi,pi);
            double test = 0.0; //Para test de convergencia cero gradiente
            double den = max(fx,1.0);
            for (size_t j = 0; j < n; ++j){
                double temp = abs(pi[j]*max(abs(xi[j]),1.0))/den;
                if (temp > test) test = temp; //test viene dado por el valor mayor de los relaciones anteriores
            }
            if (test < GTOL) return xi; //test de gradiente cero

            for(double &a : pi){
               a = -a;  //pi = -grad(f)
            }
        }
    }
};
#endif // STEEPESTDESCENT_H
