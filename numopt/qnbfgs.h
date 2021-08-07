#ifndef QNBFGS_H
#define QNBFGS_H

#include "linesearch.h"
#include "convanalizer.h"

/*Metodo BFGS uno de los metodos Quasi-Newton mas populares algoritmo
 * presentado por J Nocedal y S.J. Wright en Numerical Optimization */
template<typename T>
class QN_BFGS : public Linesearch<T>{
private:
    int iter;
    const double tol;
    double fret;
public:
    using Linesearch<T>::xi;   //punto anterior
    using Linesearch<T>::pi;   //direccion de decrecimiento
    using Linesearch<T>::funcd; //funcion objetivo
    using Linesearch<T>::linminSWC; //metodo de busqueda lineal
    using Linesearch<T>::linminMWC; //metodo de busqueda lineal
    using Linesearch<T>::xmin; //parametro encontrado
    using Linesearch<T>::df0;
    using Linesearch<T>::dfmin; //valor de la derivada direccional en el min encontrado
    using Linesearch<T>::gft; //Gradiente en xt + apt
    using Linesearch<T>::gf0; //gradiente en xt


    QN_BFGS(T &func, double toll = 1.0e-8):Linesearch<T>(func), tol(toll){}


    int IterNumber() const{return iter;}
    double FMin() const{return fret;}

    VecDoub minimize(const VecDoub &xxi){
        size_t n = xxi.size();

        MatDoub H(n,n,0.0);
        for (size_t i = 0; i < n; ++i) H(i,i) = 1.0;

        xi = xxi;
        VecDoub xi_old = xxi; //Variable para almacenar punto anterior
        pi.resize(n);
        gft.resize(n);
        gf0.resize(n);
        funcd.df(xi,pi); //Gradiente de la funcion objetivo
        VecDoub g_old = pi; //Variable para almacenar el gradiente anterior
        iter = 0; //contador de iteraciones

        VecDoub sk(n);
        VecDoub yk(n);

        MatDoub skyk;
        MatDoub yksk;
        MatDoub sksk;
        while (pi.norm() > tol){
            pi = H*pi*(-1.0);

            fret = linminSWC(1.0);

            //FIXME: No hay que volver a calcular el gradiente
            funcd.df(xi,pi);
            sk = xi - xi_old;
            yk = pi - g_old;

            if (iter == 0) H *= yk*sk/yk.norm_2();

            double rho = 1.0/(yk*sk);

            yk *= -rho;
            skyk = prod_t(sk,yk);
            yksk = prod_t(yk,sk);
            sksk = prod_t((rho)*sk,sk);

            for (size_t i = 0; i < n; ++i) {
                skyk(i,i)+=1.0;
                yksk(i,i)+=1.0;
            }

            H = skyk*H*yksk + sksk;

            //g_old no es con pi es el gradiente calculando anterior
            xi_old = xi;
            g_old = pi;
            ++iter;
        }
        return xi;
    }
    VecDoub minimize(const VecDoub &xxi,ConvAnal &conv_analizer){
        size_t n = xxi.size();

        MatDoub H(n,n,0.0);
        for (size_t i = 0; i < n; ++i) H(i,i) = 1.0;

        xi = xxi;
        VecDoub xi_old = xxi; //Variable para almacenar punto anterior
        pi.resize(n);
        gft.resize(n);
        gf0.resize(n);
        double fx = funcd(xi);
        funcd.df(xi,pi); //Gradiente de la funcion objetivo
        VecDoub g_old = pi; //Variable para almacenar el gradiente anterior
        iter = 0; //contador de iteraciones

        VecDoub sk(n);
        VecDoub yk(n);

        MatDoub skyk;
        MatDoub yksk;
        MatDoub sksk;

        conv_analizer.addIter(0.0,xi,fx); //Para el analisis de convergencia
        while (pi.norm() > tol){
            pi = H*pi*(-1.0);

            fret = linminSWC(1.0);

            conv_analizer.addIter(xmin,xi,fret); //Para el analisis de convergencia

            //FIXME: No hay que volver a calcular el gradiente
            funcd.df(xi,pi);
            sk = xi - xi_old;
            yk = pi - g_old;

            if (iter == 0) H *= yk*sk/yk.norm_2();

            double rho = 1.0/(yk*sk);

            yk *= -rho;
            skyk = prod_t(sk,yk);
            yksk = prod_t(yk,sk);
            sksk = prod_t((rho)*sk,sk);

            for (size_t i = 0; i < n; ++i) {
                skyk(i,i)+=1.0;
                yksk(i,i)+=1.0;
            }

            H = skyk*H*yksk + sksk;

            //g_old no es con pi es el gradiente calculando anterior
            xi_old = xi;
            g_old = pi;
            ++iter;
        }
        return xi;
    }
};

#endif // QNBFGS_H
