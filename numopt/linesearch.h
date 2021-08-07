#ifndef LINESEARCH_H
#define LINESEARCH_H

#include "linalg/mathmatrix.h"
#include "interpolmin.h"
#include "convanalizer.h"
#include <cmath>

using std::size_t;
using std::abs;
using std::pow;
using std::max;
using std::min;
using std::sqrt;


/* Functor intermediario entre la funcion objetivo y los
 * metodos de busqueda lineal, aqui se representa computacionalmente
 * la funcion unidimensional G(a) = f(xi + a*pi) (xi: punto de evaluacion)
 * pi: direccion de busqueda(descenso) */
template<typename T>
struct DF1dim
{
    const VecDoub &xi; //punto inicial
    const VecDoub &pi; //direccion de busqueda(descenso)
    VecDoub &dft;  //gradiente  en el nuevo punto xt
    VecDoub &dfi;  //gradiente  en el nuevo punto anterior
    T &funcd;     //Funcion objetivo
    size_t n;     //Numero de dimensiones
    VecDoub xt;   //Nuevo punto de evaluacion
    //FIXME: estas 2 variables son por referencia y pertenecen a la clase de la busqueda lineal

    DF1dim(const VecDoub &xxi, const VecDoub &ppi,VecDoub &ddft,VecDoub &ddfi,T &funcdd):xi(xxi),pi(ppi),
        dft(ddft),dfi(ddfi),funcd(funcdd),n(xi.size()),xt(n){}
    double operator()(const double &a){
        xt = xi + a*pi;     //Calcula en nuevo punto
        return funcd(xt);   //Retorna el valor de la funcion objetivo en el nuevo punto
    }
    double df(){
        funcd.df(xt,dft);   //calcula el gradiente de f y lo guarda en el vector dft
        return dft*pi; //derivada direccional
    }
    double df0(){
        funcd.df(xi,dfi);
        return dfi*pi;
    }
};
/* Clase Linesearch es la base de muchas otras clases donde se implementan metodos de
 * optimizacion aqui se implementan los metodos de busqueda lineal que son parte imporante
 * de muchos algoritmos de optimizacion*/
template<typename T>
class Linesearch{
private:

    //Rutinas auxiliares para la rutina de linminMWC
    bool MWC(const double &t, const double &ft,const double &dft, const double &DELT,
             const double &SIG,const double &EPSK){
        bool wc = (ft - f0 <= DELT*t*df0) && (dft >= SIG*df0);
        bool mwc = ((2.0*DELT - 1.0)*df0 >= dft) && (dft >= SIG*df0) && (ft <= f0 + EPSK);

        return wc || mwc;
    }
    template<typename F>
    void updateInt(double &a,double &b,double &dfa,double &dfb,double &t,double &ft,
                   double &dft, const double &THET, const double &EPSK,F &func){
        if (t <= a || t >= b) return; //Si c esta fuera del intervalo

        if (dft >= 0.0){
            b = t;
            dfb = dft;
            //t = c;
            return;
        }

        if (ft <= f0 + EPSK){
            a = t;
            dfa = dft;
            //t = c;
            return;
        }

        double aprim = a;
        double bprim = t;

        while(true){
            t = (1.0-THET)*aprim + THET*bprim;
            ft = func(t);
            dft = func.df();
            if (dft >= 0.0){
                b = t;
                a = aprim;
                dfb = dft;
                return;
            }
            if (ft <= f0 + EPSK){ aprim = t; dfa = dft;}
            else {bprim = t; dfb = dft;}
        }
    }
    template<typename F>
    void secant2(double &a, double &b,double &dfa,double &dfb,double &t,
                 double &ft, double &dft,const double &THET,const double &EPSK,F &func){

        double A = a;
        double B = b;
        double DFA = dfa;
        double DFB = dfb;
        t = secant(A,B,DFA,DFB);

        double c = t;

        ft = func(t);
        dft = func.df();

        updateInt(A,B,DFA,DFB,t,ft,dft,THET,EPSK,func);

        if (c != A && c != B){
            a = A;
            b = B;
            dfa = DFA;
            dfb = DFB;
            return;
        }

        if (c == B) t = secant(b,B,dfb,DFB);
        if (c == A) t = secant(a,A,dfa,DFA);

        a = A;
        b = B;

        ft = func(t);
        dft = func.df();

        updateInt(a,b,dfa,dfb,t,ft,dft,THET,EPSK,func);
    }
protected:
    VecDoub xi; //punto anterior
    VecDoub pi; //direccion de busqueda
    VecDoub gft; //gradiente de f en el nuevo pto
    VecDoub gf0; //fradiente en el pto anterior
    T &funcd; //funcion objetivo
    double xmin; //valor del parametro a optimo encontrado
    double fmin; //valor de la funcion objetivo en este punto
    double f0;
    double df0;
    double dfmin;
public:
    Linesearch(T &funcc) : funcd(funcc){}

    /* Algoritmo recomendado por J. Nocedal y S.J. Wright, en su libro
     * Numerical optimization 2ed 2006, y que viene descrito en el paper
     * Line search with guaranteed sufficient decrease de J.J. More y D.J.Thuente */
    double linminSWC(const double &a0,const double R = 0.1, const double M = 0.0001){
        int interval_generated_counter = 0;
        double al = 0.0;
        double au;
        double at = a0;    //Valor de prueba inicial
        double fl,fu,ft;   //Valores de la funcion f(xi + api)
        double dfl,dfu,dft; //Valores la derivada de la funcion g(xi + a*pi)
        double gl,gu,gt;  //Valores de la funcion auxiliar g(xi + a*pi) = f(xi + a*pi) - f(0) - M*a*
        double dgl,dgu,dgt; //Valores de las derivadas de la funcion auxiliar en los ptos

        double atemp;

        double ac;  //Valor auxiliar para minimizador de interpolacion cubica
        double aq, as; //Valor auxiliar para minimizador interpolacion cuadratica

        bool funtion_switcher = false; //Badera para usar una funcion u otra al interpolar
        bool U2Case = true; //Flag para si no se ha brackeado la solucion

        DF1dim<T> fd1dim(xi,pi,gft,gf0,funcd);

        f0 = funcd(xi);  //Valor de la funcion ovjetivo en el pto inicial xi
        df0 = fd1dim.df0(); //Derivada direccional de la funcion objetivo en para f0

        fl = f0;   //"al" al inicio es igual a 0.0
        dfl = df0;
        gl = 0.0;
        dgl = dfl - M*df0;

        ft = fd1dim(at);    //Valor de la funcion ojetivo en el nuevo pto xi + at*pi
        dft = fd1dim.df();  //Valor de la derivada de la funcion objetivo en el nuevo pto xi + at*pi
        gt = ft - f0 - M*at*df0;
        dgt = dft - M*df0;

        //variables auxiliares para monitorear el la longitud del intervalo
        double step_2,step_1,step_0;

        bool apply_bisec = false;    //Bandera para aplicar biseccion

        for (size_t i = 1; true; ++i){

            if ( ft <= f0 + M*at*df0 && abs(dft) <= R*abs(df0)){ //Strong Wolfe: Term Test
                xmin = at;     //paso encontrado en la busqueda lineal
                xi = fd1dim.xt; //actualizamos el punto
                dfmin = dft;
                return fmin = ft; //retornamos el valor de la funcion
            }

            if (funtion_switcher || (gt <= 0.0 && dft > 0.0 )){ //Usando la funcion f para hallar el proximo pto de prueba
                if (U2Case && ft <= fl && dft*(al - at) > 0.0){
                    atemp = at + 4.0*(at - al); //Proximo paso a evaluar
                    al = at;
                    fl = ft;
                    dfl = dft;
                    at = atemp; //nueva prueba
                    ft = fd1dim(at);    //Valor de la funcion ojetivo en el nuevo pto xi + at*pi
                    dft = fd1dim.df();  //Valor de la derivada de la funcion objetivo en el nuevo pto xi + at*pi
                    continue;
                }
                /* Caso 1: minimizador de la funcion de interpolacion cubica
                     *  entre ft,fl,dft,dfl y la funcion de interpolacion cuadratica ft,fl,dfl */
                if (ft > fl){
                    aq = min_quadratic_int(at,al,ft,fl,dfl);  //minimizador de interpolacion cuadratica
                    ac = min_cub_int(at,al,ft,fl,dft,dfl); //minimizador de interpolacion cubica
                    atemp = (abs(ac - al) < abs(aq - al)) ? ac : 0.5*(ac + aq);
                    au = at;
                    fu = ft;
                    dfu = dft;
                    U2Case = false;
                }else{
                    /* Caso 2 minimizador de la funcion de interpolacion cubica
                         * usando los valores ft,fl,dft,dfl y la funcion de interpolacion
                         * cuadratica usando los valores de ft,dft,dfl */
                    if (dft*dfl < 0.0){
                        as = secant(al,at,dfl,dft);  //secante
                        ac = min_cub_int(al,at,fl,ft,dfl,dft); //minimizador de interpolacion cubica
                        atemp = (abs(ac - at) >= abs(as - at)) ? ac : as;
                        au = al;
                        al = at;
                        fu = fl;
                        dfu = dfl;
                        fl = ft;
                        dfl = dft;
                        U2Case = false;
                    }else{
                        //BUG: Posible bug
                        /* Caso 3 minimizador de la funcion de interpolacion cubica
                             * usando los valores ft,fl,dft,dfl y la funcion de interpolacion
                             * cuadratica usando los valores de ft,dft,dfl */
                        if (abs(dft) <= abs(dfl)){
                            as = secant(al,at,dfl,dft);  //secante
                            ac = min_cub_int(al,at,fl,ft,dfl,dft); //minimizador de interpolacion cubica
                            atemp = (abs(ac - at) < abs(as - at)) ? ac : as;
                            double aux = at + 0.66*(au - at);

                            if (at > al) atemp = min(aux,atemp);
                            else atemp = max(aux,atemp);
                        }else{
                            /* Caso 4: minimizador de la funcion de
                                 * interpolacion cubica entre usando los valores
                                 * ft,fu,dft,dfu */
                            atemp = min_cub_int(au,at,fu,ft,dfu,dft);
                        }
                        al = at;
                        fl = ft;
                        dfl = dft;
                    }
                }
                ++interval_generated_counter;
                step_2 = step_1;
                step_1 = step_0;
                step_0 = abs(au - al);

                /* si no ha decrecido el intervalo aplicamos
                 *  al proximo biseccion */
                apply_bisec = (interval_generated_counter > 2 && step_0 > 0.66*step_2);

                if (apply_bisec){
                    at = (au + al)*0.5;
                    interval_generated_counter = 0; //Reseteo el contador
                }else{
                    at = atemp;
                }

                ft = fd1dim(at);    //Valor de la funcion ojetivo en el nuevo pto xi + at*pi
                dft = fd1dim.df();  //Valor de la derivada de la funcion objetivo en el nuevo pto xi + at*pi

                //siempre se va a a utilizar la funcion f
                funtion_switcher = true;
            }else{ // Usando la funcion auxiliar g
                if (U2Case && gt <= gl && dgt*(al - at) > 0.0){
                    atemp = at + 4.0*(at - al);
                    al = at;
                    fl = ft;
                    dfl = dft;
                    gl = gt;
                    dgl = dgt;
                    at = atemp; //nueva prueba
                    ft = fd1dim(at);    //Valor de la funcion ojetivo en el nuevo pto xi + at*pi
                    dft = fd1dim.df();  //Valor de la derivada de la funcion objetivo en el nuevo pto xi + at*pi
                    //Funcion auxiliar en "at"
                    gt = ft - f0 - M*at*df0;
                    dgt = dft - M*df0;

                    continue;
                }
                /* Caso 1: minimizador de la funcion de interpolacion cubica
                     *  entre gt,gl,dgt,dgl y la funcion de interpolacion cuadratica gt,gl,dgl */
                if (gt > gl){
                    aq = min_quadratic_int(at,al,gt,gl,dgl);  //minimizador de interpolacion cuadratica
                    ac = min_cub_int(at,al,gt,gl,dgt,dgl); //minimizador de interpolacion cubica
                    atemp = (abs(ac - al) < abs(aq - al)) ? ac : 0.5*(ac + aq);
                    au = at;
                    fu = ft;
                    dfu = dft;
                    gu = gt;
                    dgu = dgt;
                    U2Case = false;
                }else{
                    /* Caso 2 minimizador de la funcion de interpolacion cubica
                         * usando los valores gt,gl,dgt,dgl y la funcion de interpolacion
                         * cuadratica usando los valores de gt,dgt,dgl */
                    if (dgt*dgl < 0.0){
                        as = secant(al,at,dgl,dgt);  //secante
                        ac = min_cub_int(al,at,gl,gt,dgl,dgt); //minimizador de interpolacion cubica
                        atemp = (abs(ac - at) >= abs(as - at)) ? ac : as;
                        au = al;
                        al = at;

                        fu = fl;
                        dfu = dfl;
                        gu = gl;
                        dgu = dgl;

                        fl = ft;
                        dfl = dft;
                        gl = gt;
                        dgl = dgt;

                        U2Case = false;
                    }else{
                        //BUG: Posible bug
                        /* Caso 3 minimizador de la funcion de interpolacion cubica
                             * usando los valores gt,gl,dgt,dgl y la funcion de interpolacion
                             * cuadratica usando los valores de gt,dgt,dgl */
                        if (abs(dgt) <= abs(dgl)){
                            as = secant(al,at,dgl,dgt);  //secante
                            ac = min_cub_int(al,at,gl,gt,dgl,dgt); //minimizador de interpolacion cubica
                            atemp = (abs(ac - at) < abs(as - at)) ? ac : as;
                            double aux = at + 0.66*(au - at);

                            if (at > al) atemp = min(aux,atemp);
                            else atemp = max(aux,atemp);
                        }else{
                            /* Caso 4: minimizador de la funcion de
                                 * interpolacion cubica entre usando los valores
                                 * gt,gu,dgt,dgu */
                            atemp = min_cub_int(au,at,gu,gt,dgu,dgt);
                        }
                        al = at;

                        fl = ft;
                        dfl = dft;

                        gl = gt;
                        dgl = dgt;
                    }
                }
                ++interval_generated_counter;
                step_2 = step_1;
                step_1 = step_0;
                step_0 = abs(au - al);

                /* si no ha decrecido el intervalo aplicamos
                 *  al proximo biseccion */
                apply_bisec = (interval_generated_counter > 2 && step_0 > 0.66*step_2);

                if (apply_bisec){
                    at = (au + al)*0.5;
                    interval_generated_counter = 0; //Reseteo el contador
                }else{
                    at = atemp;
                }

                ft = fd1dim(at);    //Valor de la funcion objetivo en el nuevo pto xi + at*pi
                dft = fd1dim.df();  //Valor de la derivada direccional en el nuevo pto xi + at*pi

                //Funcion auxiliar en "at"
                gt = ft - f0 - M*at*df0;
                dgt = dft - M*df0;
            }
        }
    }
    double linminMWC(const double &a0, const double DELT = 0.1, const double SIG = 0.9,
                     const double EPS = 1.0e-6, const double THET = 0.5, const double &GAM = 0.66){
        double a = 0.0,b; //limites del intervalo donde se realiza la busqueda
        double t = a0; //Primer valor de trial
        //double fa,fb;
        double ft;
        double dfa,dft,dfb;

        DF1dim<T> df1dim(xi,pi,gft,gf0,funcd);
        f0 = funcd(xi);
        df0 = df1dim.df0();

        dfa = df0;

        double EPSK = EPS*abs(f0);

        while(true){
            ft = df1dim(t);
            dft = df1dim.df();
            if (MWC(t,ft,dft,DELT,SIG,EPSK)){
                xmin = t;     //paso encontrado en la busqueda lineal
                xi = df1dim.xt; //actualizamos el punto
                dfmin = dft;
                return fmin = ft; //retornamos el valor de la funcion
            }
            if (dft >= 0.0){
                b = t;
                //fb = ft;
                dfb = dft;
                break;

            }
            t *= 2.1;
        }
        double length = 0.0;
        for(size_t i = 0; i < 10000; ++i){
            length = b-a;
            secant2(a,b,dfa,dfb,t,ft,dft,THET,EPSK,df1dim);
            if ( b-a > GAM*(length)){
                t = (a+b)*0.5;
                ft = df1dim(t);
                dft = df1dim.df();
                updateInt(a,b,dfa,dfb,t,ft,dft,THET,EPSK,df1dim);
            }
            if (MWC(t,ft,dft,DELT,SIG,EPSK)){
                xmin = t;     //paso encontrado en la busqueda lineal
                xi = df1dim.xt; //actualizamos el punto
                dfmin = dft;
                return fmin = ft; //retornamos el valor de la funcion
            }
        }
       throw ("Linesearch<T>::linminMWC : Failed, maximo numero de iteraciones alcanzadas");
    }
};

#endif // LINESEARCH_H
