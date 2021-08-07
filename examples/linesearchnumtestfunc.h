#ifndef LINESEARCHNUMTESTFUNC_H
#define LINESEARCHNUMTESTFUNC_H

#include <cmath>


using std::size_t;
using std::abs;
using std::pow;
using std::max;
using std::min;
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
    //Calculo del minimizador de la funcion interpoladora cubica
    /*double A = ((df1 + df0)*(a1-a0) - 2.0*(f1 - f0))/pow(a1 - a0,3);
    double B1 = 3.0*(f1 - f0)*(a1 + a0)/pow(a1 - a0,3);
    double B2 = (a1*(df1 + 2.0*df0) + a0*(df0 + 2.0*df1))/pow(a1 - a0,2);
    double B = B1 - B2;

    double C = df1 - a1*(2.0*B + 3.0*a1*A);

    double D = B*B - 3.0*A*C; //Discriminante de la funcion "C + 2B*x + 3*A*x^2"

    if (D < 0.0) return 0.0; //Comprobando el valor del discriminante,
                             //si es menor que cero retorno 0

    //Valor del minimizador de la interpolacion cubica
    double ac1 = (-B + sqrt(D))/(3.0*A);
    double ac2 = (-B - sqrt(D))/(3.0*A);

    //Chequeo de que el valor encontrado de "a" se encuentre entre "a0" y "a1"
    if ((ac1 > a0 && ac1 < a1) || (ac1 < a0 && ac1 > a1)) return ac1;

    return ac2;*/

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
    //return (2.0*a0*(f1 - f0) - df0*(a1*a1 - a0*a0))/(2.0*(f1 - f0 - df0*(a1 - a0)));
}

/* funcion para encontrar en minimizador de la funcion interpoladora cuadratica entre
 * los ptos "a0" y "a1" conociendo los valores de la derivada en los puntos */
inline double secant(const double &a0, const double &a1, const double &df0, const double &df1){

    return (a0*df1 - a1*df0)/(df1 - df0);
}


/* Funciones Para prueba numerica */
struct Fun1{
    double operator()(const double &a){
        double b = 2.0;
        return -a/(a*a + b);
    }
    double df(const double &a){
        double b = 2.0;
        return (a*a - b)/pow(a*a + b,2.0);
    }
};
struct Fun2{
    double operator()(const double &a){
        double b = 0.004;
        return pow(a + b,5.0) - 2.0*pow(a + b,4.0);
    }
    double df(const double &a){
        double b = 0.004;
        return 5.0*pow(a + b,4.0) - 8.0*pow(a + b,3.0);
    }
};
struct Fun3{
    double operator()(const double &a){
        double b = 0.01;
        double l = 39.0;
        double pi = 4.0*atan(1.0);

        double aux = pow(a - 1,2.0)/(2.0*b) + 0.5*b;

        if (a <= 1-b) aux = 1 - a;
        if (a >= 1+b) aux = a - 1;

        return aux + 2.0*(1 - b)/(l*pi)*sin(l*pi*a*0.5);
    }
    double df(const double &a){
        double b = 0.01;
        double l = 39.0;
        double pi = 4.0*atan(1.0);

        double aux = (a-1)/b;

        if (a <= 1-b) aux = -1;
        if (a >= 1+b) aux = 1;

        return aux + (1 - b)*cos(l*pi*a*0.5);
    }
};

struct Fun4{
    double operator()(const double &a){
        double b1 = 0.001;
        double b2 = 0.001;

        double gam1 = sqrt(1.0 + b1*b1) - b1;
        double gam2 = sqrt(1.0 + b2*b2) - b2;

        double aux2 = sqrt(pow(1.0-a,2.0) + b2*b2);
        double aux1 = sqrt(a*a + b1*b1);

        return gam1*aux2 + gam2*aux1;
    }
    double df(const double &a){
        double b1 = 0.001;
        double b2 = 0.001;

        double gam1 = sqrt(1.0 + b1*b1) - b1;
        double gam2 = sqrt(1.0 + b2*b2) - b2;

        double aux2 = sqrt(pow(1.0-a,2.0) + b2*b2);
        double aux1 = sqrt(a*a + b1*b1);

        return a*gam2/aux1 - (1.0 - a)*gam1/aux2;

    }
};
struct Fun5{
    double operator()(const double &a){
        double b1 = 0.01;
        double b2 = 0.001;

        double gam1 = sqrt(1.0 + b1*b1) - b1;
        double gam2 = sqrt(1.0 + b2*b2) - b2;

        double aux2 = sqrt(pow(1.0-a,2.0) + b2*b2);
        double aux1 = sqrt(a*a + b1*b1);

        return gam1*aux2 + gam2*aux1;
    }
    double df(const double &a){
        double b1 = 0.01;
        double b2 = 0.001;

        double gam1 = sqrt(1.0 + b1*b1) - b1;
        double gam2 = sqrt(1.0 + b2*b2) - b2;

        double aux2 = sqrt(pow(1.0-a,2.0) + b2*b2);
        double aux1 = sqrt(a*a + b1*b1);

        return a*gam2/aux1 - (1.0 - a)*gam1/aux2;

    }
};
struct Fun6{
    double operator()(const double &a){
        double b1 = 0.001;
        double b2 = 0.01;

        double gam1 = sqrt(1.0 + b1*b1) - b1;
        double gam2 = sqrt(1.0 + b2*b2) - b2;

        double aux2 = sqrt(pow(1.0-a,2.0) + b2*b2);
        double aux1 = sqrt(a*a + b1*b1);

        return gam1*aux2 + gam2*aux1;
    }
    double df(const double &a){
        double b1 = 0.001;
        double b2 = 0.01;

        double gam1 = sqrt(1.0 + b1*b1) - b1;
        double gam2 = sqrt(1.0 + b2*b2) - b2;

        double aux2 = sqrt(pow(1.0-a,2.0) + b2*b2);
        double aux1 = sqrt(a*a + b1*b1);

        return a*gam2/aux1 - (1.0 - a)*gam1/aux2;

    }
};
template<typename T>
double linminnumtest(const double &a0, const double &R, const double &M, size_t &iter, double &amin, double &fmin, double &dfmin, T &funcd){
    int interval_generated_counter = 0;
    double al = 0.0;
    double au;
    double df0;
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
    //DF1dim<T> fd1dim(xi,pi,funcd);

    const double f0 = funcd(0.0);  //Valor de la funcion ovjetivo en el pto inicial xi
    df0 = funcd.df(0.0); //Valor de la derivada de la funcion objetivo en pto inicial xi

    fl = f0;   //"al" al inicio es igual a 0.0
    dfl = df0;
    gl = 0.0;
    dgl = dfl - M*df0;

    ft = funcd(at);    //Valor de la funcion ojetivo en el nuevo pto xi + at*pi
    dft = funcd.df(at);  //Valor de la derivada de la funcion objetivo en el nuevo pto xi + at*pi
    gt = ft - f0 - M*at*df0;
    dgt = dft - M*df0;

    double step_2,step_1,step_0;
    //step_2=step_1=step_0 = au;   //variables auxiliares para monitorear el la longitud del nuevo intervalo

    bool apply_bisec = false;    //Bandera para aplicar biseccion

    for (size_t i = 1; true; ++i){

        if ( ft <= f0 + M*at*df0 && abs(dft) <= R*abs(df0)){ //Condicion fuerte de WOLFE
            iter = i;
            amin = at;
            fmin = ft;
            dfmin = dft;
            return ft;
        }

        funtion_switcher = (funtion_switcher || (gt <= 0.0 && dft > 0.0 ));

        if (funtion_switcher){ //Usando la funcion f para hallar el proximo pto de prueba
            if (U2Case && ft <= fl && dft*(al - at) > 0.0){
                atemp = at + 4.0*(at - al); //Proximo paso a evaluar
                al = at;
                fl = ft;
                dfl = dft;
                at = atemp; //nueva prueba
                ft = funcd(at);    //Valor de la funcion ojetivo en el nuevo pto xi + at*pi
                dft = funcd.df(at);  //Valor de la derivada de la funcion objetivo en el nuevo pto xi + at*pi
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
                    /* Caso 3 minimizador de la funcion de interpolacion cubica
                         * usando los valores ft,fl,dft,dfl y la funcion de interpolacion
                         * cuadratica usando los valores de ft,dft,dfl */
                    if (abs(dft) <= abs(dfl)){
                        as = secant(al,at,dfl,dft);  //secante
                        ac = min_cub_int(al,at,fl,ft,dfl,dft); //minimizador de interpolacion cubica
                        atemp = (abs(ac - at) < abs(as - at)) ? ac : as;
                        double aux = at + 0.66*(au - at);
                        //WARNING puede ser un Bug
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
        }else{ // Usando la funcion auxiliar g para hallar el proximo pto de prueba
            if (U2Case && gt <= gl && dgt*(al - at) > 0.0){
                atemp = at + 4.0*(at - al);
                al = at;
                fl = ft;
                dfl = dft;
                gl = gt;
                dgl = dgt;
                at = atemp; //nueva prueba
                ft = funcd(at);    //Valor de la funcion ojetivo en el nuevo pto xi + at*pi
                dft = funcd.df(at);  //Valor de la derivada de la funcion objetivo en el nuevo pto xi + at*pi
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
                    /* Caso 3 minimizador de la funcion de interpolacion cubica
                         * usando los valores gt,gl,dgt,dgl y la funcion de interpolacion
                         * cuadratica usando los valores de gt,dgt,dgl */
                    if (abs(dgt) <= abs(dgl)){
                        as = secant(al,at,dgl,dgt);  //secante
                        ac = min_cub_int(al,at,gl,gt,dgl,dgt); //minimizador de interpolacion cubica
                        atemp = (abs(ac - at) < abs(as - at)) ? ac : as;
                        double aux = at + 0.66*(au - at);
                        //WARNING: Cuidado con esto posible bug, estudiar bien los casos
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
            interval_generated_counter = 0;
        }else{
            at = atemp;
        }

        ft = funcd(at);    //Valor de la funcion ojetivo en el nuevo pto xi + at*pi
        dft = funcd.df(at);  //Valor de la derivada de la funcion objetivo en el nuevo pto xi + at*pi

        //Funcion auxiliar en "at"
        gt = ft - f0 - M*at*df0;
        dgt = dft - M*df0;
    }
}
#endif // LINESEARCHNUMTESTFUNC_H
