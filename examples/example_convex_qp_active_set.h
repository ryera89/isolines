#ifndef EXAMPLE_CONVEX_QP_ACTIVE_SET_H
#define EXAMPLE_CONVEX_QP_ACTIVE_SET_H
/*Ejemplos 16.3 pp.452 y 16.4 pp.475 J. Nocedal "Numerical Optimization" 2ed
 * Aplicacion del algoritmo Active-Set */

#include "numopt/qp_lin_convex_act_set.h"
#include "numopt/qp_lin_eqc.h"
#include "numopt/qp_gpm_bc.h"

struct example16_3
{
    MatDoub G;
    VecDoub c;

    MatDoub Aeq;
    VecDoub beq;

    VecDoub x;
    example16_3():G(3,3),c(3),Aeq(2,3),beq(2),x(3) {
        G(0,0) = 6.0;
        G(0,1) = 2.0;
        G(0,2) = 1.0;

        G(1,0) = 2.0;
        G(1,1) = 5.0;
        G(1,2) = 2.0;

        G(2,0) = 1.0;
        G(2,1) = 2.0;
        G(2,2) = 4.0;

        c[0] = -8;
        c[1] = -3;
        c[2] = -3;

        Aeq(0,0) = 1.0;
        Aeq(0,1) = 0.0;
        Aeq(0,2) = 1.0;

        Aeq(1,0) = 0.0;
        Aeq(1,1) = 1.0;
        Aeq(1,2) = 1.0;

        beq[0] = 3;
        beq[1] = 0;
    }
    VecDoub& solve(){
        qp_lin_eqc qp;
        qp.qp_lin_eqc_iter_null_space_qrdcmp(G,c,Aeq,beq,x);
        return x;
    }
};

struct example16_4
{
    MatDoub G;
    VecDoub c;

    MatDoub Ainq;
    VecDoub binq;

    VecDoub x;


    example16_4():G(2,2,0.0),c(2),Ainq(5,2),binq(5,0.0),x(2) {
        G(0,0) = 2;
        G(1,1) = 2;

        Ainq(0,0) = 1;
        Ainq(0,1) = -2;
        Ainq(1,0) = -1;
        Ainq(1,1) = -2;
        Ainq(2,0) = -1;
        Ainq(2,1) = 2;
        Ainq(3,0) = 1;
        Ainq(3,1) = 0;
        Ainq(4,0) = 0;
        Ainq(4,1) = 1;

        binq[0] = -2;
        binq[1] = -6;
        binq[2] = -2;

        c[0] = -2;
        c[1] = -5;

        x(0) = 2;
        x(1) = 0;
    }

    VecDoub& solve(){
        qp_lin_convex_act_set qp_actset;
        qp_actset.qp_lin_convex_act_set_qrdcmp(G,c,MatDoub(),VecDoub(),Ainq,binq,x);
        return x;
    }
};
//Matlab example
struct matlab_example{
    MatDoub G;
    VecDoub c;
    VecDoub lb;
    VecDoub ub;
    VecDoub x;

    matlab_example():G(900,900),c(900),lb(900),ub(900,numeric_limits<double>::infinity()),x(900,0.5){
        ifstream in("example_qp_bc/H.txt");
        in >> G;
        in.close();
        in.open("example_qp_bc/C.txt");
        in >> c;
        in.close();
        in.open("example_qp_bc/LB.txt");
        in >> lb;
        in.close();
    }
    VecDoub& solve(){
        qp_gpm_bc qp_bc;
        qp_bc.minimize(G,c,lb,ub,x);
        return x;

    }
};
#endif // EXAMPLE_CONVEX_QP_ACTIVE_SET_H
