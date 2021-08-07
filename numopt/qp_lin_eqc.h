#ifndef QP_LIN_EQC_H
#define QP_LIN_EQC_H

#include "linalg/itermethods.h"
#include "linalg/qr.h"
#include "linalg/directsolvers.h"

/* Algoritmo de resolucion de problemas QP q(x) = 1/2*(x*G*x) + c*x
 * con restricciones de igualdad lineales Ax=b mediante Iterative Null
 * Space Methods Numerical Optimization 2ed 2006, J. Nocedal y S.J. Wright,
 * pp 16.1 - 16.3 */

class qp_lin_eqc
{
private:
    size_t index(const size_t &i, const size_t &j,const size_t &n)
    {
        size_t ind = (i <= j) ? (i*n + j - i*(i+1)/2) : (j*n + i - j*(j+1)/2);
        return ind;
    }
    //R se almacena en un vector de m*(m+1)/2 porque es upper triangular y asi se ahorra espacio
    void determine_Y_Z_R_QRDCMP(const MatDoub &Aeq, MatDoub &Y,MatDoub &Z,VecDoub &R){
        size_t m = Aeq.rows();
        size_t n = Aeq.cols();
        //QR DCMP
        qrdcmp qr(Aeq.transpose());
        qr.householder_qr();

        //EXTRAYENDO Y y Z
        Y = qr.Q(slice(0,n,1),slice(0,m,1));
        Z = qr.Q(slice(0,n,1),slice(m,n-m,1));

        //Almacenando qr.R en un el vector R
        for (size_t i = 0; i < m; ++i)
            for (size_t j = i; j < m; ++j)
                R[index(i,j,m)] = qr.R(i,j);
    }
public:
    qp_lin_eqc() = default;

    void qp_lin_eqc_iter_null_space_qrdcmp(const MatDoub &G, const VecDoub &c, const MatDoub &Aeq, const VecDoub &beq,VecDoub &x){
        size_t m = Aeq.rows();
        size_t n = Aeq.cols();
        VecDoub R(m*(m+1)/2); //Para almacenar R
        //determine_Y_Z_R_QRDCMP(Aeq,Y,Z,R);
        qrdcmp qr(Aeq.transpose());
        qr.householder_qr();
        for (size_t i = 0; i < m; ++i)
            for (size_t j = i; j < m; ++j)
                R[index(i,j,m)] = qr.R(i,j);
        qr.R.clear();
        MatDoub2DRef Y = qr.Q(slice(0,n,1),slice(0,m,1));
        MatDoub2DRef Z = qr.Q(slice(0,n,1),slice(m,n-m,1));

        VecDoub x_y(m);

        for (size_t k = 0; k < m; ++k){
            double sum = 0.0;
            for (size_t i = 0; i < k; ++i){
                sum += x_y[i]*R[index(i,k,m)];
            }
            x_y[k] = (beq[k] - sum)/R[index(k,k,m)];
        }
        R.clear(); //Ya no necesitamos mas a R;

        //TEst
        //VecDoub xy(2,0.0);
        //drtslv::LUdcmp LU(Aeq*Y);
        //LU.solve(beq,xy);

        VecDoub cz = Y*x_y;
        cz = G*cz;
        cz = cz*Z;
        cz += c*Z; //Z_t*c = c*Z

        VecDoub x_z(n-m,0.0);

        cg_redsys(Z,G,cz,x_z);

        x = Y*x_y;
        x += Z*x_z;
    }
};


#endif // QP_LIN_EQC_H
