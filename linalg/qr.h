#ifndef QR_H
#define QR_H


#include <set>
#include "householdervec.h"
#include "mathmatrix.h"
#include "givens_rotations.h"

using std::set;

/* Golub, G.H., and Van Loan, C.F. 1996,
 * Matrix Computations, 3rd ed. (Baltimore: Johns Hopkins University Press)
 * Algorithm 5.2.1. Factorizacion QR mediante Householder, dado una matriz A mxn
 * m >= n encuentra los H1,H2....Hn tal que Q = H1*H2...*Hn y Transp(Q)*A = R la matriz
 * R se almacena en la parte triangular superior de A y en la parte inferior de A se almacenan
 * los vectores de householder vj */
struct qrdcmp{

    MatDoub R;
    MatDoub Q;

    qrdcmp(const MatDoub &AA):R(AA){assert(AA.rows() >= AA.cols());}
    qrdcmp(MatDoub &&AA):R(AA){assert(AA.rows() >= AA.cols());}
    //Especializado para optimizacion Metodo de Active Set
    qrdcmp(const MatDoub &Aeq,const MatDoub &Ainq, const set<size_t> &W,const size_t &n):
        R(n,Aeq.rows() + W.size()){
        for (size_t i = 0; i < R.rows(); ++i){
            for (size_t j = 0; j < Aeq.rows(); ++j){
                R(i,j) = Aeq(j,i);
            }
            size_t k = Aeq.rows();
            for (auto iter = W.begin(); iter != W.end(); ++iter){
                R(i,k) = Ainq(*iter,i);
                ++k;
            }
        }
    }

    qrdcmp(const qrdcmp &other) = delete;
    qrdcmp(qrdcmp &&other) = delete;
    qrdcmp& operator =(const qrdcmp &other) = delete;
    qrdcmp& operator =(qrdcmp &&other) = delete;

    //TODO despues arreglar esta funcion para que el usuario decida si quiere calcular Q
    void householder_qr(){
        size_t m = R.rows();
        size_t n = R.cols();

        vector<double> betj(n);
        betj.shrink_to_fit();
        householdervec hhv;
        for (size_t j = 0; j < R.cols(); ++j){
            hhv(R(slice(j,m-j,1),slice(j,1,1)));
            betj[j] = hhv.bet;
            MatDoub2DRef mref_jm_jn(R(slice(j,m-j,1),slice(j,n-j,1)));
            VecDoub w = hhv.v*mref_jm_jn;
            w*=hhv.bet;
            mref_jm_jn-=prod_t(hhv.v,w);
            if (j < m-1){
                for (size_t i = j + 1; i < m; ++i) R(i,j) = hhv.v[i - j];
            }
        }
        /* Calculo de Q por acumulacion backguard algorithm 5.1.5 Golub, G.H.,
         * and Van Loan, C.F. 1996, Matrix Computations, 3rd ed */
        Q = IdentityMatrix(m);

        for (long int j = n-1; j >= 0;){
            MatDoub2DRef mref_jm_jm = Q(slice(j,m-j,1),slice(j,m-j,1));
            VecDoub w = hhv.v*mref_jm_jm;
            w*=betj[j];
            mref_jm_jm-=prod_t(hhv.v,w);
            hhv.v.push_back(0.0); //aumento en 1 la longitud de v
            --j; //disminuyo j
            for (size_t i = 1; i < hhv.v.size(); ++i) hhv.v[i] = R(m-hhv.v.size()+i,j);
            hhv.v[0] = 1.0;
        }

    }
    void givens_qr(){
        size_t m = R.rows();
        size_t n = R.cols();

        for (size_t i = m-1; i >= 1; --i){
            for (size_t j = 0; i < std::min(i,n); ++j){
                GivensRotation givens_rot(R(i-1,j),R(i,j));
                givens_rot.apply_givens_to_rows(i-1,i,j,n,R);
            }
        }
    }
};


#endif // QR_H
