#ifndef LUDCMP_H
#define LUDCMP_H

#include "mathmatrix.h"

class LUdcmp{
private:
    size_t n;
    MatDoub LU;
    std::vector<size_t> p_ind_v;
    int d;
public:
    LUdcmp(const MatDoub &A);
    /*LUdcmp(const MatDoub &Aeq,const MatDoub &Ainq,const set<size_t> &W)
        :n(Aeq.rows() + W.size()),LU(n),p_ind_v(n),d(1.0){
        for (size_t i = 0; i < n; ++i){
            for (size_t j = 0; j < Aeq.rows(); ++j){
                LU(i,j) = Aeq(j,i);
            }
            size_t k = Aeq.rows();
            for (auto iter = W.begin(); iter != W.end() ; ++iter){
                LU(i,k) = Ainq(*iter,i);
                ++k;
            }
        }
        p_ind_v.shrink_to_fit();

        double big,temp;
        size_t imax = 0;

        for (size_t k = 0; k < n; ++k){   //Lazo principal
            big = 0.0;
            for (size_t i = k; i < n; ++i){ //Busqueda del pivote
                temp = std::abs(LU(i,k));
                if (temp > big){
                    big = temp;
                    imax = i;     //guardo la fila del pivote
                }
            }
            if (k != imax){
                LU.swap_rows(k,imax); //Intercambia las filas
                d = -d; //cambia la paridad del discriminante
            }
            p_ind_v[k] = imax;   //guardo la fila del pivote para despues usar despues en solve()

            if (LU(k,k) == 0.0) throw("LUdcmp::LUdcmp Sistema singular o  mal condicionado");
            for (size_t i = k+1; i < n; ++i){
                temp = LU(i,k) /= LU(k,k);
                for (size_t j = k+1; j < n; ++j){
                    LU(i,j) -= temp*LU(k,j);
                }
            }
        }
    }*/

    void solve(const VecDoub &b,VecDoub &x);

    /*void solve_lagrange_mult(const VecDoub &beq,const VecDoub &binq, const set<size_t> &W,VecDoub &x){
        for (size_t i = 0; i < beq.size(); ++i) x[i] = beq[i];
        size_t k = beq.size();
        for (auto iter = W.begin(); iter != W.end(); ++iter){
            x[k] = binq[*iter];
            ++k;
        }
        size_t ip, ii = 0;
        double sum;

        for (size_t i = 0; i < n; ++i){ //Ly = b
            ip = p_ind_v[i];    //Registro de permutaciones de columnas
            sum = x[ip];
            x[ip] = x[i];
            if (ii != 0)
                for (size_t j = ii - 1; j < i; ++j) sum -= LU(i,j)*x[j];
            else if (sum != 0.0) ii = i + 1; //Un elemento diferente de cero se a encontrado y hay que entrar el loop de arriba

            x[i] = sum;
        }
        for (long int i = n-1; i >= 0; --i){ //Backsustitution Ux=y
            sum = x[i];
            for (size_t j = i + 1; j < n; ++j) sum -= LU(i,j)*x[j];
            x[i] = sum/LU(i,i);
        }


    }*/
    void solve(const MatDoub &B, MatDoub &X);

    double det();
};

inline LUdcmp::LUdcmp(const MatDoub &A):n(A.rows()),LU(A),p_ind_v(n),d(1.0){
    p_ind_v.shrink_to_fit(); //ahorrar memoria, no se va a agregar ningun elemento
    double big,temp;
    size_t imax = 0;

    for (size_t k = 0; k < n; ++k){   //Lazo principal
        big = 0.0;
        for (size_t i = k; i < n; ++i){ //Busqueda del pivote
            temp = std::abs(LU(i,k));
            if (temp > big){
                big = temp;
                imax = i;     //guardo la fila del pivote
            }
        }
        if (k != imax){
            LU.swap_rows(k,imax); //Intercambia las filas
            d = -d; //cambia la paridad del discriminante
        }
        p_ind_v[k] = imax;   //guardo la fila del pivote para despues usar despues en solve()

        big = LU(k,k);
        if (big == 0.0) throw("LUdcmp::LUdcmp Sistema singular o  mal condicionado");
        for (size_t i = k+1; i < n; ++i){
            temp = LU(i,k) /= big;
            for (size_t j = k+1; j < n; ++j){
                LU(i,j) -= temp*LU(k,j);
            }
        }
    }
}

inline void LUdcmp::solve(const VecDoub &b, VecDoub &x)
{
    assert(b.size() == n && x.size() == n);

    std::copy(b.begin(),b.end(),x.begin());

    size_t ip, ii = 0;
    double sum;

    for (size_t i = 0; i < n; ++i){ //Ly = b
        ip = p_ind_v[i];    //Registro de permutaciones de columnas
        sum = x[ip];
        x[ip] = x[i];
        if (ii != 0)
            for (size_t j = ii - 1; j < i; ++j) sum -= LU(i,j)*x[j];
        else if (sum != 0.0) ii = i + 1; //Un elemento diferente de cero se a encontrado y hay que entrar el loop de arriba

        x[i] = sum;
    }
    for (long int i = n-1; i >= 0; --i){ //Backsustitution Ux=y
        sum = x[i];
        for (size_t j = i + 1; j < n; ++j) sum -= LU(i,j)*x[j];
        x[i] = sum/LU(i,i);
    }

}

inline void LUdcmp::solve(const MatDoub &B, MatDoub &X)
{
    assert(B.rows() == n && X.rows() == n && B.cols() == X.cols());

    size_t m = B.cols();

    VecDoub xx(n);
    for (size_t j = 0; j < m; ++j){
        for (size_t i = 0; i < n; ++i) xx[i] = B(i,j);
        solve(xx,xx);
        for (size_t i = 0;i < n; ++i) X(i,j) = xx[i];
    }
}

inline double LUdcmp::det()
{
    double dd = d;
    for (size_t i = 0; i < n; ++i) dd*=LU(i,i);

    return dd;
}
#endif // LUDCMP_H
