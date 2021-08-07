#ifndef CHOLESKY_H
#define CHOLESKY_H

#include "mathmatrix.h"
//TODO revisar cholesky para posible optimizacion
class Cholesky
{
public:
    Cholesky(const MatDoub &A);
    ~Cholesky() = default;

    void solve(const VecDoub &b,VecDoub &x);

private:
    size_t n;  //square matrix
    VecDoub m_L;


    size_t index(const size_t &i, const size_t &j);
};

Cholesky::Cholesky(const MatDoub &A):n(A.rows()),m_L(n*(n+1)/2)  //Se raaliza la descomposicion de la matriz a A = L*trans(L)
{
    assert(A.rows() == A.cols());

    double sum = 0.0;

    for (size_t i = 0; i < n; ++i){
        for (size_t j = i; j < n; ++j){
            sum = A(i,j);
            for (size_t k = 0; k < i; ++k) sum -= m_L[index(i,k)]*m_L[index(j,k)];

            if (i == j){
                assert(sum > 0.0); //Matriz definida positiva
                m_L[index(i,i)] = std::sqrt(sum);
            }else{
                m_L[index(j,i)] = sum/m_L[index(i,i)];
            }
        }
    }
}

inline void Cholesky::solve(const VecDoub &b, VecDoub &x)
{
    double sum = 0.0;

    assert(b.size() == n && x.size() == n);

    for (size_t i = 0; i < n; ++i){  //Resolviendo Ly = b L:lower triangular
        sum = b[i];
        for (size_t j = 0; j < i; ++j) sum -= m_L[index(i,j)]*x[j];

        x[i] = sum/m_L[index(i,i)];
    }

    for (long int i = n-1; i >= 0; --i){
        sum = x[i];
        for (size_t j = i + 1; j < n; ++j) sum -= m_L[index(i,j)]*x[j];

        x[i] = sum/m_L[index(i,i)];
    }

}


inline size_t Cholesky::index(const size_t &i, const size_t &j)
{
    size_t ind = (i <= j) ? (i*n + j - i*(i+1)/2) : (j*n + i - j*(j+1)/2);

    return ind;
}
#endif // CHOLESKY_H
