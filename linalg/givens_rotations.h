#ifndef GIVENS_ROTATIONS_H
#define GIVENS_ROTATIONS_H

#include <cmath>
#include "mathmatrix.h"


using std::abs;
using std::sqrt;

class GivensRotation{
private:
    double c;
    double s;

public:
    GivensRotation(const double &a, const double &b){
        if (b == 0.0){
            c = 1.0;
            s = 0.0;
        }else{
            if (abs(b) > abs(a)){
                double r = -a/b;
                s = 1.0/sqrt(1.0+r*r);
                c = s*r;
            }else{
                double r = -b/a;
                c = 1.0/sqrt(1.0+r*r);
                s = c*r;
            }
        }
    }
    void givens(const double &a, const double &b){
        if (b == 0.0){
            c = 1.0;
            s = 0.0;
        }else{
            if (abs(b) > abs(a)){
                double r = -a/b;
                s = 1.0/sqrt(1.0+r*r);
                c = s*r;
            }else{
                double r = -b/a;
                c = 1.0/sqrt(1.0+r*r);
                s = c*r;
            }
        }
    }
    void apply_givens_to_mref(const size_t &beg,const size_t &end,
                              MatDoub1DRef &mref1,MatDoub1DRef &mref2){
        assert(mref1.size() == mref2.size());
        double r1;
        double r2;
        for (size_t j = beg; j < end; ++j){
            r1 = mref1(j);
            r2 = mref2(j);
            mref1(j) = c*r1 - s*r2;
            mref2(j) = s*r1 + c*r2;
        }
    }
    void apply_givens_to_vec(const size_t &beg,const size_t &end,
                              VecDoub &vec1,VecDoub &vec2){
        assert(vec1.size() == vec2.size());
        double r1;
        double r2;
        for (size_t j = beg; j < end; ++j){
            r1 = vec1[j];
            r2 = vec2[j];
            vec1[j] = c*r1 - s*r2;
            vec2[j] = s*r1 + c*r2;
        }
    }
    VecDoub& apply_givens_to_vector(const size_t &i, const size_t &k, VecDoub &vec){
        assert(i < k);
        double r1;
        double r2;

        r1 = vec[i];
        r2 = vec[k];
        vec[i] = c*r1 - s*r2;
        vec[k] = s*r1 + c*r2;

        return vec;
    }
    MatDoub& apply_givens_to_rows(const size_t &i, const size_t &k, MatDoub &m){
        assert(i < k);
        double r1;
        double r2;
        for (size_t j = 0; j < m.cols(); ++j){
            r1 = m(i,j);
            r2 = m(k,j);
            m(i,j) = c*r1 - s*r2;
            m(k,j) = s*r1 + c*r2;
        }
        return m;
    }
    MatDoub& apply_givens_to_rows(const size_t &i, const size_t &k,const size_t &cols_beg,
                                  const size_t &cols_end, MatDoub &m){
        assert(i < k);
        double r1;
        double r2;
        for (size_t j = cols_beg; j < cols_end; ++j){
            r1 = m(i,j);
            r2 = m(k,j);
            m(i,j) = c*r1 - s*r2;
            m(k,j) = s*r1 + c*r2;
        }
        return m;
    }
    MatDoub& apply_givens_to_cols(const size_t &i, const size_t &k, MatDoub &m){
        assert(i < k);
        double r1;
        double r2;
        for (size_t j = 0; j < m.rows(); ++j){
            r1 = m(j,i);
            r2 = m(j,k);
            m(j,i) = c*r1 - s*r2;
            m(j,k) = s*r1 + c*r2;
        }
        return m;
    }
    MatDoub& apply_givens_to_cols(const size_t &i, const size_t &k,const size_t &rows_beg,
                                  const size_t &rows_end, MatDoub &m){
        assert(i < k);
        double r1;
        double r2;
        for (size_t j = rows_beg; j < rows_end; ++j){
            r1 = m(j,i);
            r2 = m(j,k);
            m(j,i) = c*r1 - s*r2;
            m(j,k) = s*r1 + c*r2;
        }
        return m;
    }
    MatDoub2DRef& apply_givens_to_rows(const size_t &i, const size_t &k, MatDoub2DRef &m){
        assert(i < k);
        double r1;
        double r2;
        for (size_t j = 0; j < m.cols(); ++j){
            r1 = m(i,j);
            r2 = m(k,j);
            m(i,j) = c*r1 - s*r2;
            m(k,j) = s*r1 + c*r2;
        }
        return m;
    }
    MatDoub2DRef& apply_givens_to_rows(const size_t &i, const size_t &k,const size_t &cols_beg,
                                  const size_t &cols_end, MatDoub2DRef &m){
        assert(i < k);
        double r1;
        double r2;
        for (size_t j = cols_beg; j < cols_end; ++j){
            r1 = m(i,j);
            r2 = m(k,j);
            m(i,j) = c*r1 - s*r2;
            m(k,j) = s*r1 + c*r2;
        }
        return m;
    }
    MatDoub2DRef& apply_givens_to_cols(const size_t &i, const size_t &k, MatDoub2DRef &m){
        assert(i < k);
        double r1;
        double r2;
        for (size_t j = 0; j < m.rows(); ++j){
            r1 = m(j,i);
            r2 = m(j,k);
            m(j,i) = c*r1 - s*r2;
            m(j,k) = s*r1 + c*r2;
        }
        return m;
    }
    MatDoub2DRef& apply_givens_to_cols(const size_t &i, const size_t &k,const size_t &rows_beg,
                                  const size_t &rows_end, MatDoub2DRef &m){
        assert(i < k);
        double r1;
        double r2;
        for (size_t j = rows_beg; j < rows_end; ++j){
            r1 = m(j,i);
            r2 = m(j,k);
            m(j,i) = c*r1 - s*r2;
            m(j,k) = s*r1 + c*r2;
        }
        return m;
    }
    //Construye la matriz de givens con los valores de c y s calculados
    MatDoub transpose_givens_matrix(const size_t &i, const size_t &k, const size_t &n){
        assert(i < k);
        MatDoub res(IdentityMatrix(n));

        res(i,i) = c;
        res(i,k) = -s;
        res(k,i) = s;
        res(k,k) = c;

        return res;
    }
    double C() const{
        return c;
    }
    double S() const{
        return s;
    }
};

#endif // GIVENS_ROTATIONS_H
