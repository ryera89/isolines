#ifndef MATHMATRIX_H
#define MATHMATRIX_H

#include "mathvector.h"
#include "matrix_ref2d.h"
#include <cassert>
#include <algorithm>
#include <ostream>
#include <fstream>
#include <iomanip>

using std::vector;
using std::size_t;
using std::for_each;
using std::ostream;
using std::ifstream;
using std::ofstream;



template<typename T>
class matrix
{
public:
    explicit matrix():m_rows(0),m_cols(0){}
    explicit matrix(const size_t &m,const size_t &n):m_rows(m),m_cols(n),m_ptr(m*n){
        if (m==0 && n == 0){
            m_rows = 0;
            m_cols = 0;
        }
    }
    explicit matrix(const size_t &m,const size_t &n,const T &val):m_rows(m),m_cols(n),m_ptr(m*n,val){
        if (m==0 && n == 0){
            m_rows = 0;
            m_cols = 0;
        }
    }
    matrix(const matrix<T> &m) = default;
    matrix(matrix<T> &&m) = default;
    matrix<T>& operator =(const matrix<T> &m) = default;
    matrix<T>& operator =(matrix<T> &&m) = default;
    matrix(const matrix_ref2D<T> &mref):m_rows(mref.rows()),m_cols(mref.cols()),
        m_ptr(m_rows*m_cols){
        //m_ptr.shrink_to_fit();
        for (size_t i = 0; i < m_rows; ++i)
            for (size_t j = 0; j < m_cols; ++j)
                m_ptr[i*m_cols + j] = mref(i,j);
    }
    //WARNING ver que esto no se puede usar con el mismo la misma matrix eg Z = Z(slice,slice)
    matrix<T>& operator =(const matrix_ref2D<T> &mref){
           m_rows = mref.rows();
           m_cols = mref.cols();
           m_ptr.resize(m_rows*m_cols);
           for (size_t i = 0; i < m_rows; ++i)
               for (size_t j = 0; j < m_cols; ++j)
                   m_ptr[i*m_cols + j] = mref(i,j);

           return *this;
    }

    matrix_ref1D<T> row(const size_t &i){
        assert(i < m_rows);
        slice s(i*m_cols,m_cols,1);
        return matrix_ref1D<T>(s,m_ptr.data());
    }
    matrix_ref1D<const T> row(const size_t &i) const{
        assert(i < m_rows);
        slice s(i*m_cols,m_cols,1);
        return matrix_ref1D<const T>(s,data());
    }

    matrix_ref1D<T> column(const size_t &i){
        assert(i < m_cols);
        slice s(i,m_rows,m_cols);
        return matrix_ref1D<T>(s,m_ptr.data());
    }
    matrix_ref1D<const T> column(const size_t &i) const{
        assert(i < m_cols);
        slice s(i,m_rows,m_cols);
        return matrix_ref1D<const T>(s,data());
    }

    T& operator ()(const size_t &i, const size_t &j){return m_ptr[m_cols*i + j];}
    const T& operator ()(const size_t &i, const size_t &j) const{return m_ptr[m_cols*i + j];}

    matrix_ref2D<T> operator ()(const slice &s1,const slice &s2){

        if (s1.m_size == 0 || s2.m_size == 0 ) return matrix_ref2D<T>();

        size_t row_begin = s1.m_start;
        size_t row_end = row_begin + s1.m_size;
        size_t col_begin = s2.m_start;
        size_t col_end = col_begin + s2.m_size;

        assert(row_end <= m_rows && col_end <= m_cols);

        matrix_slice mslice;
        mslice.m_start = row_begin * m_cols + col_begin;
        mslice.m_size = s1.m_size * s2.m_size;
        mslice.m_rows = s1.m_size;
        mslice.m_cols = s2.m_size;
        mslice.m_stride = m_cols;

        return matrix_ref2D<T>(mslice,data());
    }
    matrix_ref2D<const T> operator ()(const slice &s1,const slice &s2) const{
        if (s1.m_size == 0 || s2.m_size == 0 ) return matrix_ref2D<const T>();

        size_t row_begin = s1.m_start;
        size_t row_end = row_begin + s1.m_size;
        size_t col_begin = s2.m_start;
        size_t col_end = col_begin + s2.m_size;

        assert(row_end <= m_rows && col_end <= m_cols);

        matrix_slice mslice;
        mslice.m_start = row_begin * m_cols + col_begin;
        mslice.m_size = s1.m_size * s2.m_size;
        mslice.m_rows = s1.m_size;
        mslice.m_cols = s2.m_size;
        mslice.m_stride = m_cols;

        return matrix_ref2D<const T>(mslice,data());
    }

    size_t rows() const{return m_rows;}
    size_t cols() const{return m_cols;}

    matrix<T> transpose() const{
        if (m_cols == m_rows){
            matrix<T> res(*this);
            for (size_t i = 0; i < m_rows; ++i)
                for (size_t j = i+1; j < m_cols; ++j)
                    std::swap(res(i,j),res(j,i));

            return res;
        }else{
           //TODO Implementer un algoritmo mas optimo
            matrix<T> res(m_cols,m_rows);
            for (size_t i = 0; i < m_rows; ++i)
                for (size_t j = 0; j < m_cols; ++j)
                    res(j,i) = m_ptr[i*m_cols + j];

            return res;
        }
    }

    //Operaciones Aritmeticas
    template<typename Scalar>
    matrix<T>& operator +=(const Scalar &val){
        if (val == 0) return *this;
        for_each(m_ptr.begin(),m_ptr.end(),[&val](T &elem){elem+=val;});
        return *this;
    }
    template<typename Scalar>
    matrix<T>& operator -=(const Scalar &val){
        if (val == 0) return *this;
        for_each(m_ptr.begin(),m_ptr.end(),[&val](T &elem){elem-=val;});
        return *this;
    }
    template<typename Scalar>
    matrix<T>& operator *=(const Scalar &val){
        if (val == 1) return *this;
        for_each(m_ptr.begin(),m_ptr.end(),[&val](T &elem){elem*=val;});
        return *this;
    }
    template<typename Scalar>
    matrix<T>& operator /=(const Scalar &val){
        if (val == 1) return *this;
        for_each(m_ptr.begin(),m_ptr.end(),[&val](T &elem){elem/=val;});
        return *this;
    }
    template<typename T2>
    matrix<T>& operator +=(const matrix<T2> &m1){
        assert(m_rows == m1.rows() && m_cols == m1.cols());
        auto it = m1.cbegin();
        for_each(m_ptr.begin(),m_ptr.end(),[&m1,&it](T &elem){elem += *it; ++it;});
        return *this;
    }
    template<typename T2>
    matrix<T>& operator -=(const matrix<T2> &m1){
        assert(m_rows == m1.rows() && m_cols == m1.cols());
        auto it = m1.cbegin();
        for_each(m_ptr.begin(),m_ptr.end(),[&m1,&it](T &elem){elem -= *it; ++it;});
        return *this;
    }
    template<typename T2>
    matrix<T>& operator +=(const matrix_ref2D<T2> &m1){
        assert(m_rows == m1.rows() && m_cols == m1.cols());
        for (size_t i = 0; i < m_rows; ++i)
            for (size_t j = 0; j < m_cols; ++j) m_ptr[i*m_cols+j] += m1(i,j);

        return *this;
    }
    template<typename T2>
    matrix<T>& operator -=(const matrix_ref2D<T2> &m1){
        assert(m_rows == m1.rows() && m_cols == m1.cols());
        for (size_t i = 0; i < m_rows; ++i)
            for (size_t j = 0; j < m_cols; ++j) m_ptr[i*m_cols+j] -= m1(i,j);
        return *this;
    }
    matrix<T>& swap_rows(const size_t &i, const size_t &j){
        assert(i < m_rows && j < m_rows);
        if (i == j) return *this;
        for (size_t k = 0; k < m_cols; ++k) std::swap(m_ptr[m_cols*i + k],m_ptr[m_cols*j + k]);
        return *this;
    }
    matrix<T>& prepend_row(const matrix_ref1D<T> &new_row){
        assert(new_row.size() == m_cols);
        vector<T> temp(m_cols);
        for (size_t i = 0; i < m_cols; ++i) temp[i] = new_row[i];

        m_ptr.insert(m_ptr.begin(),temp.cbegin(),temp.cend());

        ++m_rows;

        return *this;
    }

    matrix<T>& append_row(const matrix_ref1D<T> &new_row){
        assert(new_row.size() == m_cols);
        //++m_rows;
        m_ptr.resize(size()+m_cols);

        for (size_t j = 0; j < m_cols; ++j){
            m_ptr[m_cols*(m_rows) + j] = new_row[j];
        }
        ++m_rows;
        return *this;
    }
    matrix<T>& prepend_column(const matrix_ref1D<T> &new_col){
        assert(new_col.size() == m_rows);
        size_t n = m_cols+1;
        matrix<T> temp(m_rows,n);

        for (size_t i = 0; i < m_rows; ++i){
            temp(i,0) = new_col[i];
            for (size_t j = 0; j < m_cols; ++j){
                temp(i,j+1) = m_ptr[i*m_cols+j];
            }

        }
        *this = temp;

        return *this;
    }

    matrix<T>& append_column(const matrix_ref1D<T> &new_col){
        assert(new_col.size() == m_rows);
        size_t n = m_cols+1;
        matrix<T> temp(m_rows,n);

        for (size_t i = 0; i < m_rows; ++i){
            for (size_t j = 0; j < m_cols; ++j){
                temp(i,j) = m_ptr[i*m_cols+j];
            }
            temp(i,m_cols) = new_col[i];
        }
        *this = temp;

        return *this;
    }
    matrix<T>& remove_column(const size_t &col){
        assert(col < m_cols);
        size_t n = m_cols-1;
        matrix<T> temp(m_rows,n);

        for (size_t i = 0; i < m_rows; ++i){
            for (size_t j = 0; j < col; ++j){
                temp(i,j) = m_ptr[i*m_cols+j];
            }
            for (size_t j = col+1; j < m_cols; ++j){
                temp(i,j-1) = m_ptr[i*m_cols+j];
            }
        }
        *this = temp;

        return *this;
    }

    void resize(const size_t &nrow,const size_t &ncols){
        if (ncols == 0 || nrow == 0){
            m_rows = 0;
            m_cols = 0;
        }else{
            m_rows = nrow;
            m_cols = ncols;
        }

        m_ptr.resize(nrow*ncols);
    }
    size_t size() const{return m_ptr.size();}
    T* data() {return m_ptr.data();}
    const T* data() const {return m_ptr.data();}


    typename vector<T>::const_iterator cbegin() const{return m_ptr.cbegin();}
    typename vector<T>::const_iterator cend() const{return m_ptr.cend();}

    typename vector<T>::iterator begin() {return m_ptr.begin();}
    typename vector<T>::iterator end() {return m_ptr.end();}

    void setToZero(){for (T &v : m_ptr) v = 0;}
    void clear(){
        m_rows = 0;
        m_cols = 0;
        m_ptr.clear();
    }
private:
    size_t m_rows;
    size_t m_cols;
    std::vector<T> m_ptr;

};
template<typename T>
ostream& operator <<(ostream& os,const matrix<T> &m){
    std::ios_base::fmtflags ff = std::ios::scientific;
    ff |= std::ios::showpos;
    os.setf(ff);
    for (size_t i = 0; i < m.rows(); ++i){
        for (size_t j = 0; j < m.cols(); ++j)
            os << m(i,j) << '\t' ;
        os << std::endl;
    }
    os.unsetf(ff);
    return os;
}
template<typename T>
ifstream& operator >> (ifstream &in,matrix<T> &m){
    if (in.is_open()){
        for (size_t i = 0; i < m.rows(); ++i){
            for (size_t j = 0; j < m.cols(); ++j)
            in >> m(i,j);
        }
    }
    return in;
}
template<typename T>
ofstream& operator << (ofstream &ofs,const matrix<T> &m){
    if (ofs.is_open()){
        std::ios_base::fmtflags ff = std::ios::scientific;
        ff |= std::ios::showpos;
        ofs.setf(ff);
        for (size_t i = 0; i < m.rows(); ++i){
            for (size_t j = 0; j < m.cols(); ++j)
                ofs << m(i,j) << '\t' ;
            ofs << std::endl;
        }
        ofs.unsetf(ff);
    }
    return ofs;

}
//Matrix Aritmetic Ops ----------------------------------------------------
template<typename T>
inline matrix<T> sum(const matrix<T> &m1,const matrix<T> &m2){
    matrix<T> res(m1);
    res+=m2;
    return res;
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline matrix<RT> sum(const matrix<T1> &m1, const matrix<T2> &m2){
    assert(m1.rows() == m2.rows() && m1.cols() == m2.cols());
    matrix<RT> res(m1.rows(),m1.cols());
    auto it1 = m1.cbegin();
    auto it2 = m2.cbegin();
    for_each(res.begin(),res.end(),[&it1,&it2](RT &elem){elem = *it1 + *it2; ++it1; ++it2;});
    return res;
}
template<typename T>
inline matrix<T> sub(const matrix<T> &m1,const matrix<T> &m2){
    matrix<T> res(m1);
    res-=m2;
    return res;
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline matrix<RT> sub(const matrix<T1> &m1, const matrix<T2> &m2){
    assert(m1.rows() == m2.rows() && m1.cols() == m2.cols());
    matrix<RT> res(m1.rows(),m1.cols());
    auto it1 = m1.cbegin();
    auto it2 = m2.cbegin();
    for_each(res.begin(),res.end(),[&it1,&it2](RT &elem){elem = *it1 - *it2; ++it1; ++it2;});
    return res;
}
template<typename T>
inline matrix<T> mult(const T &val, const matrix<T> &m){
    if (val == 0){
        matrix<T> res(m.rows(),m.cols(),0.0);
        return res;
    }
    matrix<T> res(m);
    for_each(res.begin(),res.end(),[&val](T &elem){elem = val*elem;});
    return res;
}
template<typename T,typename S,typename RT = Common_type<T,S>>
inline matrix<RT> mult(const S &val, const matrix<T> &m){
    if (val == 0){
        matrix<RT> res(m.rows(),m.cols(),0.0);
        return res;
    }
    matrix<RT> res(m.rows(),m.cols());
    auto it = m.cbegin();
    for_each(res.begin(),res.end(),[&val,&it](RT &elem){elem = val*(*it); ++it;});
    return res;
}
template<typename T>
inline matrix<T> mult(const matrix<T> &m,const T &val){
    if (val == 0){
        matrix<T> res(m.rows(),m.cols(),0.0);
        return res;
    }
    matrix<T> res(m);
    res*=val;
    return res;
}
template<typename T,typename S,typename RT = Common_type<T,S>>
inline matrix<RT> mult(const matrix<T> &m,const S &val){
    if (val == 0){
        matrix<RT> res(m.rows(),m.cols(),0.0);
        return res;
    }
    matrix<RT> res(m.rows(),m.cols());
    auto it = m.cbegin();
    for_each(res.begin(),res.end(),[&val,&it](RT &elem){elem = (*it)*val; ++it;});
    return res;
}
template<typename T>
inline matrix<T> div(const matrix<T> &m,const T &val){
    matrix<T> res(m);
    res/=val;
    return res;
}
template<typename T,typename S,typename RT = Common_type<T,S>>
inline matrix<RT> div(const matrix<T> &m,const S &val){
    matrix<RT> res(m.rows(),m.cols());
    auto it = m.cbegin();
    for_each(res.begin(),res.end(),[&val,&it](RT &elem){elem = (*it)/val; ++it;});
    return res;
}
template<typename T>
inline mathvector<T> mult(const mathvector<T> &v,const matrix<T> &m){
    if (v.size() == 0 || m.rows() == 0 || m.cols() == 0){
        return mathvector<T>();
    }
    assert(v.size() == m.rows());
    mathvector<T> res(m.cols(),0.0);
    for (size_t i = 0; i < m.rows(); ++i) mult(v[i],m.row(i),res);
    return res;
}
template<typename T>
inline mathvector<T> mult(const matrix<T> &m, const mathvector<T> &vec){
    assert(m.cols() == vec.size());
    mathvector<T> res(m.rows());
    for (size_t i = 0; i < m.rows(); ++i){
        double temp = 0.0;
        for (size_t j = 0; j < m.cols(); ++j)
            temp += m(i,j)*vec[j];
        res[i] = temp;
    }
    return res;
}
template<typename T1,typename T2, typename RT = Common_type<T1,T2>>
inline matrix<RT> mult(const matrix<T1> &m1, const matrix<T2> &m2){
    if (m1.cols() == 0 || m1.rows() == 0 || m2.cols() == 0 || m2.rows() == 0){
        return matrix<RT>();
    }
    assert(m1.cols() == m2.rows());
    matrix<RT> res(m1.rows(),m2.cols(),0.0);
    for (size_t i = 0; i < m1.rows(); ++i){
        matrix_ref1D<RT> ref = res.row(i);
        for (size_t j = 0; j < m2.rows(); ++j){
            mult(m1(i,j),m2.row(j),ref);
        }
    }
    return res;
}
template<typename T>
inline matrix<T> operator +(const matrix<T> &m1,const matrix<T> &m2){
    return sum(m1,m2);
}
template<typename T1, typename T2,typename RT = Common_type<T1,T2>>
inline matrix<RT> operator +(const matrix<T1> &m1,const matrix<T2> &m2){
    return sum(m1,m2);
}
template<typename T>
inline matrix<T> operator -(const matrix<T> &m1,const matrix<T> &m2){
    return sub(m1,m2);
}
template<typename T1, typename T2,typename RT = Common_type<T1,T2>>
inline matrix<RT> operator -(const matrix<T1> &m1,const matrix<T2> &m2){
    return sub(m1,m2);
}
template<typename T>
inline matrix<T> operator *(const T &val,const matrix<T> &m){
    return mult(val,m);
}
template<typename T1, typename T2,typename RT = Common_type<T1,T2>>
inline matrix<RT> operator *(const T1 &val,const matrix<T2> &m){
    return mult(val,m);
}
template<typename T>
inline matrix<T> operator *(const matrix<T> &m,const T &val){
    return mult(m,val);
}
template<typename T1, typename T2,typename RT = Common_type<T1,T2>>
inline matrix<RT> operator *(const matrix<T1> &m,const T2 &val){
    return mult(m,val);
}
template<typename T>
inline matrix<T> operator /(const matrix<T> &m,const T &val){
    return div(m,val);
}
template<typename T1, typename T2,typename RT = Common_type<T1,T2>>
inline matrix<RT> operator /(const matrix<T1> &m,const T2 &val){
    return div(m,val);
}
template<typename T>
inline mathvector<T> operator *(const mathvector<T> &v, const matrix<T> &m){
    return mult(v,m);
}
template<typename T>
inline mathvector<T> operator *(const matrix<T> &m,const mathvector<T> &v){
    return mult(m,v);
}
template<typename T1, typename T2,typename RT = Common_type<T1,T2>>
inline matrix<RT> operator *(const matrix<T1> &m1,const matrix<T2> &m2){
    return mult(m1,m2);
}
//--------------------------------------------------------------------------------
inline matrix<double> prod_t(const mathvector<double> &vec1,const mathvector<double> &vec2_t){
    matrix<double> res(vec1.size(),vec2_t.size());
    for (size_t i = 0; i < res.rows(); ++i) res.row(i) = vec1[i]*vec2_t;
    return res;
}
//---------------------------------------------------------------------------------
//matrixref_2D Aritmetic ops-------------------------------------------------------
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline matrix<RT> sum(const matrix_ref2D<T1> &m1, const matrix_ref2D<T2> &m2){
    assert(m1.rows() == m2.rows() && m1.cols() == m2.cols());
    matrix<RT> res(m1.rows(),m1.cols());
    for (size_t i = 0; i < res.rows(); ++i)
        for(size_t j = 0; j < res.cols(); ++j)
            res(i,j) = m1(i,j) + m2(i,j);
    return res;
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline matrix<RT> sub(const matrix_ref2D<T1> &m1, const matrix_ref2D<T2> &m2){
    assert(m1.rows() == m2.rows() && m1.cols() == m2.cols());
    matrix<RT> res(m1.rows(),m1.cols());
    for (size_t i = 0; i < res.rows(); ++i)
        for(size_t j = 0; j < res.cols(); ++j)
            res(i,j) = m1(i,j) - m2(i,j);
    return res;
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline matrix<RT> operator +(const matrix_ref2D<T1> &m1, const matrix_ref2D<T2> &m2){
    return sum(m1,m2);
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline matrix<RT> operator -(const matrix_ref2D<T1> &m1, const matrix_ref2D<T2> &m2){
    return sub(m1,m2);
}
typedef matrix<double> MatDoub;


template<typename T>
inline matrix_ref2D<T>& operator +=(matrix_ref2D<T> &ref, const matrix<T> &m1){
    assert(ref.rows() == m1.rows() && ref.cols() == m1.cols());
    for (size_t i = 0; i < ref.rows(); ++i)
        for (size_t j = 0; j < ref.cols(); ++j)
            ref(i,j) += m1(i,j);
    return ref;
}
template<typename T>
inline matrix_ref2D<T>& operator -=(matrix_ref2D<T> &ref, const matrix<T> &m1){
    assert(ref.rows() == m1.rows() && ref.cols() == m1.cols());
    for (size_t i = 0; i < ref.rows(); ++i)
        for (size_t j = 0; j < ref.cols(); ++j)
            ref(i,j) -= m1(i,j);
    return ref;
}
//TODO mejorar esto
inline MatDoub IdentityMatrix(const size_t &n){
    MatDoub I(n,n,0.0);
    for (size_t i = 0; i < n; ++i) I(i,i) = 1.0;
    return I;
}
//matrix matrix_ref1D ops
template<typename T>
inline mathvector<T> operator *(const matrix_ref1D<T> &v, const matrix<T> &m){
    if (v.size() == 0 || m.rows() == 0 || m.cols() == 0){
        return mathvector<T>();
    }
    assert(v.size() == m.rows());
    mathvector<T> res(m.cols(),0.0);
    for (size_t i = 0; i < m.rows(); ++i) mult(v[i],m.row(i),res);
    return res;
}
template<typename T>
inline mathvector<T> operator *(const matrix<T> &m,const matrix_ref1D<T> &v){
    assert(m.cols() == v.size());
    mathvector<T> res(m.rows());
    for (size_t i = 0; i < m.rows(); ++i){
        double temp = 0.0;
        for (size_t j = 0; j < m.cols(); ++j)
            temp += m(i,j)*v[j];
        res[i] = temp;
    }
    return res;
}

#endif // MATHMATRIX_H
