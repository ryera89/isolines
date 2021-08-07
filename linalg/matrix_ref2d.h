#ifndef MATRIX_REF2D_H
#define MATRIX_REF2D_H

#include "matrix_ref1d.h"
#include "mathvector.h"

struct matrix_slice{
    matrix_slice():m_start(0),m_size(0),m_rows(0),m_cols(0),m_stride(0){}
    //Contructor del slice a partir de las filas y las columnas
    matrix_slice(const size_t &rows,const size_t &cols):m_start(0),m_size(rows*cols),
        m_rows(rows),m_cols(cols),m_stride(cols){}
    //Constructor del slice a partir del start y extents
    matrix_slice(const size_t &start,const size_t &rows,const size_t &cols):m_start(start),
        m_size(rows*cols),m_rows(rows),m_cols(cols),m_stride(m_cols){}
    //Se utiliza otro stride
    matrix_slice(const size_t &start,const size_t &rows,const size_t &cols,const size_t &stride):
        m_start(start),m_size(rows*cols),m_rows(rows),m_cols(cols),m_stride(stride){}

    size_t operator()(const size_t &i, const size_t &j)const{
        return i*m_stride + j;
    }
    size_t m_start;
    size_t m_size;
    //extents
    size_t m_rows;
    size_t m_cols;
    //stride
    size_t m_stride;
};

template<typename T>
class matrix;

template<typename T>
class matrix_ref2D{
public:
    matrix_ref2D():m_ptr(nullptr){}
    matrix_ref2D(const matrix_slice &s, T *p):m_desc(s),m_ptr(p + s.m_start){}
    matrix_ref2D(const matrix_ref2D<T> &mref) = default;
    matrix_ref2D(matrix_ref2D<T> &&mref) = default;
    matrix_ref2D<T>& operator=(const matrix_ref2D<T> &mref) = default;
    matrix_ref2D<T>& operator =(matrix_ref2D<T> &&mref) = default;

    matrix_ref2D<T>& operator =(const matrix<T> &m);

    matrix_ref2D<T>& operator =(const T &val){
        for (size_t i = 0; i < rows(); ++i){
            for (size_t j = 0; j < cols(); ++j){
                m_ptr[m_desc(i,j)] = val;
            }
        }
        return *this;
    }

    matrix_ref1D<T> row(const size_t &i){
        assert(i < m_desc.m_rows); //Verificar que se esta accediendo al lugar correcto de memoria
        slice s(i*m_desc.m_stride,m_desc.m_cols,1);
        return matrix_ref1D<T>(s,data());
    }
    matrix_ref1D<const T> row(const size_t &i) const{
        assert(i < m_desc.m_rows); //Verificar que se esta accediendo al lugar correcto de memoria
        slice s(i*m_desc.m_stride,m_desc.m_cols,1);
        return matrix_ref1D<const T>(s,data());
    }
    matrix_ref1D<T> column(const size_t &i){
        assert(i < m_desc.m_cols);
        slice s(i,m_desc.m_rows,m_desc.m_cols);
        return matrix_ref1D<T>(s,data());
    }
    matrix_ref1D<const T> column(const size_t &i) const{
        assert(i < m_desc.m_cols);
        slice s(i,m_desc.m_rows,m_desc.m_cols);
        return matrix_ref1D<const T>(s,data());
    }

    T& operator ()(const size_t &i, const size_t &j){return m_ptr[m_desc(i,j)];}
    const T& operator ()(const size_t &i, const size_t &j) const{return m_ptr[m_desc(i,j)];}

    //Operaciones Aritmeticas
    template<typename Scalar>
    matrix_ref2D<T>& operator +=(const Scalar &val){
        if (val == 0) return *this;
        for (size_t i = 0; i < m_desc.m_rows; ++i)
            for (size_t j = 0; j < m_desc.m_cols; ++j)
                m_ptr[m_desc(i,j)] += val;
        return *this;
    }
    template<typename Scalar>
    matrix_ref2D<T>& operator -=(const Scalar &val){
        if (val == 0) return *this;
        for (size_t i = 0; i < m_desc.m_rows; ++i)
            for (size_t j = 0; j < m_desc.m_cols; ++j)
                m_ptr[m_desc(i,j)] -= val;
        return *this;
    }
    template<typename Scalar>
    matrix_ref2D<T>& operator *=(const Scalar &val){
        if (val == 1) return *this;
        for (size_t i = 0; i < m_desc.m_rows; ++i)
            for (size_t j = 0; j < m_desc.m_cols; ++j)
                m_ptr[m_desc(i,j)] *= val;
        return *this;
    }
    template<typename Scalar>
    matrix_ref2D<T>& operator /=(const Scalar &val){
        if (val == 1) return *this;
        for (size_t i = 0; i < m_desc.m_rows; ++i)
            for (size_t j = 0; j < m_desc.m_cols; ++j)
                m_ptr[m_desc(i,j)] /= val;
        return *this;
    }
    template<typename T2>
    matrix_ref2D<T>& operator +=(const matrix_ref2D<T2> &m1){
        assert(m_desc.m_rows == m1.rows() && m_desc.m_cols == m1.cols());
        for (size_t i = 0; i < m_desc.m_rows; ++i)
            for (size_t j = 0; j < m_desc.m_cols; ++j)
                m_ptr[m_desc(i,j)] += m1(i,j);
        return *this;
    }
    template<typename T2>
    matrix_ref2D<T>& operator -=(const matrix_ref2D<T2> &m1){
        assert(m_desc.m_rows == m1.rows() && m_desc.m_cols == m1.cols());
        for (size_t i = 0; i < m_desc.m_rows; ++i)
            for (size_t j = 0; j < m_desc.m_cols; ++j)
                m_ptr[m_desc(i,j)] -= m1(i,j);
        return *this;
    }

    size_t rows() const{return m_desc.m_rows;}
    size_t cols() const{return m_desc.m_cols;}

    size_t size() const{return m_desc.m_size;}
    T* data() {return m_ptr;}
    const T* data() const {return m_ptr;}

    void setToZero(){
        for (size_t i = 0; i < m_desc.m_rows; ++i)
            for (size_t j = 0; j < m_desc.m_cols; ++j)
                m_ptr[m_desc(i,j)] = 0;
    }
private:
    matrix_slice m_desc;
    T *m_ptr;
};
template<typename T>
ostream& operator <<(ostream& os,const matrix_ref2D<T> &m){
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
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline void sum(const matrix_ref2D<T1> &m1, const matrix_ref2D<T2> &m2,matrix_ref2D<RT> &res){
    assert(m1.rows() == m2.rows() && m1.cols() == m2.cols() && m1.rows() == res.rows() && m2.cols() == res.cols());
    for (size_t i = 0; i < res.rows(); ++i)
        for(size_t j = 0; j < res.cols(); ++j)
            res(i,j) = m1(i,j) + m2(i,j);
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline void sub(const matrix_ref2D<T1> &m1, const matrix_ref2D<T2> &m2,matrix_ref2D<RT> &res){
    assert(m1.rows() == m2.rows() && m1.cols() == m2.cols() && m1.rows() == res.rows() && m2.cols() == res.cols());
    for (size_t i = 0; i < res.rows(); ++i)
        for(size_t j = 0; j < res.cols(); ++j)
            res(i,j) = m1(i,j) - m2(i,j);
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline mathvector<RT> mult(const mathvector<T1> &v, const matrix_ref2D<T2> &mref){
    assert(v.size() == mref.rows());
    mathvector<RT> res(mref.cols(),0.0);
    for (size_t i = 0; i < mref.rows(); ++i) mult(v[i],mref.row(i),res);
    return res;
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline mathvector<RT> mult(const matrix_ref1D<T1> &v, const matrix_ref2D<T2> &mref){
    assert(v.size() == mref.rows());
    mathvector<RT> res(mref.cols(),0.0);
    for (size_t i = 0; i < mref.rows(); ++i) mult(v[i],mref.row(i),res);
    return res;
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline mathvector<RT> operator*(const mathvector<T1> &v, const matrix_ref2D<T2> &mref){
    return mult(v,mref);
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline mathvector<RT> operator*(const matrix_ref1D<T1> &v, const matrix_ref2D<T2> &mref){
    return mult(v,mref);
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline mathvector<RT> mult(const matrix_ref2D<T1> &mref,const mathvector<T2> &v){
    assert(mref.cols() == v.size());
    mathvector<RT> res(mref.rows(),0.0);
    for (size_t i = 0; i < mref.rows(); ++i){
        double temp = 0.0;
        for (size_t j = 0; j < mref.cols(); ++j)
            temp += mref(i,j)*v[j];
        res[i] = temp;
    }
    return res;
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline mathvector<RT> mult(const matrix_ref2D<T1> &mref,const matrix_ref1D<T2> &v){
    assert(mref.cols() == v.size());
    mathvector<RT> res(mref.rows(),0.0);
    for (size_t i = 0; i < mref.rows(); ++i){
        double temp = 0.0;
        for (size_t j = 0; j < mref.cols(); ++j)
            temp += mref(i,j)*v[j];
        res[i] = temp;
    }
    return res;
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline mathvector<RT> operator*(const matrix_ref2D<T1> &mref,const mathvector<T2> &v){
    return mult(mref,v);
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline mathvector<RT> operator*(const matrix_ref2D<T1> &mref,const matrix_ref1D<T2> &v){
    return mult(mref,v);
}

typedef matrix_ref2D<double> MatDoub2DRef;
#endif // MATRIX_REF2D_H
