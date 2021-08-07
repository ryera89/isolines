#ifndef MATRIX_REF1D_H
#define MATRIX_REF1D_H

#include <cassert>
#include <ostream>
#include <iomanip>
#include "my_concepts.h"

using std::ostream;
using mstdc::Common_type;
using mstdc::Convertible;

struct slice
{
    slice() = default;
    slice(const size_t &start,const size_t &size, const size_t &stride):
        m_start(start),m_size(size),m_stride(stride){}

    size_t operator()(const size_t &i,const size_t &j)const{
        return i*m_stride + j;
    }

    size_t m_start;
    size_t m_size;
    size_t m_stride;
};

template <typename T>
class mathvector;

template<typename T>
class matrix_ref1D{
public:
    matrix_ref1D():m_ptr(nullptr){}
    matrix_ref1D(const  slice &s, T* p):m_desc(s),m_ptr(p + s.m_start){}
    matrix_ref1D(const matrix_ref1D<T> &mref) = default;
    matrix_ref1D(matrix_ref1D<T> &&mref) = default;
    matrix_ref1D<T>& operator = (const matrix_ref1D<T> &mref) = default;
    matrix_ref1D<T>& operator = (matrix_ref1D<T> &&mref) = default;

    matrix_ref1D<T>& operator =(const mathvector<T> &vec);

    T& operator ()(const size_t &i){return m_ptr[i*m_desc.m_stride];}
    const T& operator()(const size_t &i) const{return m_ptr[i*m_desc.m_stride];}

    T& operator [](const size_t &i){return m_ptr[i*m_desc.m_stride];}
    const T& operator[](const size_t &i) const{return m_ptr[i*m_desc.m_stride];}

    //Operaciones Aritmeticas
    template<typename Scalar>
    matrix_ref1D<T>& operator +=(const Scalar &val){
        if (val == 0) return *this;
        for (size_t i = 0; i < m_desc.m_size; ++i) m_ptr[i*m_desc.m_stride] += val;
        return *this;
    }
    template<typename Scalar>
    matrix_ref1D<T>& operator -=(const Scalar &val){
        if (val == 0) return *this;
        for (size_t i = 0; i < m_desc.m_size; ++i) m_ptr[i*m_desc.m_stride] -= val;
        return *this;
    }
    template<typename Scalar>
    matrix_ref1D<T>& operator *=(const Scalar &val){
        if (val == 1) return *this;
        for (size_t i = 0; i < m_desc.m_size; ++i) m_ptr[i*m_desc.m_stride] *= val;
        return *this;
    }
    template<typename Scalar>
    matrix_ref1D<T>& operator /=(const Scalar &val){
        if (val == 1) return *this;
        for (size_t i = 0; i < m_desc.m_size; ++i) m_ptr[i*m_desc.m_stride] /= val;
        return *this;
    }
    template<typename T2>
    matrix_ref1D<T>& operator +=(const matrix_ref1D<T2> &m1){
        assert(size() == m1.size());
        for (size_t i = 0; i < m_desc.m_size; ++i) m_ptr[i*m_desc.m_stride] += m1(i);
        return *this;
    }
    template<typename T2>
    matrix_ref1D<T>& operator -=(const matrix_ref1D<T2> &m1){
        assert(size() == m1.size());
        for (size_t i = 0; i < m_desc.m_size; ++i) m_ptr[i*m_desc.m_stride] -= m1(i);
        return *this;
    }
    void setToZero(){
         for (size_t i = 0; i < m_desc.m_size; ++i) m_ptr[i*m_desc.m_stride] = 0;
    }

    size_t size() const{return m_desc.m_size;}
private:
    slice m_desc;
    T *m_ptr;

};
template<typename T>
ostream& operator <<(ostream& os,const matrix_ref1D<T> &v){
    std::ios_base::fmtflags ff = std::ios::scientific;
    ff |= std::ios::showpos;
    os.setf(ff);
    for (size_t i = 0; i < v.size(); ++i)
        os << v(i) << std::endl;

    os.unsetf(ff);

    return os;
}

template<typename T1,typename T2, typename RT = Common_type<T1,T2>>
inline void mult(const T1 &val, const matrix_ref1D<T2> &mref, matrix_ref1D<RT> &res){
    if (val == 0) return;
    for (size_t i = 0; i < mref.size(); ++i) res(i) += val*mref(i);
}
template<typename T1,typename T2, typename RT = Common_type<T1,T2>>
inline RT mult(const matrix_ref1D<T1> &mref1,const matrix_ref1D<T2> &mref2){
    assert(mref1.size() == mref2.size());
    RT res = 0;
    for (size_t i = 0; i < mref1.size(); ++i) res += mref1(i)*mref2(i);
    return res;
}
template<typename T1,typename T2, typename RT = Common_type<T1,T2>>
inline RT operator *(const matrix_ref1D<T1> &mref1,const matrix_ref1D<T2> &mref2){
    return mult(mref1,mref2);
}
typedef matrix_ref1D<double> MatDoub1DRef;


#endif // MATRIX_REF1D_H
