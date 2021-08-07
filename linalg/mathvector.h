#ifndef MATHVECTOR_H
#define MATHVECTOR_H

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <ostream>
#include <fstream>
#include <iomanip>
#include "my_concepts.h"
#include "matrix_ref1d.h"

using std::for_each;
using std::inner_product;
using std::vector;
using std::size_t;
using std::sqrt;
using std::ostream;
using std::ifstream;
using std::ofstream;
using mstdc::Common_type;
using mstdc::Convertible;

template<typename T>
class mathvector : public vector<T>{
public:
    explicit mathvector() : vector<T>(){}
    explicit mathvector(const size_t &count) : vector<T>(count){}
    explicit mathvector(const size_t &count, const T& value): vector<T>(count,value){}
    mathvector(const mathvector<T> &other) : vector<T>(other){}
    mathvector(mathvector<T> &&other) : vector<T>(other){}

    /*mathvector(const MatDoub1DRef &mref): vector<double>(mref.size()){
        for (size_t i = 0; i < size(); ++i) operator [](i) = mref(i);
    }*/

    mathvector<T>& operator =(const mathvector &other) = default;
    mathvector<T>& operator =(mathvector &&other) = default;
    mathvector<T>& operator =(const matrix_ref1D<T> &mref){
        resize(mref.size());
        for (size_t i = 0; i < this->size(); ++i) this->operator [](i) = mref(i);

        return *this;
    }

    matrix_ref1D<T> operator()(const size_t &start,const size_t &length,const size_t stride = 1){
        assert(start <= this->size() && length <= this->size()-start);

        slice s(start,length,stride);

        return matrix_ref1D<T>(s,this->data());
    }
    matrix_ref1D<const T> operator()(const size_t &start,const size_t &length,const size_t stride = 1)const{
        assert(start < this->size() && length <= this->size()-start);

        slice s(start,length,stride);

        return matrix_ref1D<const T>(s,this->data());
    }
    mathvector<T>& append(const mathvector<T> &v1){
        this->insert(this->end(),v1.begin(),v1.end());
    }
    mathvector<T>& append(const matrix_ref1D<T> &v1){
        for (size_t i = 0; i < v1.size(); ++i)
            this->push_back(v1[i]);
    }

    T& operator ()(const size_t &i){return this->operator [](i);}
    const T& operator ()(const size_t &i) const{return this->operator [](i);}

    double norm() const{
        double val = inner_product(this->begin(),this->end(),this->begin(),0.0);
        return sqrt(val);
    }
    double norm_2() const{
        return inner_product(this->begin(),this->end(),this->begin(),0.0);
    }
    double norm_l1() const{
        double res = 0.0;
        for_each(this->begin(),this->end(),[&res](const T &elem){res += std::abs(elem);});

        return res;
    }
    double norm_inf() const{
        double res = 0.0;
        for_each(this->begin(),this->end(),[&res](const T &elem){if (res < std::abs(elem)) res = elem;});

        return res;
    }
    template<typename T1>
    mathvector<T>& operator*=(const T1 &scalar){
        if (scalar == 1.0) return *this;
        if (scalar == -1.0){
            for_each(this->begin(),this->end(),[&scalar](double &val){val = -val;});
            return *this;
        }
        for_each(this->begin(),this->end(),[&scalar](double &val){val*=scalar;});
        return *this;
    }
    template<typename T1>
    mathvector<T>& operator/=(const T1 &scalar){
        for_each(this->begin(),this->end(),[&scalar](double &val){val/=scalar;});
        return *this;
    }
    template<typename T1>
    mathvector<T>& operator +=(const T1 &s){
        for_each(this->begin(),this->end(),[&s](double &val){val += s;});
        return *this;
    }
    template<typename T1>
    mathvector<T>& operator +=(const mathvector<T1> &other){
        assert(this->size() == other.size());
        size_t i = 0;
        for_each(this->begin(),this->end(),[&other,&i](double &val){val += other[i]; ++i;});
        return *this;
    }
    template<typename T1>
    mathvector<T>& operator +=(const matrix_ref1D<T1> &other){
        assert(this->size() == other.size());
        size_t i = 0;
        for_each(this->begin(),this->end(),[&other,&i](double &val){val += other[i]; ++i;});
        return *this;
    }
    template<typename T1>
    mathvector<T>& operator -=(const mathvector<T1> &other){
        if (other.size() == 0){
            return *this;
        }
        assert(this->size() == other.size());
        size_t i = 0;
        for_each(this->begin(),this->end(),[&other,&i](double &val){val -= other[i]; ++i;});
        return *this;
    }
    template<typename T1>
    mathvector<T>& operator -=(const matrix_ref1D<T1> &other){
        assert(this->size() == other.size());
        size_t i = 0;
        for_each(this->begin(),this->end(),[&other,&i](double &val){val -= other[i]; ++i;});
        return *this;
    }
    template<typename T1>
    mathvector<T>& operator -=(const T1 &s){
        for_each(this->begin(),this->end(),[&s](double &val){val -= s;});
        return *this;
    }

    double c_angle(const mathvector<T>& other){
        assert(this->size() == other.size());
        double res = inner_product(this->begin(),this->end(),other.begin(),0.0);
        return res /= norm()*other.norm();
    }
    mathvector<double> unitary(){
        mathvector<double> res(this->size());
        std::copy(this->begin(),this->end(),res.begin());
        return res/=norm();
    }
};
template<typename T>
inline mathvector<T> operator+(const mathvector<T> &v1,const mathvector<T> &v2){
    mathvector<T> res(v1);
    res+=v2;
    return res;
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline mathvector<RT> operator+(const mathvector<T1> &v1,const mathvector<T2> &v2){
    assert(v1.size() == v2.size());
    mathvector<RT> res(v1.size());
    auto it1 = v1.cbegin();
    auto it2 = v2.cbegin();
    for_each(res.begin(),res.end(),[&it1,&it2](RT &elem){elem = *it1 + *it2; ++it1; ++it2;});
    return res;
}
template<typename T>
inline mathvector<T> operator-(const mathvector<T> &v1,const mathvector<T> &v2){
    mathvector<T> res(v1);
    res-=v2;
    return res;
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline mathvector<RT> operator-(const mathvector<T1> &v1,const mathvector<T2> &v2){
    assert(v1.size() == v2.size());
    mathvector<RT> res(v1.size());
    auto it1 = v1.cbegin();
    auto it2 = v2.cbegin();
    for_each(res.begin(),res.end(),[&it1,&it2](RT &elem){elem = *it1 - *it2; ++it1; ++it2;});
    return res;
}
template<typename S>
inline mathvector<S> operator*(const S &val,const mathvector<S> &v){
    if (val == 0){
        return mathvector<S>(v.size(),0);
    }
    mathvector<S> res(v.size());
    auto it = v.cbegin();
    for_each(res.begin(),res.end(),[&val,&it](S &elem){elem = val*(*it); ++it;});
    return res;
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline mathvector<RT> operator*(const T1 &val,const mathvector<T2> &v){
    if (val == 0){
        return mathvector<RT>(v.size(),0);
    }
    mathvector<RT> res(v.size());
    auto it = v.cbegin();
    for_each(res.begin(),res.end(),[&val,&it](RT &elem){elem = val*(*it); ++it;});
    return res;
}
template<typename S>
inline mathvector<S> operator*(const mathvector<S> &v,const S &val){
    if (val == 0){
        return mathvector<S>(v.size(),0);
    }
    mathvector<S> res(v);
    return res*=val;
}
template<typename T1, typename T2, typename RT = Common_type<T1,T2>>
inline mathvector<RT> operator*(const mathvector<T2> &v,const T1 &val){
    if (val == 0){
        return mathvector<RT>(v.size(),0);
    }
    mathvector<RT> res(v.size());
    auto it = v.cbegin();
    for_each(res.begin(),res.end(),[&val,&it](RT &elem){elem = (*it)*val; ++it;});
    return res;
}
template<typename T>
inline T operator*(const mathvector<T> &v1,const mathvector<T> &v2){
    assert(v1.size() == v2.size());
    T res = inner_product(v1.begin(),v1.end(),v2.begin(),0.0);
    return res;
}
template<typename S>
inline mathvector<S> operator/(const mathvector<S> &v,const S &val){
    mathvector<S> res(v);
    return res/=val;
}
template<typename S,typename T,typename RT = Common_type<S,T>>
inline mathvector<S> operator/(const mathvector<T> &v,const S &val){
    mathvector<RT> res(v.size());
    auto it = v.cbegin();
    for_each(res.begin(),res.end(),[&val,&it](RT &elem){elem = (*it)/val; ++it;});
    return res;
}
inline double c_angle(const mathvector<double> &v1, const mathvector<double> &v2){
    return v1*v2/(v1.norm()*v2.norm());
}
template<typename T>
ostream& operator <<(ostream& os,const mathvector<T> &v){
    std::ios_base::fmtflags ff = std::ios::scientific;
    ff |= std::ios::showpos;
    os.setf(ff);
    for_each(v.begin(),v.end(),[&os](const double &elem){os << elem << std::endl;});
    os.unsetf(ff);
    return os;
}
template<typename T>
ifstream& operator >> (ifstream &in,mathvector<T> &vec){
    if (in.is_open()){
        for (size_t i = 0; i < vec.size(); ++i){
            in >> vec[i];
        }
    }
    return in;
}
template<typename T>
ofstream& operator << (ofstream &ofs,const mathvector<T> &vec){
    if (ofs.is_open()){
        std::ios_base::fmtflags ff = std::ios::scientific;
        ff |= std::ios::showpos;
        ofs.setf(ff);
        for (size_t i = 0; i < vec.size(); ++i){
            ofs << vec[i] << std::endl;
        }
        ofs.unsetf(ff);
    }

    return ofs;
}
template<typename T1,typename T2, typename RT = Common_type<T1,T2>>
inline void mult(const T1 &val,const matrix_ref1D<T2> &mref,mathvector<RT> &vec){
    assert(mref.size() == vec.size());

    if (val == 0.0) return;

    size_t i = 0;
    for_each(vec.begin(),vec.end(),[&val,&i,&mref](double &elem){elem += val*mref(i); ++i;});
}
template<typename T1,typename T2, typename RT = Common_type<T1,T2>>
inline RT operator*(const matrix_ref1D<T1> &mref,const mathvector<T2> &v){
    assert(mref.size() == v.size());
    RT res = 0.0;
    for (size_t i = 0; i < v.size(); ++i) res+=mref[i]*v[i];

    return res;
}
template<typename T1,typename T2, typename RT = Common_type<T1,T2>>
inline RT operator*(const mathvector<T1> &v,const matrix_ref1D<T2> &mref){
    assert(mref.size() == v.size());
    RT res = 0.0;
    for (size_t i = 0; i < v.size(); ++i) res+=mref[i]*v[i];

    return res;
}
template<typename T1,typename T2, typename RT = Common_type<T1,T2>>
inline mathvector<RT> operator+(const mathvector<T1> &v1,const matrix_ref1D<T2> mref){
    assert (v1.size() == mref.size());

    mathvector<RT> res(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) res[i] = v1[i]+mref[i];

    return res;
}
template<typename T1,typename T2, typename RT = Common_type<T1,T2>>
inline mathvector<RT> operator-(const mathvector<T1> &v1,const matrix_ref1D<T2> mref){
    assert (v1.size() == mref.size());

    mathvector<RT> res(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) res[i] = v1[i]-mref[i];

    return res;
}
template<typename T1,typename T2, typename RT = Common_type<T1,T2>>
inline mathvector<RT> operator+(const matrix_ref1D<T2> mref,const mathvector<T1> &v1){
    assert (v1.size() == mref.size());

    mathvector<RT> res(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) res[i] = mref[i]+v1[i];

    return res;
}
template<typename T1,typename T2, typename RT = Common_type<T1,T2>>
inline mathvector<RT> operator-(const matrix_ref1D<T2> mref,const mathvector<T1> &v1){
    assert (v1.size() == mref.size());

    mathvector<RT> res(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) res[i] = mref[i]-v1[i];

    return res;
}
template<typename T1,typename T2, typename RT = Common_type<T1,T2>>
mathvector<RT> operator*(const T1 &s,const matrix_ref1D<T2> &mref){
    mathvector<RT> res(mref.size());
    for (size_t i = 0; i < res.size(); ++i) res[i] = s*mref(i);

    return res;
}
template<typename T1,typename T2, typename RT = Common_type<T1,T2>>
mathvector<RT> operator*(const matrix_ref1D<T2> &mref,const T1 &s){
    mathvector<RT> res(mref.size());
    for (size_t i = 0; i < res.size(); ++i) res[i] = mref(i)*s;

    return res;
}
typedef mathvector<double> VecDoub;

#endif // MATHVECTOR_H
