#include "matrix_ref1d.h"
#include "mathvector.h"

template<typename T>
matrix_ref1D<T>& matrix_ref1D<T>::operator =(const mathvector<T> &vec){
    assert(this->size() == vec.size());
    for (size_t i = 0; i < m_desc.m_size; ++i) m_ptr[i*m_desc.m_stride] = vec[i];
    return *this;
}
