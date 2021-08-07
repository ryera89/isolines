#include "matrix_ref2d.h"
#include "mathmatrix.h"

template<typename T>
inline matrix_ref2D<T>& matrix_ref2D<T>::operator =(const matrix<T> &m)
{
    assert(m.rows() == rows() && m.cols() == cols());

    for (size_t i = 0; i < rows(); ++i){
        for (size_t j = 0; j < cols(); ++j){
            m_ptr[m_desc(i,j)] = m(i,j);
        }
    }
    return *this;
}
