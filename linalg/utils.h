#ifndef UTILS_H
#define UTILS_H
/* retorna el valor absoluto de a con el signo de b */
#include <cmath>
namespace mutils {
inline double Sign(const double &a, const double &b){
    return (b < 0) ? -(std::abs(a)):std::abs(a);
}
inline int Sign(const int &a, const int &b){
    return (b < 0) ? -(std::abs(a)):std::abs(a);
}
}
static bool abs_compare(int a, int b)
{
    return (std::abs(a) < std::abs(b));
}
#endif // UTILS_H
