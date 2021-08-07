#ifndef MY_CONCEPTS_H
#define MY_CONCEPTS_H

#include <type_traits>

namespace mstdc {
template<bool B,typename T = void>
using Enable_if = typename std::enable_if<B,T>::type;

template<typename T1,typename T2>
using Common_type = typename std::common_type<T1,T2>::type;

template<typename F,typename T>
constexpr bool Convertible()
{
       return std::is_convertible<F,T>::value;
}

constexpr bool All(){return true;}

template<typename... Args>
constexpr bool All(bool b, Args... args)
{
    return b && All(args...);
}
template<typename T, typename U>
constexpr bool Same()
{
    return std::is_same<T,U>::value;
}

constexpr bool Some() {return false;}

template<typename... Args>
constexpr bool Some(bool b, Args... args)
{
    return b || Some(args...);
}

}
#endif // MY_CONCEPTS_H
