#ifndef CONVANALIZER_H
#define CONVANALIZER_H

#include "linalg/mathvector.h"
#include <ostream>
#include <iomanip>

using std::vector;
using std::size_t;
using std::ostream;

/*Clase para analizar la velocidad de convergencia de algunos metodos empleados
 * aqui se almacenas los valores de la funcion en cada iteracion el parametro encontrado
 * en la busqueda lineal y el punto de evaluacion. */
struct ConvAnal{
    size_t n;
    VecDoub alfas; //Parametro de busqueda lineal en cada iteracion
    vector<VecDoub> points; //puntos de evaluacion en cada iteracion
    VecDoub FuncVals; //Valores de la funcion en cada iteracion
    VecDoub conv_const; //Valor de la constante de convergencia para cada iteracion

    ConvAnal():n(0),alfas(0),points(0),FuncVals(0),conv_const(0){}

    void addIter(const double &alfa, const VecDoub &point, const double &funval){
        alfas.push_back(alfa);
        points.push_back(point);
        FuncVals.push_back(funval);
        ++n;
    }

    void CalcConvConsts(){
        conv_const.resize(n);
        conv_const[0] = 1.0;
        if (n > 1){
            VecDoub x_sol = points[n-1];
            for (size_t i = 1; i < n-1; ++i){
                conv_const[i] = (points[i] - x_sol).norm()/(points[i-1] - x_sol).norm();
            }
        }

    }
    void clear(){
        alfas.clear();
        points.clear();
        FuncVals.clear();
        conv_const.clear();
        n = 0;
    }
};
ostream& operator << (ostream &os, const ConvAnal &anal){
    os << std::setw(5) << "Iter" << std::setw(3) << "   "
       << std::setw(20) << "Eval Point" << std::setw(3) << "   "
       << std::setw(18) << "L Search" << std::setw(3) << "   "
       << std::setw(14) << "F Val" << std::setw(3) << "   "
       << std::setw(14) << "Conv Const" << std::endl;
    for (size_t i = 0; i < anal.n; ++i){
        os <<std::scientific << std::showpos << std::setw(5) << i << std::setw(3) << "   "
           << "(" << std::setw(7) << std::setprecision(4) << anal.points[i][0] << ","
           << std::setw(7) << std::setprecision(4) << anal.points[i][1] << ")" << std::setw(3) << "   " << std::setprecision(6)
           << std::setw(10) << anal.alfas[i] << std::setw(3) << "   "
           << std::setw(10) << anal.FuncVals[i] << std::setw(3) << "   "
           << std::setw(10) << anal.conv_const[i] << std::endl;

    }
    return os;
}
#endif // CONVANALIZER_H
