#include <QApplication>
#include <iostream>
#include "examples/functors.h"
#include "numopt/cgfrpr.h"
#include "isolineswindow.h"

/*void determineEjem1604Isolines(contourWidget &cw){
    VecDoub x(101);
    VecDoub y(101);
    for (size_t i = 0; i <= 100; ++i){
        x[i] = -2.5 + 0.1*i;
        y[i] = -4.0 + 0.1*i;
    }
    MatDoub Z(101,101);

    ejemplo_16_04 func;
    VecDoub xi(2);

    for (size_t i = 0; i <= 100; ++i){
        for (size_t j = 0; j <= 100; ++j){
            xi[0] = x[j];
            xi[1] = y[i];
            Z(i,j) = func(xi);
        }
    }
    cw.contour(x,y,Z,5);
    cw.contour(x,y,Z,10);
    cw.contour(x,y,Z,20);
    cw.contour(x,y,Z,30);
    cw.contour(x,y,Z,40);
    cw.contour(x,y,Z,50);
    cw.contour(x,y,Z,60);
    cw.contour(x,y,Z,70);
    cw.contour(x,y,Z,80);
}*/
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    QLocale::setDefault(QLocale::C);
    isolinesWindow w;
    w.show();
    //contourWidget cw;
    //determineEjem1604Isolines(cw);
    //cw.show();
    return a.exec();
}
