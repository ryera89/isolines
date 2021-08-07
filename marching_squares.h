#ifndef MARCHING_SQUARES_H
#define MARCHING_SQUARES_H

#include "linalg/mathmatrix.h"
#include <QPointF>
#include <QList>

class marching_squares{
private:
    enum LUT{EMPTY, BLEFT, BRIGHT, BRIGHT_BLEFT, URIGHT, URIGHT_BLEFT,
            URIGHT_BRIGHT, URIGHT_BRIGHT_BLEFT, ULEFT, ULEFT_BLEFT, ULEFT_BRIGHT,
            ULEFT_BRIGHT_BLEFT,ULEFTH_URIGHT,ULEFT_URIGHT_BLEFT,ULEFT_URIGHT_BRIGHT,FULL};

    matrix<uint8_t> binary_image;
    matrix_ref2D<uint8_t> bin_ref;

    int getCellIndex(){
       int index = 0;
       if (bin_ref(0,0)) index|=8;
       if (bin_ref(0,1)) index|=4;
       if (bin_ref(1,1)) index|=2;
       if (bin_ref(1,0)) index|=1;

       return index;
    }
    double lininterp(const double &x,const double x1,const double x0, const double &y1,const double &y0){
        return y0 + (x-x0)*(y1-y0)/(x1-x0);
    }

    QPair<QPointF,QPointF> getIsolineCellCoordinates(const double &isoline,const int &index, const size_t &i,
                                             const size_t j,const MatDoub &Z,const VecDoub &x, const VecDoub &y){


        //std::cout <<  index << std::endl;
        switch (index) {
        case BLEFT:{
            double y0 = y[i];
            double y1 = y[i+1];
            double x0 = Z(i,j);
            double x1 = Z(i+1,j);
            double yy = lininterp(isoline,x1,x0,y1,y0);
            QPointF p1(x[j],yy);
            y0 = x[j];
            y1 = x[j+1];
            x0 = Z(i+1,j);
            x1 = Z(i+1,j+1);
            double xx = lininterp(isoline,x1,x0,y1,y0);
            QPointF p2(xx,y[i+1]);
            //std::cout << p1.x() << " " << p1.y() << " " << p2.x() << " " << p2.y() << std::endl;
            return QPair<QPointF,QPointF>(p1,p2);
        }
        case BRIGHT:{
            double y0 = x[j];
            double y1 = x[j+1];
            double x0 = Z(i+1,j);
            double x1 = Z(i+1,j+1);
            double xx = lininterp(isoline,x1,x0,y1,y0);
            QPointF p1(xx,y[i+1]);
            y0 = y[i];
            y1 = y[i+1];
            x0 = Z(i,j+1);
            x1 = Z(i+1,j+1);
            double yy = lininterp(isoline,x1,x0,y1,y0);
            QPointF p2(x[j+1],yy);
            //std::cout << p1.x() << " " << p1.y() << " " << p2.x() << " " << p2.y() << std::endl;
            return QPair<QPointF,QPointF>(p1,p2);
        }
        case BRIGHT_BLEFT:{
            double y0 = y[i];
            double y1 = y[i+1];
            double x0 = Z(i,j);
            double x1 = Z(i+1,j);
            double yy = lininterp(isoline,x1,x0,y1,y0);
            QPointF p1(x[j],yy);
            y0 = y[i];
            y1 = y[i+1];
            x0 = Z(i,j+1);
            x1 = Z(i+1,j+1);
            yy = lininterp(isoline,x1,x0,y1,y0);
            QPointF p2(x[j+1],yy);
            //std::cout << p1.x() << " " << p1.y() << " " << p2.x() << " " << p2.y() << std::endl;
            return QPair<QPointF,QPointF>(p1,p2);
        }
        case  URIGHT:{
            double y0 = x[j];
            double y1 = x[j+1];
            double x0 = Z(i,j);
            double x1 = Z(i,j+1);
            double xx = lininterp(isoline,x1,x0,y1,y0);
            QPointF p1(xx,y[i]);
            y0 = y[i];
            y1 = y[i+1];
            x0 = Z(i,j+1);
            x1 = Z(i+1,j+1);
            double yy = lininterp(isoline,x1,x0,y1,y0);
            QPointF p2(x[j+1],yy);
            //std::cout << p1.x() << " " << p1.y() << " " << p2.x() << " " << p2.y() << std::endl;
            return QPair<QPointF,QPointF>(p1,p2);
        }
        case URIGHT_BRIGHT:{
            double y0 = x[j];
            double y1 = x[j+1];
            double x0 = Z(i,j);
            double x1 = Z(i,j+1);
            double xx = lininterp(isoline,x1,x0,y1,y0);
            QPointF p1(xx,y[i]);
            y0 = x[j];
            y1 = x[j+1];
            x0 = Z(i+1,j);
            x1 = Z(i+1,j+1);
            xx = lininterp(isoline,x1,x0,y1,y0);
            QPointF p2(xx,y[i+1]);
            //std::cout << p1.x() << " " << p1.y() << " " << p2.x() << " " << p2.y() << std::endl;
            return QPair<QPointF,QPointF>(p1,p2);
        }
        case URIGHT_BRIGHT_BLEFT:{
            double y0 = y[i];
            double y1 = y[i+1];
            double x0 = Z(i,j);
            double x1 = Z(i+1,j);
            double yy = lininterp(isoline,x1,x0,y1,y0);
            QPointF p1(x[j],yy);

            y0 = x[j];
            y1 = x[j+1];
            x0 = Z(i,j);
            x1 = Z(i,j+1);
            double xx = lininterp(isoline,x1,x0,y1,y0);
            QPointF p2(xx,y[i]);
            //std::cout << p1.x() << " " << p1.y() << " " << p2.x() << " " << p2.y() << std::endl;
            return QPair<QPointF,QPointF>(p1,p2);
        }
        case ULEFT:{
            double y0 = y[i];
            double y1 = y[i+1];
            double x0 = Z(i,j);
            double x1 = Z(i+1,j);
            double yy = lininterp(isoline,x1,x0,y1,y0);
            QPointF p1(x[j],yy);

            y0 = x[j];
            y1 = x[j+1];
            x0 = Z(i,j);
            x1 = Z(i,j+1);
            double xx = lininterp(isoline,x1,x0,y1,y0);
            QPointF p2(xx,y[i]);
            //std::cout << p1.x() << " " << p1.y() << " " << p2.x() << " " << p2.y() << std::endl;
            return QPair<QPointF,QPointF>(p1,p2);
        }
        case ULEFT_BLEFT:{
            double y0 = x[j];
            double y1 = x[j+1];
            double x0 = Z(i,j);
            double x1 = Z(i,j+1);
            double xx = lininterp(isoline,x1,x0,y1,y0);
            QPointF p1(xx,y[i]);
            y0 = x[j];
            y1 = x[j+1];
            x0 = Z(i+1,j);
            x1 = Z(i+1,j+1);
            xx = lininterp(isoline,x1,x0,y1,y0);
            QPointF p2(xx,y[i+1]);
            //std::cout << p1.x() << " " << p1.y() << " " << p2.x() << " " << p2.y() << std::endl;
            return QPair<QPointF,QPointF>(p1,p2);
        }
        case ULEFT_BRIGHT_BLEFT:{
            double y0 = x[j];
            double y1 = x[j+1];
            double x0 = Z(i,j);
            double x1 = Z(i,j+1);
            double xx = lininterp(isoline,x1,x0,y1,y0);
            QPointF p1(xx,y[i]);
            y0 = y[i];
            y1 = y[i+1];
            x0 = Z(i,j+1);
            x1 = Z(i+1,j+1);
            double yy = lininterp(isoline,x1,x0,y1,y0);
            QPointF p2(x[j+1],yy);
            //std::cout << p1.x() << " " << p1.y() << " " << p2.x() << " " << p2.y() << std::endl;
            return QPair<QPointF,QPointF>(p1,p2);
        }
        case ULEFTH_URIGHT:{
            double y0 = y[i];
            double y1 = y[i+1];
            double x0 = Z(i,j);
            double x1 = Z(i+1,j);
            double yy = lininterp(isoline,x1,x0,y1,y0);
            QPointF p1(x[j],yy);
            y0 = y[i];
            y1 = y[i+1];
            x0 = Z(i,j+1);
            x1 = Z(i+1,j+1);
            yy = lininterp(isoline,x1,x0,y1,y0);
            QPointF p2(x[j+1],yy);
            //std::cout << p1.x() << " " << p1.y() << " " << p2.x() << " " << p2.y() << std::endl;
            return QPair<QPointF,QPointF>(p1,p2);
        }
        case ULEFT_URIGHT_BLEFT:{
            double y0 = x[j];
            double y1 = x[j+1];
            double x0 = Z(i+1,j);
            double x1 = Z(i+1,j+1);
            double xx = lininterp(isoline,x1,x0,y1,y0);
            QPointF p1(xx,y[i+1]);
            y0 = y[i];
            y1 = y[i+1];
            x0 = Z(i,j+1);
            x1 = Z(i+1,j+1);
            double yy = lininterp(isoline,x1,x0,y1,y0);
            QPointF p2(x[j+1],yy);
            //std::cout << p1.x() << " " << p1.y() << " " << p2.x() << " " << p2.y() << std::endl;
            return QPair<QPointF,QPointF>(p1,p2);
        }
        case ULEFT_URIGHT_BRIGHT:{
            double y0 = y[i];
            double y1 = y[i+1];
            double x0 = Z(i,j);
            double x1 = Z(i+1,j);
            double yy = lininterp(isoline,x1,x0,y1,y0);
            QPointF p1(x[j],yy);
            y0 = x[j];
            y1 = x[j+1];
            x0 = Z(i+1,j);
            x1 = Z(i+1,j+1);
            double xx = lininterp(isoline,x1,x0,y1,y0);
            QPointF p2(xx,y[i+1]);
            //std::cout << p1.x() << " " << p1.y() << " " << p2.x() << " " << p2.y() << std::endl;
            return QPair<QPointF,QPointF>(p1,p2);
        }
        default:
            break;
        }
    }
    void getBinaryRepresentation(const MatDoub &Z,const double &isovalue){
        //obtengo la representacion binaria
        for (size_t i = 0; i < Z.rows(); ++i){
            for (size_t j = 0; j < Z.cols(); ++j){
                if (Z(i,j) > isovalue) binary_image(i,j) = 1;
                else binary_image(i,j) = 0;
            }
        }
    }

public:
    marching_squares() = default;

    QVector<QPair<QPointF,QPointF>> getContour(const VecDoub &x, const VecDoub &y,const MatDoub &Z, const double &isovalue);

};

inline QVector<QPair<QPointF,QPointF>> marching_squares::getContour(const VecDoub &x, const VecDoub &y,const MatDoub &Z, const double &isovalue)
{
    assert(x.size() == Z.cols() && y.size() == Z.rows());

    size_t m = Z.rows();
    size_t n = Z.cols();
    binary_image.resize(m,n);
    getBinaryRepresentation(Z,isovalue);

    QVector<QPair<QPointF,QPointF>> points;
    //Recorro las celdas
    for (size_t i = 0; i < m-1; ++i){
        for (size_t j = 0; j < n-1; ++j){
            bin_ref = binary_image(slice(i,2,1),slice(j,2,1));

            int index = getCellIndex();
            //if (index == 5 || index == 10) std::cout << index << std::endl;
            switch (index) {
            case EMPTY:
                break;
            case FULL:
                break;
            case URIGHT_BLEFT:
                break;
            case ULEFT_BRIGHT:
                break;
            default:
                points.append(getIsolineCellCoordinates(isovalue,index,i,j,Z,x,y));
                break;
            }
        }
    }
     return points;
}
#endif // MARCHING_SQUARES_H
