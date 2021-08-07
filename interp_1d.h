#ifndef INTERP_1D_H
#define INTERP_1D_H

#include "linalg/mathvector.h"

struct base_interp{
    size_t n, mm, jsav, cor, dj;
    const double *xx, *yy;
    base_interp(const VecDoub &x,const double *y,size_t m):
        n(x.size()),mm(m),jsav(0),cor(0),xx(&x[0]),yy(y){
        dj = std::min(1,(int)pow((double)n,0.25));
    }

    double interp(double x){
        int jlo = cor ? hunt(x) : locate(x);
        return rawinterp(jlo,x);
    }
    int locate(const double x);
    int hunt(const double x);

    double virtual rawinterp(int jlo,double x) = 0;
};


inline int base_interp::locate(const double x)
{
    int ju,jm,jl;

    if (n < 2 || mm < 2 || mm > n) throw("locate error size");

    bool ascnd = (xx[n-1] >= xx[0]); //monotony ascendent=true decrecent=false

    jl = 0;
    ju = n-1;
    while (ju-jl > 1){
        jm = (ju+jl) >> 1;
        if (x >= xx[jm] == ascnd) jl = jm;
        else ju = jm;
    }
    cor = std::abs(jl - jsav) > dj ? 0 : 1;
    jsav = jl;
    return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
}

inline int base_interp::hunt(const double x)
{
    int jl=jsav, jm, ju, inc=1;
    if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
    bool ascnd=(xx[n-1] >= xx[0]); //True if ascending order of table, false otherwise.

    if (jl < 0 || jl > n-1) { //Input guess not useful. Go immediately to bisection.
        jl=0;
        ju=n-1;
    } else {
        if (x >= xx[jl] == ascnd) { //Hunt up:

            for (;;) {
                ju = jl + inc;
                if (ju >= n-1) { ju = n-1; break;} //Off end of table.

                else if (x < xx[ju] == ascnd) break; //Found bracket.

                else { //Not done, so double the increment and try again.

                    jl = ju;
                    inc += inc;
                }
            }
        } else { //Hunt down:

            ju = jl;
            for (;;) {
                jl = jl - inc;
                if (jl <= 0) { jl = 0; break;}  //Off end of table.

                else if (x >= xx[jl] == ascnd) break; //Found bracket.

                else { //Not done, so double the increment and try again.

                    ju = jl;
                    inc += inc;

                }
            }
        }
    }

    while (ju-jl > 1) { //Hunt is done, so begin the final bisection phase:

        jm = (ju+jl) >> 1;
        if (x >= xx[jm] == ascnd)
            jl=jm;
        else
            ju=jm;
    }
    cor = abs(jl-jsav) > dj ? 0 : 1; //Decide whether to use hunt or locate next time.

    jsav = jl;

    return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));

}
struct linear_interp : base_interp
{

    linear_interp(const VecDoub &xv,const VecDoub &yv)
        : base_interp(xv,&yv[0],2) {}
    double rawinterp(int j, double x) {
        if (xx[j]==xx[j+1]) return yy[j]; //Table is defective, but we can recover.
        else return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
    }


};
#endif // INTERP_1D_H

