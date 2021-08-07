#ifndef TRUSTREGION_H
#define TRUSTREGION_H

#include "linalg/directsolvers.h"
#include <string>

using std::string;

template<typename T>
class Trust_region{
private:
    T &funcd; //funcion objetivo
    double m_oreg; //region promedio
    double m_reg; //region inicial
    double m_sig; //parametro de aceptacion/rechazo de paso
    double fmin;
    size_t m_iter;

    //Solo se usa cuando B es definido positivo
    void Dogleg(const MatDoub &B, const VecDoub &g,VecDoub &p){
        Cholesky ch_solver(B);
        ch_solver.solve(g,p);
        //ch_solver.~Cholesky(); //Llamada explicita al destructor

        p*=-1.0;

        if (p.norm() <= m_reg) return;

        double temp1 = g*(B*g);
        double gnorm = g.norm();
        double gnorm2 = gnorm*gnorm;
        VecDoub pu = (-gnorm2/temp1)*g;
        //pu*=-1.0;

        double c = pu.norm_2() - m_reg*m_reg;
        p-=pu; //pb-bu

        double b = pu*p;
        double a = p.norm_2();

        double D = b*b - a*c; //Discriminante

        if (D < 0.0){ //usamos punto de cauchy
            p = (-m_reg/gnorm)*g;
            if (temp1 > 0.0){
                double t = gnorm*gnorm2/(m_reg*temp1);
                t = std::min(t,1.0);
                p*=t;
            }
        }else{

            double bet1 = sqrt(D);
            double bet2 = -bet1;

            bet1 -= b;
            bet2 -= b;

            bet1/=a;
            bet2/=a;

            bet1+=1.0;
            bet2+=1.0;

            //Resultado dentro de [0.2]
            if ((bet1 < 0.0 || bet1 > 2.0) && (bet2 < 0.0 || bet2 > 2.0)) throw(string("Trust_region : Dogleg : Error "
                "parametor t fuera del intervalo [0,2]"));

            double t;
            //t pertenece al intervalo [0,2]
            if (bet1 >= 0.0 && bet1 <= 2.0) t = bet1;
            else t = bet2;

            if (t <= 1.0){
                p = t*pu;
            }else{
                p*=(t-1);
                p+=pu;
            }
        }

    }
public:
    explicit Trust_region(T &func, const double &oreg, const double &reg, const double &sig):
        funcd(func),m_oreg(oreg),m_reg(reg),m_sig(sig){}

    size_t IterNumber() const{return m_iter;}
    double FMin() const{return fmin;}
    //Solo se usa cuando B es definido positivo
    void mininimize_dogleg(VecDoub &x){
        VecDoub p(x.size());
        fmin = funcd(x); //Valor de la funcion objetivo en x
        double ftemp;  // valor de la funcion en x+p
        VecDoub g;
        MatDoub B;
        funcd.d1f(x,g); //Gradiente de la funcion objetivo en x
        funcd.d2f(x,B); //Hessiano de la funcion objetivo en x
        bool XMOD = false; //Bandera true si x fue modificado ie x = x+p para recalcular g y B
        for (m_iter = 0 ; ; ++m_iter){

            if (XMOD){
                funcd.d1f(x,g); //Gradiente de la funcion objetivo en x
                funcd.d2f(x,B); //Hessiano de la funcion objetivo en x
            }

            if (g.norm() <= 1e-8) return;

            Dogleg(B,g,p);

            VecDoub xtemp = x+p;
            ftemp = funcd(xtemp);
            double rho = -(fmin - ftemp)/(g*p + 0.5*(p*B*p));

            if (rho < 0.25){
                m_reg*=0.25;
            }else{
                if (rho > 0.75 && p.norm() == m_reg){
                    //m_reg = std::min(2.0*m_reg,m_oreg);
                    m_reg*=2.0;
                }
            }

            if (rho > m_sig){
                x=xtemp;
                fmin = ftemp;
                XMOD = true;
            }else{
                XMOD = false;
            }
        }
    }

};

#endif // TRUSTREGION_H
