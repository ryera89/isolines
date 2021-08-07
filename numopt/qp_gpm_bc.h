#ifndef QP_GPM_BC_H
#define QP_GPM_BC_H

#include "linalg/mathvector.h"
#include <numeric>
#include <set>
#include "linalg/mathmatrix.h"

using std::set;

using std::numeric_limits;
/* Metodo de Gradiente proyectado para resolver problemas
 * de optimizacion tipo QP q(x) = 1/2*(x*G*x) + c*x con restricciones de caja
 * Numerical Optimization 2ed 2006, J. Nocedal y S.J. Wright,
 * pp 16.7 */
//template<class T>
class qp_gpm_bc{
private:
    size_t its;

    void detCauchyPoint(const MatDoub &G,const VecDoub &c,const VecDoub &g,
                        VecDoub &x,const VecDoub &lb, const VecDoub &ub){
        size_t n = x.size();
        //vector de breakpoints
        VecDoub tb(n);

        //Breakpoint
        for (size_t i = 0; i < n; ++i){
            if (g[i] < 0.0 && ub[i] < numeric_limits<double>::infinity()){
                tb[i] = (x[i] - ub[i])/g[i];
            }else{
                if (g[i] > 0.0 && lb[i] > numeric_limits<double>::lowest()){
                    tb[i] = (x[i] - lb[i])/g[i];
                }else{
                    tb[i] = numeric_limits<double>::max();
                }
            }
        }
        //Elimina los duplicados y ordena el set 0<t1<t2<t3<...<tl  l <= n
        set<double> intervals(tb.begin(),tb.end());
        intervals.erase(0.0);

        VecDoub p(n);

        double prev_t = 0.0;
        //for (size_t i = 0; i < n; ++i) p[i] = -g[i];

        //TODO: Optimizar la estrategia de update de p
        for (auto it = intervals.begin(); it != intervals.end();++it){
            for (size_t i = 0; i < n; ++i){
                if (prev_t < tb[i]) p[i] = -g[i];
                else p[i] = 0;
            }
            VecDoub vaux(G*p);
            double fp = c*p + x*vaux;
            if (fp > 0.0) return;
            double fpp = p*vaux;
            double delta_t = -fp/fpp;
            if (delta_t < *(it) - prev_t){
                x+=delta_t*p;
                return;
            }
            prev_t = *(it);
        }
    }
    void cg_redsys_qp_bc(const MatDoub &G, VecDoub &c, VecDoub &x, const vector<size_t> &vperm, const size_t &n,
                         const size_t &m,const VecDoub &lb, const VecDoub &ub, const double tol = 1.0e-8){
        //VecDoub vaux(n-m);
        for (size_t i = m; i < n; ++i){
            double temp = 0.0;
            for (size_t j = m; j < n; ++j){
                temp += G(vperm[i],vperm[j])*x[vperm[j]];
            }
            c[i-m] += temp;
        }
        double sq_prev_norm = c.norm_2();

        VecDoub p(c.size());

         for(size_t i = 0; i < p.size(); ++i) p[i] = -c[i];

         VecDoub vaux(c.size());
         //TODO cambiar condicion de parada
         bool bound_reached = false;
         while(!bound_reached && sqrt(sq_prev_norm) > tol){
             //Z_t*P*G*P_t*Z*p P: permutaciones
             for (size_t i = m; i < n; ++i){
                 double temp = 0.0;
                 for (size_t j = m; j < n; ++j){
                     temp += G(vperm[i],vperm[j])*p[j-m];
                 }
                 vaux[i-m] = temp;
             }
             double a = sq_prev_norm/(p*vaux);

             for (size_t i = m; i < n; ++i){
                 double temp = x[vperm[i]] + a*p[i-m];
                 if (temp >= ub[vperm[i]]){
                     x[vperm[i]] = ub[vperm[i]];
                     bound_reached = true;
                     continue;
                     //break;
                 }
                 if (temp <= lb[vperm[i]]){
                    x[vperm[i]] = lb[vperm[i]];
                    bound_reached = true;
                    continue;
                    //break;
                 }
                 x[vperm[i]] = temp;
             }

             c += a*vaux;
             double sq_new_norm = c.norm_2();
             double bet = sq_new_norm/sq_prev_norm;
             sq_prev_norm = sq_new_norm;
             p*=bet;
             p-=c;
         }
    }
    void create_pvec(const VecDoub &x, const VecDoub &lb, const VecDoub &ub,vector<size_t> &vperm,
                     const size_t &n,size_t &m){
        //size_t n = x.size();
        vperm.clear();
        vperm.reserve(n);
        vector<size_t> v_aux;
        v_aux.reserve(n);
        for (size_t i = 0; i < n; ++i){
            if (x[i] == lb[i] || x[i] == ub[i]){
                //xy.push_back(x[i]);
                vperm.push_back(i);
            }else{
                v_aux.push_back(i);
            }
        }
        m = vperm.size();
        vperm.insert(vperm.end(),v_aux.begin(),v_aux.end());
    }
public:
    qp_gpm_bc(){}
    size_t iter()const{return its;}

    //g: gradiente en x   g = Gx + c
    //lb: vector de limites inferiores si xi no tiene lbi se pasa como -lowest
    //ub: vector de limites superiores si xi no tiene ubi se pasa como infinity
    void minimize(const MatDoub &G, const VecDoub &c,const VecDoub &lb,
                  const VecDoub &ub,VecDoub &x){

        bool KKT_COND = false;
        size_t n = x.size();
        vector<size_t> vperm(n);
        for (size_t i = 0; i < n; ++i) vperm[i] = i;
        //vperm.reserve(n);
        size_t m = 0;
        VecDoub lamda;
        lamda.reserve(n);
        for (its = 0; ;++its){
            VecDoub g = G*x + c;
            //if (g.norm() < 1.0e-8) break;
            //create_pvec(x,lb,ub,vperm,n,m);
            lamda.resize(m,0.0);
            for (size_t i = 0; i < m; ++i){
                double temp = g[vperm[i]];
                if (temp >= 0.0){
                    lamda[i] = temp;
                    KKT_COND = true;
                }
                else{
                    KKT_COND = false;
                    break;
                }
            }
            for (size_t i = m; i < n; ++i){
                if (abs(g[vperm[i]]) > 1.0e-8){
                    KKT_COND = false;
                    break;
                }else{
                    KKT_COND = true;
                }
            }
            if (KKT_COND) break;
            //TODO: Verificar KKT
            //VecDoub xy;
            //vector<size_t> v_aux;
            //v_aux.reserve(n);
            detCauchyPoint(G,c,g,x,lb,ub);
            //Determinacion de vector de permutaciones
            create_pvec(x,lb,ub,vperm,n,m);
            //size_t m = xy.size();
            //size_t m = vperm.size();
            //Inserto el resto de las permutaciones
            //vperm.insert(vperm.end(),v_aux.begin(),v_aux.end());
            //TODO Liberar v_aux;
            VecDoub cz(n-m);

            //Determinando cz
            for (size_t i = m; i < n; ++i){
                double temp = 0.0;
                for (size_t j = 0; j < m; ++j){
                    temp += G(vperm[i],vperm[j])*x[vperm[j]];
                }
                temp += c[vperm[i]];
                cz[i-m] = temp;
            }

            cg_redsys_qp_bc(G,cz,x,vperm,n,m,lb,ub);
        }
    }

};

#endif // QP_GPM_BC_H
