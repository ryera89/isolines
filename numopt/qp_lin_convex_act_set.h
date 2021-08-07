#ifndef QP_LIN_CONVEX_ACT_SET_H
#define QP_LIN_CONVEX_ACT_SET_H

#include "linalg/qr.h"
#include <set>
#include "linalg/itermethods.h"

using std::set;
using std::pair;

class qp_lin_convex_act_set{
  private:
    size_t index(const size_t &i, const size_t &j,const size_t &n)
    {
        size_t ind = (i <= j) ? (i*n + j - i*(i+1)/2) : (j*n + i - j*(j+1)/2);
        return ind;
    }
    void computeFeasibleStart(){}
    void computeW0(set<size_t> &W,set<size_t> &Win,const MatDoub &Ainq,const VecDoub &binq,
                   const VecDoub &x){
        VecDoub temp(Ainq*x);
        for (size_t i = 0; i < temp.size(); ++i){
            if (temp[i] == binq[i]) W.insert(i);
            else Win.insert(i);
        }
    }
    void determine_Y_Z_R_QRDCMP(const MatDoub &Aeq, MatDoub &Y,MatDoub &Z,VecDoub &R){
        size_t m = Aeq.rows();
        size_t n = Aeq.cols();
        //QR DCMP
        qrdcmp qr(Aeq.transpose());
        qr.householder_qr();

        //EXTRAYENDO Y y Z
        Y = qr.Q(slice(0,n,1),slice(0,m,1));
        Z = qr.Q(slice(0,n,1),slice(m,n-m,1));

        //Almacenando qr.R en un el vector R
        for (size_t i = 0; i < m; ++i)
            for (size_t j = i; j < m; ++j)
                R[index(i,j,m)] = qr.R(i,j);
    }
    void determinePy(const VecDoub &R,const size_t &meq,const size_t &minq_act, VecDoub &p_y){
        size_t m = meq + minq_act;
        p_y.resize(m);
        for (size_t k = 0; k < m; ++k){
            double sum = 0.0;
            for (size_t i = 0; i < k; ++i){
                sum += p_y[i]*R[index(i,k,m)];
            }
                p_y[k] = -sum/R[index(k,k,m)];
        }
    }
    void determinePyOnConstraintDeactivation(const size_t &c_pos, const VecDoub &R, VecDoub &p_y){
        p_y.pop_back();
        size_t m = p_y.size();
        for (size_t k = c_pos; k < m; ++k){
            double sum = 0.0;
            for (size_t i = 0; i < k; ++i){
                sum += p_y[i]*R[index(i,k,m)];
            }
                p_y[k] = -sum/R[index(k,k,m)];
        }

    }
    void determinePyOnConstraintActivation(const VecDoub &R, VecDoub &p_y){
        size_t k = p_y.size();
        size_t m = k+1;
        double sum = 0.0;
        for (size_t i = 0; i < k; ++i){
            sum += p_y[i]*R[index(i,k,m)];
        }
        p_y.push_back(-sum/R[index(k,k,m)]);

    }
    //k: restriccion activada
    void updateFactorizationOnConstrainActivation(const MatDoub &Ainq,const size_t &k
                                                  ,MatDoub &Q,MatDoub2DRef &Y,MatDoub2DRef &Z,VecDoub &R){

        //mathvector a(Ainq.row(k).size());
        // for (size_t i = 0; i < )
        //MatDoub1DRef refaux = Ainq.row(k);
        VecDoub w(Ainq.row(k)*Z);

        size_t nn = Q.rows();
        size_t mm = Y.cols();
        size_t n = w.size(); //nn-mm

        //vector<GivensRotation> giv_rot_vec(n-1);
        //vector<pair<size_t,size_t>> index_vec(n-1);
        for (size_t i = n-1; i > 0; --i){
            GivensRotation giv_rot(w[i-1],w[i]);
            giv_rot.apply_givens_to_vector(i-1,i,w);
            //giv_rot_vec[i-1] = giv_rot; //Storing the givens rotations
            //index_vec[i-1] = pair<size_t,size_t>(i-1,i); //Storing indexes
            giv_rot.apply_givens_to_cols(i-1,i,Z);

        }
        //MatDoub Ztemp(Z);

        /*for (long int i = n-2; i >= 0; --i){
            pair<size_t,size_t> ptemp = index_vec[i];
            giv_rot_vec[i].apply_givens_to_cols(ptemp.first,ptemp.second,Z);
        }*/

        size_t m = mm+1;
        VecDoub Rtemp(m*(m+1)/2);
        if (mm){
            VecDoub temp(Ainq.row(k)*Y);
            //size_t mtemp = m-1;
            for (size_t i = 0; i < m; ++i){
                for (size_t j = i; j < m; ++j){
                    if (j < mm){
                        Rtemp[index(i,j,m)] = R[index(i,j,mm)];
                    }else{
                        if (i < mm) Rtemp[index(i,j,m)] = temp[i];
                        else Rtemp[index(i,j,m)] = w[0];
                    }
                }
            }
        }else{
          Rtemp[0] = w[0];
        }
        R = Rtemp;

        //Y.append_column(Z(0));
        //Z = Z(slice(0,nn,1),slice(1,n-1,1));
        //Y.append_column(Z(0));

        Y = Q(slice(0,nn,1),slice(0,m,1));
        Z = Q(slice(0,nn,1),slice(m,n-1,1));

    }
    //k restriccion desactivada
    void updateFactorizationOnConstrainDeactivation(const size_t &k,MatDoub &Q,
                                                    MatDoub2DRef &Y,MatDoub2DRef &Z,VecDoub &R){
        size_t nn = Q.rows();
        size_t m = Y.cols();

        //vector<GivensRotation> giv_rot_vec(n-1);
        //vector<pair<size_t,size_t>> index_vec(n-1);
        size_t mtemp = m-1;
        for(size_t i = k; i < mtemp; ++i){
            GivensRotation giv_rot(R[index(i,i+1,m)],R[index(i+1,i+1,m)]);
            for (size_t j = i+1; j < m;++j){

                double r1 = R[index(i,j,m)];
                double r2 = R[index(i+1,j,m)];

                R[index(i,j,m)] = giv_rot.C()*r1 - giv_rot.S()*r2;
                R[index(i+1,j,m)] = giv_rot.S()*r1 + giv_rot.C()*r2;
            }
            giv_rot.apply_givens_to_cols(i,i+1,Y);
        }
        VecDoub Rtemp(mtemp*m/2);
        for (size_t i = 0; i < mtemp; ++i){
            for (size_t j = i; j < k; ++j){
                Rtemp[index(i,j,mtemp)] = R[index(i,j,m)];
            }
            for (size_t j = k+1; j < m; ++j){
                Rtemp[index(i,j-1,mtemp)] = R[index(i,j,m)];
            }
        }
        R = Rtemp;
        Y = Q(slice(0,nn,1),slice(0,mtemp,1));
        Z = Q(slice(0,nn,1),slice(mtemp,nn-mtemp,1));
    }
    void computeP(const MatDoub &G,const VecDoub &gk,
                  const MatDoub2DRef &Z, VecDoub &p){
           //size_t m = Y.cols();
           //size_t n = Y.rows();
           //VecDoub cz = Y*p_y;
           //cz = G*cz;
           //cz = cz*Z;
           VecDoub cz = gk*Z; //Z_t*c = c*Z

           size_t mm = cz.size();

           VecDoub p_z(mm,0.0);

           cg_redsys(Z,G,cz,p_z);

           //p = Y*p_y;
           p = Z*p_z;
    }
    void determineLagrangeMult(const MatDoub &Q, const VecDoub &R, const VecDoub &gk, VecDoub &lamda){
        VecDoub temp = gk*Q;
        size_t m = lamda.size();
        for (long int i = m-1; i >= 0; --i){
            double sum = 0.0;
            for (long int j = m-1; j > i; --j){
                sum += lamda[j]*R[index(i,j,m)];
            }
            lamda[i] = (temp[i] - sum)/R[index(i,i,m)];
        }

    }
public:
    VecDoub lamda;
    qp_lin_convex_act_set() = default;

    void qp_lin_convex_act_set_qrdcmp(const MatDoub &G, const VecDoub &c, const MatDoub &Aeq,
                                      const VecDoub &beq,const MatDoub &Ainq,
                                      const VecDoub &binq,VecDoub &x){
        //TODO: PHASE I encontrar un punto de inicio valido
        size_t meq = Aeq.rows();
        size_t minq = Ainq.rows();
        size_t n = x.size();
        lamda.reserve(meq+minq);
        set<size_t> W;   //Working set
        set<size_t> Win; //Inactive set;
        computeW0(W,Win,Ainq,binq,x);
        size_t minq_act = W.size();
        size_t m = meq + minq_act;
        //QR dcmp solo una ves
        qrdcmp qr(Aeq,Ainq,W,n);
        qr.householder_qr();
        //EXTRAYENDO Y y Z como slices en qr.Q
        MatDoub2DRef Y = qr.Q(slice(0,n,1),slice(0,m,1));
        MatDoub2DRef Z = qr.Q(slice(0,n,1),slice(m,n-m,1));

        VecDoub R(m*(m+1)/2); //Vector para almacenar la matriz qr.R aprovechando la forma de esta
        //Almacenando qr.R en R
        for (size_t i = 0; i < m; ++i)
            for (size_t j = i; j < m; ++j)
                R[index(i,j,m)] = qr.R(i,j);
        qr.R.clear(); //liberando qr.R

        //VecDoub p_y(m);
        //determinePy(R,meq,minq_act,p_y);
        VecDoub p(x.size(),0.0);
        bool calc_gk = true;
        VecDoub gk;
        for (size_t its = 0;;++its){
            if (calc_gk){
              gk = G*x;
              gk+=c;
              calc_gk = false;
            }

            if (m < n)
                computeP(G,gk,Z,p);
            else
               p.assign(n,0.0);

            if (p.norm() < 1.0e-8){ //Calcular lagrange mult
                lamda.resize(m);
                determineLagrangeMult(qr.Q,R,gk,lamda);
                double lamda_min = 0.0;
                auto iter = W.begin();
                size_t min_pos = meq;
                if (lamda.size() > meq){
                    lamda_min = lamda[meq];
                    auto iter_aux = W.begin();
                    for (size_t i = meq+1; i < lamda.size(); ++i){
                        ++iter_aux;
                        if (lamda[i] < lamda_min){
                            min_pos = i;
                            lamda_min = lamda[i];
                            iter = iter_aux;//Posicion de la minima lamda
                        }
                    }
                }
                if (lamda_min >= 0.0) return; //Se cumplen KKT cond
                Win.insert(*iter); //Agrega la restriccion que se desactivara al set inactivo
                W.erase(iter);//Elimina la restricion con menor lamda
                updateFactorizationOnConstrainDeactivation(min_pos,qr.Q,Y,Z,R);
                --m;
                //determinePyOnConstraintDeactivation(min_pos,R,p_y);
            }else{
                double ak = 1.0;
                auto iter_rest = Win.begin(); //Posicion de la posible restriccion a activar
                size_t constr_to_active = 0; //Index de la posible restriccion a activar
                for (auto iter = Win.begin();iter!=Win.end();++iter){
                    double den = (Ainq.row(*iter))*p;
                    if (den < 0.0){
                        double aktemp = binq[*iter] - Ainq.row(*iter)*x;
                        aktemp/=den;
                        if (aktemp < ak){
                            ak = aktemp;
                            iter_rest = iter;
                            constr_to_active = *iter;
                        }
                    }
                }
                p*=ak; //ap
                x+=p; //new x = x + ap
                calc_gk = true;
                if (ak < 1.0){ //Existen restricciones de bloqueo
                    W.insert(*iter_rest);
                    Win.erase(iter_rest);
                    updateFactorizationOnConstrainActivation(Ainq,constr_to_active,qr.Q,Y,Z,R);
                    ++m;
                    //determinePyOnConstraintActivation(R,p_y);
                }
            }
        }

    }
};

#endif // QP_LIN_CONVEX_ACT_SET_H
