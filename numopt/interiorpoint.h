#ifndef INTERIORPOINT_H
#define INTERIORPOINT_H

#include "linalg/lbldcmp.h"
#include "linalg/ludcmp.h"
#include "linalg/matrix_ref2d.cpp"
#include "linalg/matrix_ref1d.cpp"
#include <iostream>

using std::cout;
using std::endl;
using std::log;
//Interfaz para las resticciones de igualdad
struct constraints_cont
{
    double (*pFunc)(const VecDoub &x); //Puntero a funcion
    void (*pFuncD1)(const VecDoub &x,VecDoub &g); //Puntero a funcion que calcula en gradiente
    void (*pFuncD2)(const VecDoub &x,MatDoub &H); //Puntero a funcion que calcula el hessiano

    constraints_cont():pFunc(NULL),pFuncD1(NULL),pFuncD2(NULL){}
    constraints_cont(double func(const VecDoub&), void func1D(const VecDoub&,VecDoub&),void func2D(const VecDoub&,MatDoub&))
        :pFunc(func),pFuncD1(func1D),pFuncD2(func2D){}
};
//Interfaz para las resticciones de desigualdad
/*struct inq_const_cont
{
    double (*pFunc)(const VecDoub &x); //Puntero a funcion
    void (*pFuncD1)(const VecDoub &x,VecDoub &g); //Puntero a funcion que calcula en gradiente
    void (*pFuncD2)(const VecDoub &x,MatDoub &H); //Puntero a funcion que calcula el hessiano

    inq_const_cont(double func(const VecDoub&), void func1D(const VecDoub&,VecDoub&),void func2D(const VecDoub&,MatDoub&))
        :pFunc(func),pFuncD1(func1D),pFuncD2(func2D){}
};*/
template<typename T>
class interior_point{
private:
    T &funcd; //Funtor, funcion objetivo, gradiente y hessiano
    vector<constraints_cont> &v_ec; //direccion de las restricciones de igualdad
    vector<constraints_cont> &v_ic; //direccion de las restricciones de desigualdad
    double BETA; //Parametros para algoritmo de regularizacion BETA < 1.0
    double NU;   //Parametros para algoritmo de regularizacion NU > 0
    double DELTA;//delta parametro para control de inercia
    double MU; //Parametro de barrera
    double TAU;
    double V; //parametro de barrera, funcion de merito
    double alfa_s_max;
    double alfa_z_max;
    bool RECALCULATE_MERIT_FUN_FLAG;
    bool FIRST_TIME_MERIT_FUN_FLAG;
    double mfx;

    size_t l; //numero de resticciones de igualdad
    size_t m; //numero de restriciones de desigualdad

    //VecDoub xx;
    VecDoub ss;
    VecDoub yy;
    VecDoub zz;
    size_t n; //dimensiones del problema inicial
    MatDoub KKT;
    VecDoub p; //pasos de las variables
    matrix_inertia inert_ref;

    VecDoub Ce; //valores de las resiricciones de igualdad en x
    VecDoub Ci; //valores de las restricciones de desigualdad en x
    VecDoub g; //gradiente de la funcion objetivo
    double Dfx; //Valor de la derivada direccional de f en el punto x
    double Dfs; //Valor de la derivada direccional del termino de barrera sum[log(si)] i=1...m
    double Dfm; //Valor de la derivada direccional de la funcion de merito
    MatDoub Aet; //Matriz jacobiana transpueta de restricciones de igualdad
    MatDoub Ait; //Matriz jacobiana transpuesta de restricciones de desigualdad
    MatDoub1DRef px; //Referencias a partes del vector p
    MatDoub1DRef ps;
    MatDoub1DRef py;
    MatDoub1DRef pz;
    MatDoub Lxx;

    VecDoub rhs1;
    //VecDoub rhs2;
    //VecDoub rhs3;
    VecDoub rhs4;

    size_t its;



    struct meritFunction
    {
        meritFunction() {}

        double operator()(interior_point *p,const VecDoub &x,const VecDoub &s){
            double fx = p->funcd(x);
            double mu = p->MU; //parametro de barrera
            double v = p->V;
            VecDoub Ce; //Valores de las funciones de igualdad
            p->calculateEqConstFuncVal(x,Ce);
            VecDoub Ci; //Valores de las funciones de desigualdad
            p->calculateInqConstFuncVal(x,Ci);

            double l_ssum = 0.0;
            for (size_t i = 0; i < s.size(); ++i) l_ssum += log(s[i]);

            double norm_Ce = Ce.norm_l1();
            double norm_Ci_s = (Ci - s).norm_l1();

            double res = fx - mu*l_ssum + v*(norm_Ce + norm_Ci_s);
            return res;
        }
        double Dmf(interior_point *p){ //Derivada direccional
            double Dfx = p->Dfx;
            double Dfs = p->Dfs;
            double mu = p->MU; //parametro de barrera
            double v = p->V;

            VecDoub temp = p->px*p->Aet;
            double res = Dfx - mu*Dfs;
            for (size_t i = 0; i < p->l; ++i){
                if (p->Ce[i] >= 0.0) res += v*temp[i];
                else res -= v*temp[i];
            }
            temp = p->px*p->Ait;
            temp-=p->ps;

            VecDoub Ci_s = p->Ci - p->ss;
            for (size_t i = 0; i < p->m; ++i){
                if (Ci_s[i] >= 0.0) res += v*temp[i];
                else res -= v*temp[i];
            }
            return res;
        }
    };


    meritFunction mfunc;

    void modifySystemForInertiaControl(MatDoub &K,const double &GAMMA){
        for (size_t i = 0; i < n; ++i){
            K(i,i)+=DELTA;
        }
        if (GAMMA != 0.0){
            MatDoub2DRef ref = K(slice(n+m,l,1),slice(n+m,l,1));
            for (size_t i = 0; i < l ; ++i){
                ref(i,i) = -GAMMA;
            }
        }
    }
    bool factorizeSystemInertiaCheckAndSolve(MatDoub &K,VecDoub &pp,matrix_inertia &inert){
        //double GAMMA = 0.0;
        lbldcmp lbl(K);
        inert = lbl.inertia();
        if (inert_ref == inert){
            //TODO usar la factorizacion para esto, hay que ver como lo hago
            LUdcmp lu(K);
            lu.solve(pp,pp);
            return true;
        }

        return false;
    }

    void inercorr_reg(MatDoub &K,VecDoub &pp){
        matrix_inertia inert;
        double GAMMA = 0.0;
        if (factorizeSystemInertiaCheckAndSolve(K,pp,inert)) return;

        if (inert.c_eigen != 0) GAMMA = 1e-9*NU*pow(MU,BETA);

        if (DELTA == 0.0){
            DELTA = 1e-5;
        }else{
            DELTA/=2.0;
        }
        modifySystemForInertiaControl(K,GAMMA);
        //cout << KKT << endl;
        while (!(factorizeSystemInertiaCheckAndSolve(K,pp,inert))){
            DELTA*=10;
            modifySystemForInertiaControl(K,0.0);
            //cout << KKT << endl;
        }
    }
    void calculateInqConstFuncVal(const VecDoub &x,VecDoub &ci){
        ci.resize(m);
        for (size_t i = 0; i < m; ++i) ci[i] = v_ic[i].pFunc(x);
    }
    void calculateEqConstFuncVal(const VecDoub &x,VecDoub &ce){
        ce.resize(l);
        for (size_t i = 0; i < l; ++i) ce[i] = v_ec[i].pFunc(x);
    }
    void determineInitSlacks(const VecDoub &ci){
        for (size_t i = 0; i < m; ++i){
            if (ci[i] > 0.0) ss[i] = ci[i];
            else ss[i] = 1.0; //TODO hay que ver esto porque se x seria un punto que viola las restricciones
        }
    }
    /*void determineInitZ(){
        for (size_t i = 0; i < m; ++i) zz[i] = MU/ss[i];
    }*/
    //Aet transpuesta de Ae
    void determineAet(const VecDoub &x,MatDoub &Aet){
        Aet.resize(n,l);
        VecDoub g;
        for (size_t i = 0; i < l; ++i){
            v_ec[i].pFuncD1(x,g);
            Aet.column(i) = g;
        }
    }
    void determineAit(const VecDoub &x,MatDoub &Ait){
        Ait.resize(n,m);
        VecDoub g;
        for (size_t i = 0; i < m; ++i){
            v_ic[i].pFuncD1(x,g);
            Ait.column(i) = g;
        }
    }
    void determineInitLagMult(){
        size_t nn = l+m;
        MatDoub Atemp(nn,nn);

        MatDoub2DRef mref = Atemp(slice(0,l,1),slice(0,l,1));

        MatDoub A = Aet.transpose();
        //MatDoub Ai = Ait.transpose();

        mref = A*Aet;

        mref = Atemp(slice(0,l,1),slice(l,m,1));
        mref = A*Ait;

        A = Ait.transpose();

        mref = Atemp(slice(l,m,1),slice(0,l,1));
        mref = A*Aet;


        mref = Atemp(slice(l,m,1),slice(l,m,1));
        mref = A*Ait;

        for (size_t i = 0; i < m; ++i){
            mref(i,i)+= ss[i]*ss[i];
        }

        VecDoub b;
        b.reserve(nn);
        b = g*Aet;
        VecDoub btemp = g*Ait;
        btemp+=MU*ss;
        b.insert(b.end(),btemp.begin(),btemp.end());

        LUdcmp lu(Atemp);
        lu.solve(b,b);

        for (size_t i = 0; i < l; ++i){
            yy[i] = b[i];
        }
        for (size_t i = l; i < nn; ++i){
            if (b[i] >= 0.0){
                zz[i-l] = b[i];
            }else{
                zz[i-l] = std::min(1e-4,MU/ss[i-l]);
            }
        }

    }
    double meritFuncX(const double &ce_norm,const double &ci_s_norm,const VecDoub &x){
        double fx = funcd(x);
        double t1 = 0.0;
        for (size_t i = 0; i < m; ++i) t1 += std::log(ss[i]);
        t1*=MU;

        FIRST_TIME_MERIT_FUN_FLAG = false;

        return fx - t1 + V*(ce_norm + ci_s_norm);

    }
    double meritFuncXnew(const VecDoub &ce_new,const VecDoub &ci_new, const VecDoub &x_new,const VecDoub &s_new){
        double fx_new = funcd(x_new);
        double t1 = 0.0;
        for (size_t i = 0; i < m; ++i) t1 += std::log(s_new[i]);
        t1*=MU;

        return fx_new - t1 + V*(ce_new.norm_l1() + (ci_new - s_new).norm_l1());

    }
    void backtrackingLinSearch(VecDoub &x){
        //MatDoub1DRef px = p(0,n,1);
        //MatDoub1DRef ps = p(n,m,1);
        double alfa = alfa_s_max;
        VecDoub xtemp = x + alfa*px;
        VecDoub stemp = ss + alfa*ps;

        //VecDoub ce_temp(l);
        //calculateEqConstFuncVal(xtemp,ce_temp);
        //VecDoub ci_temp(m);
        //calculateInqConstFuncVal(xtemp,ci_temp);

        //double fx_new = meritFuncXnew(ce_temp,ci_temp,xtemp,stemp);
        double mfx_new = mfunc(this,xtemp,stemp);

        while (mfx_new > mfx + 0.0001*alfa*Dfm){
            //alfa = - Dfx*alfa*alfa/(2.0*(fx_new - mfx - Dfx*alfa));
            alfa*=0.9;
            xtemp = x + alfa*px;
            stemp = ss + alfa*ps;
            //calculateEqConstFuncVal(xtemp,ce_temp);
            //calculateInqConstFuncVal(xtemp,ci_temp);

            //mfx_new = meritFuncXnew(ce_temp,ci_temp,xtemp,stemp);
            mfx_new = mfunc(this,xtemp,stemp);
        }
        alfa_s_max = alfa;
        x = xtemp;
        ss = stemp;
        mfx = mfx_new;

    }
    double meritFuncDirDeriv(const double &gp, const double &sp,const double &norm1,const double norm2){
        return gp - sp + V*(norm1 + norm2);
    }
    void determineStepL(VecDoub &x){
        //MatDoub1DRef px = p(0,n,1);
        //MatDoub1DRef ps = p(n,m,1);

        VecDoub pAet = px*Aet;
        //double pAet_norm = pAet.norm_l1();
        double mp = (pAet + Ce).norm_l1();
        //VecDoub ci_s = ci-ss; //ci(x) - s
        VecDoub pAit_ps = px*Ait - ps;
        //double pAit_ps_norm = pAit_ps.norm_l1();
        mp += (pAit_ps + rhs4).norm_l1();

        double ce_norm = Ce.norm_l1();
        double ci_s_norm = rhs4.norm_l1();
        double m0 = ce_norm + ci_s_norm;

        Dfx = g*px; //grad(f)*px
        double pLp = 0.5*px*Lxx*px;
        Dfs = 0.0;
        for (size_t i = 0 ; i < m; ++i){
            Dfs += ps[i]/ss[i];
        }
        //Dfs*=MU;

        MatDoub2DRef S_1Z = KKT(slice(n,m,1),slice(n,m,1));
        double pS_1Zp = 0.0;
        for (size_t i = 0; i < m; ++i){
            pS_1Zp += ps[i]*S_1Z(i,i)*ps[i];
        }
        pS_1Zp*=0.5;

        double nom = Dfx + pLp + MU*Dfs + pS_1Zp;

        determineV(nom,m0);

        //Dfx = meritFuncDirDeriv(gp,sp,pAet_norm,pAit_ps_norm);
        Dfm = mfunc.Dmf(this);

        /*if (RECALCULATE_MERIT_FUN_FLAG || FIRST_TIME_MERIT_FUN_FLAG){
            mfx = meritFuncX(ce_norm,ci_s_norm,x);
        }*/
        mfx = mfunc(this,x,ss);
        backtrackingLinSearch(x);

    }
    //Parametro de funcion de merito
    void determineV(const double &nom, const double &m0,const double RHO = 0.5){
        /*double v = nom/((1-RHO)*m0);
        if (v > V){
            V = 1.5*v;
            RECALCULATE_MERIT_FUN_FLAG = true;
        }else{
            RECALCULATE_MERIT_FUN_FLAG = false;
        }*/
        double v = 0.0;
        for (size_t i = 0; i < l; ++i){
            double temp = std::abs(yy[i]);
            if (v < temp) v = temp;
        }
        for (size_t i = 0; i < m; ++i){
            double temp = std::abs(zz[i]);
            if (v < temp) v = temp;
        }

        if (V <= v){
            v*=1.5;
            V = v;
        }

    }
    double calcE(bool FLAG){
        VecDoub sz(m);
        if (FLAG)
            for (size_t i = 0; i < m; ++i) sz[i] = ss[i]*zz[i] - MU;
        else
            for (size_t i = 0; i < m; ++i) sz[i] = ss[i]*zz[i];

        double max = rhs1.norm_l1();
        double temp = sz.norm_l1();
        if (max < temp){
            max = temp;
        }
        temp = Ce.norm_l1();
        if (max < temp){
            max = temp;
        }
        temp = rhs4.norm_l1();
        if (max < temp){
            max = temp;
        }
        return max;
    }
    void determineSz(VecDoub &v2,bool FLAG){
        if (FLAG)
            for (size_t i = 0; i < m; ++i) v2[i] = ss[i]*zz[i] - MU;
        else
            for (size_t i = 0; i < m; ++i) v2[i] = ss[i]*zz[i];
    }
    void conformRHS(){
        for (size_t i = 0; i < n; ++i) p[i] = -rhs1[i];
        for (size_t i = 0; i < m; ++i) p[i+n] = -(zz[i] - MU/ss[i]);
        for (size_t i = 0; i < l; ++i) p[i+n+m] = -Ce[i];
        for (size_t i = 0; i < m; ++i) p[i+n+m+l] = -rhs4[i];

        //cout << p << endl;
    }
    void conformKKTMatrix(const VecDoub &x){
        MatDoub2DRef mref = KKT(slice(0,n,1),slice(0,n,1));
        //MatDoub Htemp;
        funcd.df2(x,Lxx);

        for (size_t i = 0; i < l; ++i){
            MatDoub Htemp1;
            v_ec[i].pFuncD2(x,Htemp1);
            Htemp1*=yy[i];
            Lxx -= Htemp1;
        }
        for (size_t j = 0; j < m; ++j){
            MatDoub Htemp1;
            v_ic[j].pFuncD2(x,Htemp1);
            Htemp1*=zz[j];
            Lxx -= Htemp1;
        }
        mref = Lxx;

        mref = KKT(slice(0,n,1),slice(n+m,l,1));
        mref = Aet;
        mref = KKT(slice(0,n,1),slice(n+m+l,m,1));
        mref = Ait;

        mref = KKT(slice(n,m,1),slice(n,m,1));
        for (size_t i = 0 ; i < m; ++i){
            mref(i,i) = zz[i]/ss[i];
        }
        mref = KKT(slice(n,m,1),slice(n+m+l,m,1));
        for (size_t i = 0 ; i < m; ++i){
            mref(i,i) = -1.0;
        }
        mref = KKT(slice(n+m,l,1),slice(0,n,1));

        for (size_t i = 0 ; i < l; ++i){
            for (size_t j = 0; j < n; ++j){
                mref(i,j) = Aet(j,i);
            }
        }

        mref = KKT(slice(n+m+l,m,1),slice(0,n,1));

        for (size_t i = 0 ; i < m; ++i){
            for (size_t j = 0; j < n; ++j){
                mref(i,j) = Ait(j,i);
            }
        }
        mref = KKT(slice(n+m+l,m,1),slice(n,m,1));
        for (size_t i = 0 ; i < m; ++i){
            mref(i,i) = -1.0;
        }
        //cout << KKT << endl;
    }
    void determineAlfasMax(){

        /*double alfa = 1.0;

        if (m == 0){
            alfa_s_max = 1.0;
            alfa_z_max = 1.0;
            return;
        }*/
        for (size_t i = 0; i < m; ++i){
            double atemp_s = -TAU*ss[i]/p[n+i]; //ps p[n:n+m)
            double atemp_z = -TAU*zz[i]/p[n+m+l+i]; //pz p[n+m+l:n+m+l+m)

            if (atemp_s < alfa_s_max && atemp_s > 0.0) alfa_s_max = atemp_s;
            if (atemp_z < alfa_z_max && atemp_z > 0.0) alfa_z_max = atemp_z;
        }
        //if (alfa_s_max == 0.0) alfa_s_max = 1.0;
        //if (alfa_z_max == 0.0) alfa_z_max = 1.0;
    }
public:
    interior_point(T &func,vector<constraints_cont> &eqc,vector<constraints_cont> &inqc,const double nu = 0.5,const double beta = 0.5,const double tau = 0.995):
        funcd(func),v_ec(eqc),v_ic(inqc),NU(nu),BETA(beta),DELTA(0.0),MU(2.0),TAU(tau),V(1.0),l(v_ec.size()),
        m(v_ic.size()),RECALCULATE_MERIT_FUN_FLAG(true),FIRST_TIME_MERIT_FUN_FLAG(true){}

    void solve(VecDoub &x,const double toll = 1e-8){
        its = 0;
        double mu_tol = MU;
        n = x.size();
        inert_ref.p_eigen = n+m;
        inert_ref.n_eigen = l+m;
        inert_ref.c_eigen = 0;
        ss.resize(m); //Slacks
        yy.resize(l); //Mult Lag igualdad
        zz.resize(m); //Mult Lag des
        size_t nn = n+m+l+m;
        p.resize(nn);
        px = p(0,n,1);
        ps = p(n,m,1);
        py = p(n+m,l,1);
        pz = p(n+m+l,m,1);

        //TODO Arreglar todo eso, hay miembros de la clase para esto
        //Ci(x)
        //VecDoub ci(m);
        calculateInqConstFuncVal(x,Ci);

        determineInitSlacks(Ci);

        //TODO Arreglar todo esto, hay miembros de la clase para esto
        //VecDoub g(n);  //Gradiente de la funcion objetivo
        funcd.df(x,g); //Calcula el gradiente de la funcion objetivo

        //MatDoub Aet;
        determineAet(x,Aet);
        //MatDoub Ait;
        determineAit(x,Ait);

        determineInitLagMult();

        //Ce(x)
        //VecDoub ce(l);
        calculateEqConstFuncVal(x,Ce);

        rhs1 = g - Aet*yy - Ait*zz;
        //VecDoub v2(m);
        //determineSz(v2,false);
        rhs4 = Ci - ss;
        double E = calcE(0);

        //determineSz(v2,true);
        double Emu = calcE(1);

        KKT.resize(nn,nn);
        KKT.setToZero();

        //MatDoub1DRef py = p(n+m,l,1);
        //MatDoub1DRef pz = p(n+m+l,m,1);

        while (E > toll){
            while (Emu > mu_tol){
                conformRHS();
                conformKKTMatrix(x);
                inercorr_reg(KKT,p);
                for (size_t i = n+m; i < nn; ++i) p[i] = -p[i];
                alfa_s_max = 1.0;
                alfa_z_max = 1.0;
                determineAlfasMax();
                //determineStepL(x);
                x += alfa_s_max*px;
                ss += alfa_s_max*ps;
                yy += alfa_z_max*py;
                zz += alfa_z_max*pz;

                //updates
                funcd.df(x,g);
                determineAet(x,Aet);
                determineAit(x,Ait);
                calculateEqConstFuncVal(x,Ce);
                calculateInqConstFuncVal(x,Ci);
                rhs1 = g - Aet*yy - Ait*zz;
                //determineSz(v2,false);
                rhs4 = Ci - ss;
                E = calcE(0);
                //determineSz(v2,true);
                Emu = calcE(1);
                ++its;

                //cout << x << endl;
                //cout << Emu << endl;
                //cout << E << endl;

            }
            MU*=0.1; //TODO: cambiar esto por otra estrategia
            mu_tol = MU;
        }
    }
    size_t IterNumber() const{return its;}

};


#endif // INTERIORPOINT_H
