#ifndef LBLDCMP_H
#define LBLDCMP_H

#include "mathmatrix.h"
#include "gauss2x2.h"

struct index_pair{
    size_t i;
    size_t j;

    index_pair():i(0),j(0){}
    index_pair(const size_t &ii, const size_t &jj):i(ii),j(jj){}
};

//Inercia de una matrix numero de autovalores positivos, negativos y 0
struct matrix_inertia{
    size_t p_eigen;
    size_t n_eigen;
    size_t c_eigen;

    matrix_inertia():p_eigen(0),n_eigen(0),c_eigen(0){}
    matrix_inertia(const size_t &p,const size_t &n,const size_t &c):
        p_eigen(p),n_eigen(n),c_eigen(c){}
};
inline bool operator==(const matrix_inertia &in1,const matrix_inertia &in2){
    return (in1.p_eigen == in2.p_eigen) && (in1.n_eigen == in2.n_eigen) && (in1.c_eigen == in2.c_eigen);
}
inline bool operator!=(const matrix_inertia &in1,const matrix_inertia &in2){
    return !(in1==in2);
}

/*Descomposicion PAP^t = LBL^t*/
class lbldcmp{
private:
    size_t nn;
    vector<size_t> P; //Para las permutaciones
    VecDoub L; //Vector donde se almacena la matriz triangular unitaria L
    vector<MatDoub> B; //Vector donde se almacenan los bloques de la matriz B 
    matrix_inertia inert;
    double mu0;
    double mu1;
    index_pair mu0_indxs;
    size_t mu1_indx;

    void setPermVec(vector<size_t> &p){
        for (size_t i = 0; i < p.size(); ++i) p[i] = i;
    }

    size_t index(const size_t &i, const size_t &j)
    {
        size_t ind = (i <= j) ? (i*nn + j - i*(i+1)/2) : (j*nn + i - j*(j+1)/2);

        return ind;
    }

    void determinePiv(const MatDoub &A){

        mu0 = 0.0;
        mu1 = 0.0;

        size_t n = A.rows();

        for (size_t i = 0; i < n; ++i){
            double temp = abs(A(i,i));
            if (mu1 < temp){
                mu1 = temp;
                mu1_indx = i;
            }
            for (size_t j = i; j < n; ++j){
                temp = abs(A(i,j));
                if (mu0 < temp){
                    mu0 = temp;
                    mu0_indxs.i = i;
                    mu0_indxs.j = j;
                }
            }
        }
    }
    void swapPermVecPos(const size_t &indx, const size_t &n,vector<size_t> &p){
        size_t k = indx+n;
        for(size_t i = k; i > n; --i) std::swap(p[i],p[i-1]);
    }
public:
    lbldcmp(const MatDoub &A):nn(A.rows()),P(nn),L(nn*(nn+1)/2){
        MatDoub Atemp(A);
        //size_t n = Atemp.rows();
        //MatDoub Atemp1(nn,nn);
        B.reserve(nn);
        const double alfa = (1+sqrt(17))/8;
        setPermVec(P);
        vector<size_t> Ptemp(P);

        size_t m = Atemp.rows();
        size_t k = 0;

        while (m > 0){
            MatDoub E(2,2);

            determinePiv(Atemp);

            if (mu1 >= alfa*mu0){

                //size_t n1 = Atemp1.rows()-1;
                //Atemp1.resize(n1,n1);
                E.resize(1,1);
                double e = Atemp(mu1_indx,mu1_indx);
                E(0,0) = e;
                swapPermVecPos(mu1_indx,k,P);
                swapPermVecPos(mu1_indx,0,Ptemp);
                B.push_back(E);

                L(index(k,k)) = 1.0;

                VecDoub C(m-1);

                for (size_t i = 1; i < m; ++i){
                    C[i-1] = Atemp(mu1_indx,Ptemp[i]);
                    //C(i-2,1) = Atemp(mu0_indxs.j,Ptemp[i]);
                    L[index(k,k+i)] = C[i-1]/e;
                    //L[index(k+1,k+i)] = Atemp(mu0_indxs.i,Ptemp[i])*Einv(1,0) + Atemp(mu0_indxs.j,Ptemp[i])*Einv(1,1);
                }

                for (size_t i = 1; i < m; ++i){
                    for (size_t j = 1; j < m; ++j){
                       Atemp(i-1,j-1) = Atemp(Ptemp[i],Ptemp[j]) - L[index(k,k+i)]*C[j-1];
                    }
                }
                MatDoub2DRef mref = Atemp(slice(0,m-1,1),slice(0,m-1,1));
                Atemp = mref;
                //Atemp.resize(m-1,m-1);
                Ptemp.resize(m-1);
                setPermVec(Ptemp);

                m = Atemp.rows();
                k = nn - m;

                if (e > 0.0){
                    ++(inert.p_eigen);
                    continue;
                }
                if (e < 0.0){
                    ++(inert.n_eigen);
                    continue;
                }
                if (e == 0.0){
                    ++(inert.c_eigen);
                    continue;
                }

            }else{
                E(0,0) = Atemp(mu0_indxs.i,mu0_indxs.i);
                E(0,1) = Atemp(mu0_indxs.i,mu0_indxs.j);
                E(1,0) = Atemp(mu0_indxs.i,mu0_indxs.j);
                E(1,1) = Atemp(mu0_indxs.j,mu0_indxs.j);
                B.push_back(E);
                MatDoub Einv(2,2);
                gauss2x2 g2x2;
                g2x2.inverse(E,Einv);

                //j > i orden importante
                //TODO acer una funcion que permute con 2 valores
                swapPermVecPos(mu0_indxs.j,k,P);
                swapPermVecPos(mu0_indxs.i+1,k,P);
                swapPermVecPos(mu0_indxs.j,0,Ptemp);
                swapPermVecPos(mu0_indxs.i+1,0,Ptemp);

                L[index(k,k)] = 1.0;
                L[index(k,k+1)] = 0.0;
                L[index(k+1,k+1)] = 1.0;

                MatDoub C(m-2,2);
                for (size_t i = 2; i < m; ++i){
                    double tmp1 = Atemp(mu0_indxs.i,Ptemp[i]);
                    double tmp2 = Atemp(mu0_indxs.j,Ptemp[i]);
                    C(i-2,0) = tmp1;
                    C(i-2,1) = tmp2;
                    L[index(k,k+i)] = tmp1*Einv(0,0) + tmp2*Einv(1,0);
                    L[index(k+1,k+i)] = tmp1*Einv(0,1) + tmp2*Einv(1,1);
                }

                for (size_t i = 2; i < m; ++i){
                    for (size_t j = 2; j < m; ++j){
                       Atemp(i-2,j-2) = Atemp(Ptemp[i],Ptemp[j]) - L[index(k,k+i)]*C(j-2,0) -  L[index(k+1,k+i)]*C(j-2,1);
                    }
                }
                MatDoub2DRef mref = Atemp(slice(0,m-2,1),slice(0,m-2,1));
                Atemp = mref;
                Ptemp.resize(m-2);
                setPermVec(Ptemp);

                ++(inert.p_eigen);
                ++(inert.n_eigen);

                m = Atemp.rows();
                k = nn - m;
            }

        }

    }
    ~lbldcmp(){}

    matrix_inertia inertia() const{return inert;}

};



#endif // LBLDCMP_H
