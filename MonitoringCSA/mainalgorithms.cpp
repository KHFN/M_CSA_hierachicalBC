#include <mcl/bls12_381.hpp>
#include <iostream>
#include <limits>
#include <cmath>
#include <vector>
#include <malloc.h>
#include "primitiveroots.hpp"
#include "ntt.hpp"
#include "polyonfr.hpp"

using namespace mcl;
using namespace mcl::bn;
using namespace std;
struct pcs;



template <typename T>
void putG(vector<T> a){
    for(uint32_t i=0;i<a.size();i++){
        cout<<"<"<<a[i].getStr(10).c_str()<<">"<<endl;
    }
}

GT pairing_eq(G1 a, G2 b){

    GT res;
    pairing(res, a, b);
    return res;

}

struct IDset{
    Fr* idset;
    FrVec idVec;
    void init(uint32_t n)
    {
        FrVec ids;
        ids.resize(n,0);
   
        for(int i=0;i<n;i++){
            Fr s;
            s.setByCSPRNG();
            ids[i]=s;
        }
        idVec = ids;
        idset = &ids[0];
    }
    
};

struct PCS{

    
    struct pp{
        vector<G1> g1si;
        vector<G2> g2si;
        vector<G1> h1si;
        vector<G2> h2si;
    } pp;
    Fr s;

    void setup(uint32_t D){
        Fr s;
        s.setByCSPRNG();
        this->s = s;

        G1 g1;
        G2 g2;
        G1 h1;
        G2 h2;
        hashAndMapToG1(g1,"g1");
        hashAndMapToG2(g2,"g2");
        hashAndMapToG1(h1,"h1");
        hashAndMapToG2(h2,"h2");

        vector<G1> g1si(D+1);
        vector<G2> g2si(D+1);
        vector<G1> h1si(D+1);
        vector<G2> h2si(D+1);

        g1si[0]=g1;
        g2si[0]=g2;
        h1si[0]=h1;
        h2si[0]=h2;
    
        for(int i=1;i<=D;i++){
        
            g1si[i]=g1si[i-1]*s;
            g2si[i]=g2si[i-1]*s;
            h1si[i]=h1si[i-1]*s;
            h2si[i]=h2si[i-1]*s;

        }

        this->pp.g1si=g1si;
        this->pp.g2si=g2si;
        this->pp.h1si=h1si;
        this->pp.h2si=h2si;

    }

    G2 commit(FrVec p){

        G2 tmp;
        G2 res;
        res.clear();
        G2::mulVec(res, &pp.g2si[0], &p[0], p.size());

        return res;
    }
    G2 commit(Fr* p, uint32_t n){

        G2 tmp;
        G2 res;
        res.clear();
        G2::mulVec(res, &pp.g2si[0], p, n);

        return res;
    }

    G1 commit_G1(FrVec p){

        G1 tmp;
        G1 res;
        res.clear();
        G1::mulVec(res, &pp.g1si[0], &p[0], p.size());
        return res;
    }

    G1 commit_h1(FrVec p){

        G1 tmp;
        G1 res;
        res.clear();
        G1::mulVec(res, &pp.h1si[0], &p[0], p.size());
        return res;
    }

    G2 commit_h2(FrVec p){

        G2 tmp;
        G2 res;
        res.clear();
        G2::mulVec(res, &pp.h2si[0], &p[0], p.size());
        return res;
    }



    G2 Eval(G2 C, Fr z, Fr y, FrVec p){
        
        FrVec pz ={z};
        FrVec pmz=PolySub(p,pz);
        FrVec xmy={-y,1};
        FrVec w=PolyLongDiv(pmz,xmy);

        G2 pi_open=commit(w);  

        return pi_open;      
    }

    G1 Eval(G1 C, Fr z, Fr y, FrVec p){
        
        FrVec pz ={z};
        FrVec pmz=PolySub(p,pz);
        FrVec xmy={-y,1};
        FrVec w=PolyLongDiv(pmz,xmy);

        G1 pi_open=commit_G1(w);  

        return pi_open;      
    }



    bool verify(G2 C, Fr z, Fr y, G2 pi_open){
        bool flag=false;
        GT res1;
        GT res2;

        pairing(res1, this->pp.g1si[0],C+(this->pp.g2si[0]*(-z)));
        pairing(res2, (this->pp.g1si[1])+(this->pp.g1si[0]*(-y)),pi_open);

        if(res1==res2){
            flag = true;
        }

        return flag;       
    }

    bool verify(G1 C, Fr z, Fr y, G1 pi_open){
        bool flag=false;
        GT res1;
        GT res2;

        pairing(res1, C+(this->pp.g1si[0]*(-z)),this->pp.g2si[0]);
        pairing(res2, pi_open, (this->pp.g2si[1])+(this->pp.g2si[0]*(-y)));

        if(res1==res2){
            flag = true;
        }

        return flag;       
    }
};

struct PCS_h{

    
    struct pp{
        vector<G1> g1si;
        vector<G2> g2si;
    } pp;

    void setup(uint32_t D){
        Fr s;
        s.setByCSPRNG();

        G1 g1;
        G2 g2;
        hashAndMapToG1(g1,"h1");
        hashAndMapToG2(g2,"h2");

        vector<G1> g1si(D+1);
        vector<G2> g2si(D+1);

        g1si[0]=g1;
        g2si[0]=g2;
    
        for(int i=1;i<=D;i++){
        
            g1si[i]=g1si[i-1]*s;
            g2si[i]=g2si[i-1]*s;

        }

        this->pp.g1si=g1si;
        this->pp.g2si=g2si;

    }

    G2 commit(FrVec p){

        G2 tmp;
        G2 res;
        res.clear();
        G2::mulVec(res, &pp.g2si[0], &p[0], p.size());

        return res;
    }
    G2 commit(Fr* p, uint32_t n){

        G2 tmp;
        G2 res;
        res.clear();
        G2::mulVec(res, &pp.g2si[0], p, n);

        return res;
    }

    G1 commit_G1(FrVec p){

        G1 tmp;
        G1 res;
        res.clear();
        G1::mulVec(res, &pp.g1si[0], &p[0], p.size());
        return res;
    }

    G2 Eval(G2 C, Fr z, Fr y, FrVec p){
        
        FrVec pz ={z};
        FrVec pmz=PolySub(p,pz);
        FrVec xmy={-y,1};
        FrVec w=PolyLongDiv(pmz,xmy);

        G2 pi_open=commit(w);  

        return pi_open;      
    }

    G1 Eval(G1 C, Fr z, Fr y, FrVec p){
        
        FrVec pz ={z};
        FrVec pmz=PolySub(p,pz);
        FrVec xmy={-y,1};
        FrVec w=PolyLongDiv(pmz,xmy);

        G1 pi_open=commit_G1(w);  

        return pi_open;      
    }



    bool verify(G2 C, Fr z, Fr y, G2 pi_open){
        bool flag=false;
        GT res1;
        GT res2;

        pairing(res1, this->pp.g1si[0],C+(this->pp.g2si[0]*(-z)));
        pairing(res2, (this->pp.g1si[1])+(this->pp.g1si[0]*(-y)),pi_open);

        if(res1==res2){
            flag = true;
        }

        return flag;       
    }

    bool verify(G1 C, Fr z, Fr y, G1 pi_open){
        bool flag=false;
        GT res1;
        GT res2;

        pairing(res1, C+(this->pp.g1si[0]*(-z)),this->pp.g2si[0]);
        pairing(res2, pi_open, (this->pp.g2si[1])+(this->pp.g2si[0]*(-y)));

        if(res1==res2){
            flag = true;
        }

        return flag;       
    }
};


struct zkbpacc_setup {
    G1* g1si;
    G2* g2si;
    G1  Frakg;
    G1  Frakh;
    G2  h2;

    void init(uint32_t n)
    {   
        G1 g1;
        G2 g2;

    /*Getting generator*/
        hashAndMapToG1(g1,"g1");
        hashAndMapToG2(g2,"g2");
        hashAndMapToG1(Frakg,"Frakg");
        hashAndMapToG1(Frakh,"Frakh");
        hashAndMapToG2(h2,"h2");


        Fr s;
        s.setByCSPRNG();

        g1si=(G1*) malloc(sizeof(G1)*(n+1));
        g2si=(G2*) malloc(sizeof(G2)*(n+1));

        g1si[0]=g1;
        g2si[0]=g2;  
    
        for(int i=1;i<=n;i++){
        
            g1si[i]=g1si[i-1]*s;
            g2si[i]=g2si[i-1]*s;

        }
    }
};

struct zkbpacc_setup_pcs {
    
    PCS& pcs;
    G1  Frakg;
    G1  Frakh;
    G2  h2;

    zkbpacc_setup_pcs(PCS& _pcs) : pcs(_pcs) {
        init();
    }

    void init()
    {   
        hashAndMapToG1(Frakg,"Frakg");
        hashAndMapToG1(Frakh,"Frakh");
        hashAndMapToG2(h2,"h2");
    }

};

template <typename T>
T* Vec2arr(vector<T> a){

    uint32_t len=a.size();
    T* res=(T*) malloc((sizeof(T))*len);

    for(uint32_t i=0;i<len;i++){
        res[i]=a[i];
    }

    return res;
}

template <typename T>
vector<T> arr2Vec(T* a, uint32_t n){

    vector<T> dest(n);
    memcpy(&dest[0],&a[0],n*sizeof(T));

    return dest;
}

template <typename T>
vector<T> convG(vector<T> a, vector<T> b){

    if(a.size()!=b.size()){
        throw runtime_error("Two vectors must have same length");
    }

    uint32_t n = a.size();
    vector<T> res;
    for(uint32_t i=0;i<n;i++){

        T element = a[i]+b[i];
        res.insert(res.begin(),element);

    }
    return res;
    
}

/*
a in GG1^n, b in GG2
return \prod_{i=0}^{n-1} e(a[i],b[i]);
*/
GT MultiPairing(vector<G1> a, vector<G2> b){
    if(a.size()!=b.size()){
        throw runtime_error("Two vectors must have same length");
    }

    uint32_t n=a.size();
    G1* a_arr=&a[0]; G2* b_arr=&b[0];
    GT c;
    millerLoopVec(c,a_arr,b_arr,n);
    finalExp(c,c);

    return c;
}

struct IPPproof{

        vector<G2> gg;
        GT P;
        GT Q;
        vector<GT> L;
        vector<GT> R;
        vector<Fr> x;
        G1 w;
        G2 g;

        IPPproof(vector<G2> gg_0, GT P){
                this->gg=gg_0;
                this->P=P;
                this->L.clear();
                this->R.clear();
                this->x.clear();
        }

        void Prove(vector<G2> gg, GT P, vector<G1> ww){
            if((IsPowerOfTwo(ww.size())==false) || (IsPowerOfTwo(gg.size())==false)){
                throw runtime_error("gg and ww need to be power of two");
            }

            if(ww.size()==1){
                this->w = ww[0];
                this->g = gg[0];
                
            }else{

                vector<G2> gg_L(gg.begin(), gg.begin()+gg.size()/2);
                vector<G2> gg_R(gg.begin()+gg.size()/2, gg.end());
                
                vector<G1> ww_L(ww.begin(), ww.begin()+ww.size()/2);
                vector<G1> ww_R(ww.begin()+ww.size()/2, ww.end());
                
                G2* pt1=&gg_L[0];

                GT l = MultiPairing(ww_R,gg_L);
                GT r = MultiPairing(ww_L,gg_R);

                L.insert(L.begin(),l);
                R.insert(R.begin(),r);

                string buf = this->P.getStr()+l.getStr()+r.getStr();
                Fr c;
                c.setHashOf(buf);

                x.insert(x.begin(),1/c);
                vector<G1> ww_hat(ww.size()/2);
                vector<G2> gg_hat(gg.size()/2);

                
                for(uint32_t i=0;i<ww_R.size();i++){

                    ww_hat[i]=ww_L[i]+(ww_R[i]*(c));
                    gg_hat[i]=gg_L[i]+(gg_R[i]*(1/c));
                    
                }
                
                GT Lx;
                GT::pow(Lx,l,c);
                GT Rinvx;
                GT::pow(Rinvx,r,(1/c));

                GT P_hat=Lx*P*Rinvx;

                Prove(gg_hat,P_hat,ww_hat);
                
            }
        }   
};

bool IPPverify(IPPproof pi){

    uint32_t l = pi.L.size();
    FrVec x;

    GT P_v=pi.P;
    GT Lx;
    GT Rxinv;

    for(uint32_t i=0;i<l;i++){
        string buf = pi.P.getStr()+pi.L[i].getStr()+pi.R[i].getStr();
        Fr c;
        c.setHashOf(buf);
        x.push_back(1/c);
        GT::pow(Lx,pi.L[i],c);
        GT::pow(Rxinv,pi.R[i],1/c);
        P_v=Lx*P_v*Rxinv;
    }

    uint32_t n=pi.gg.size();
    uint32_t n_l=nextPowOf2(log2(n));
    vector<FrVec> M(n_l, FrVec((n>>1)+1,0));
    for(uint32_t i=0;i<n_l;i++){
        if(i<log2(n)){
            uint32_t j=1<<(i);
            M[i][0]= 1;
            M[i][j]= x[i];
        }else{
            M[i].resize(1);
            M[i][0]=1;
        }
    }

    FrVec u;
    FrVec v;
    int64_t index;
    index=0;
    uint8_t tmp=log2(n_l);
    for(uint8_t i=1;i<=tmp;i++){
        uint32_t L = (1<<(tmp-i));
        vector<FrVec> m(L, FrVec((n>>1)+1,0));
        for(uint64_t j=0;j<L;j++){
            u=M[index];
            index++;
            v=M[index];
            index++;
            m[j]=PolyMul(u,v);
        }
        index = 0;
        M=m;

    }

    G2* ptgg=&pi.gg[0];
    Fr* ptM=&M[0][0];

    
    G2 g_v;
    G2::mulVec(g_v,ptgg,ptM,pi.gg.size());

    GT res;
    pairing(res,pi.w,g_v);

    bool flag=(P_v==res);

    return flag;

}

struct zkIPPproof{
    
    IPPproof pi;
    GT Q;
    GT R;

    void Prove(vector<G2> gg, GT P, vector<G1> ww){
        
        vector<G1> vv(gg.size());
        for(uint32_t i=0;i<ww.size();i++){

            Fr c;
            c.setByCSPRNG();
            G1 v_i;
            hashAndMapToG1(v_i,c.getStr());
            vv[i]=v_i;

        }
        this->Q=MultiPairing(vv,gg);
        Fr chal;
        chal.setHashOf(Q.getStr());
        vector<G1> uu(ww.size());
        
        for(uint32_t i=0;i<ww.size();i++){
            uu[i]=(ww[i]*chal)+vv[i];
        }

        GT P_c;
        GT::pow(P_c,P,chal);
        this->R=P_c*(this->Q);
        
        IPPproof pi_zk(gg,P);
        pi_zk.Prove(gg,R,uu);
        this->pi=pi_zk;
        pi_zk.P=P;

    }

};

IPPproof zkIPPprove(vector<G2> gg, GT P, vector<G1> ww)
{
    vector<G1> vv(gg.size());
        for(uint32_t i=0;i<ww.size();i++){

            Fr c;
            c.setByCSPRNG();
            G1 v_i;
            hashAndMapToG1(v_i,c.getStr());
            vv[i]=v_i;

        }
        GT Q=MultiPairing(vv,gg);
        Fr chal;
        chal.setHashOf(Q.getStr());
        vector<G1> uu(ww.size());
        
        for(uint32_t i=0;i<ww.size();i++){
            uu[i]=(ww[i]*chal)+vv[i];
        }
        
        GT P_c;
        GT::pow(P_c,P,chal);
        GT R=P_c*Q;
        IPPproof pi_zk(gg,R);
        pi_zk.Q=Q;
        pi_zk.Prove(gg,R,uu);
        pi_zk.P=P;

        return pi_zk;
}

bool zkIPPverify(IPPproof pi){

    Fr chal;
    chal.setHashOf(pi.Q.getStr());
    GT P_c;
    GT::pow(P_c,pi.P,chal);
    GT R=P_c*pi.Q;
    pi.P=R;

    bool flag=IPPverify(pi);
    return flag;
}


struct ZKMP {

    zkbpacc_setup setup;
    G2 C_I;
    G2 C_A;

    ZKMP(zkbpacc_setup setup){
        this->setup=setup;
    }

    struct msg{

        G1 P_1;
        G1 P_2;
        G1 R_1;
        G1 R_2;
        GT R_3;

    }msg;

    struct response{

        Fr s_r_I;
        Fr s_r_A;
        Fr s_tau_1;
        Fr s_tau_2;
        Fr s_delta_1;
        Fr s_delta_2;

    }response;

    void prove(G2 C_I, G2 C_A, G1 pi_I, Fr r_I, Fr r_A){

        this->C_I=C_I;
        this->C_A=C_A;

        Fr tau_1;
        Fr tau_2;
        Fr r_r_I;
        Fr r_r_A;
        Fr r_tau_1;
        Fr r_tau_2;
        Fr r_delta_1;
        Fr r_delta_2;

        tau_1.setByCSPRNG();
        tau_2.setByCSPRNG();
        r_r_I.setByCSPRNG();
        r_r_A.setByCSPRNG();
        r_tau_1.setByCSPRNG();
        r_tau_2.setByCSPRNG();
        r_delta_1.setByCSPRNG();
        r_delta_2.setByCSPRNG();

        Fr delta_1 = r_I*tau_1;
        Fr delta_2 = r_I*tau_2;

        this->msg.P_1 = ((this->setup.g1si[0])*tau_1)+((this->setup.Frakg)*tau_2);
        this->msg.P_2 = pi_I+((this->setup.Frakg)*tau_1);
        this->msg.R_1 = ((this->setup.g1si[0])*r_tau_1)+((this->setup.Frakg)*r_tau_2);
        this->msg.R_2 = ((this->msg.P_1)*r_r_I)+((this->setup.g1si[0])*(-r_delta_1))+((this->setup.Frakg)*(-r_delta_2));
        
        GT R_3_1;
        GT R_3_2;
        GT R_3_3;
        GT R_3_4;

        pairing(R_3_1, this->setup.Frakg*r_tau_1, C_I);
        pairing(R_3_2, this->msg.P_2*r_r_I, this->setup.h2);
        pairing(R_3_3, this->setup.Frakg*r_delta_1, this->setup.h2);
        pairing(R_3_4, this->setup.g1si[0]*r_r_A, this->setup.h2);

        this->msg.R_3=R_3_1*R_3_2/R_3_3/R_3_4;

        string buf=this->msg.P_1.getStr()+this->msg.P_2.getStr()+this->msg.R_1.getStr()+this->msg.R_2.getStr()+this->msg.R_3.getStr();
        Fr c;
        c.setHashOf(buf);


        this->response.s_r_I=r_r_I+c*r_I;
        this->response.s_r_A=r_r_A+c*r_A;
        this->response.s_tau_1=r_tau_1+c*tau_1;
        this->response.s_tau_2=r_tau_2+c*tau_2;
        this->response.s_delta_1=r_delta_1+c*delta_1;
        this->response.s_delta_2=r_delta_2+c*delta_2;

    }

};

bool ZKMP_verify(ZKMP pi){

    string buf=pi.msg.P_1.getStr()+pi.msg.P_2.getStr()+pi.msg.R_1.getStr()+pi.msg.R_2.getStr()+pi.msg.R_3.getStr();
    Fr c;
    c.setHashOf(buf);


    G1 R_1p=pi.msg.P_1*(-c)+pi.setup.g1si[0]*pi.response.s_tau_1+pi.setup.Frakg*pi.response.s_tau_2;
    G1 R_2p=pi.msg.P_1*pi.response.s_r_I+pi.setup.g1si[0]*(-pi.response.s_delta_1)+pi.setup.Frakg*(-pi.response.s_delta_2);

    GT R_3_1p;
    GT R_3_2p;
    GT R_3_3p;
    GT R_3_4p;

    pairing(R_3_1p, pi.setup.Frakg*pi.response.s_tau_1, pi.C_I);
    pairing(R_3_2p, pi.msg.P_2*pi.response.s_r_I, pi.setup.h2);
    pairing(R_3_3p, pi.setup.Frakg*pi.response.s_delta_1, pi.setup.h2);
    pairing(R_3_4p, pi.setup.g1si[0]*pi.response.s_r_A, pi.setup.h2);

    GT R_3p=R_3_1p*R_3_2p/R_3_3p/R_3_4p;


    GT R_3_1o;
    GT R_3_2o;

    pairing(R_3_1o,pi.msg.P_2,pi.C_I*c);
    pairing(R_3_2o,pi.setup.g1si[0],pi.C_A*c);

    GT R_3o=pi.msg.R_3*(R_3_1o/R_3_2o);

    bool flag = false;
    flag = (pi.msg.R_1 == R_1p);
    flag = (pi.msg.R_2 == R_2p);
    flag = (R_3o==R_3p);

    return flag;
};

struct ZKNMP {
    zkbpacc_setup setup;
    G2 C_I;
    G2 C_A;

    ZKNMP(zkbpacc_setup setup){
        this->setup=setup;
    }

    struct msg{

        G1 P_1;
        G1 P_2;
        G1 Q_1;
        G1 Q_2;
        G1 R_1;
        G1 R_2;
        G1 R_3;
        G1 R_4;
        GT R_5;

    }msg;

    struct response{

        Fr s_r_I;
        Fr s_r_A;
        Fr s_tau_i[4];
        Fr s_delta_i[4];

    }response;

    void prove(G2 C_I, G2 C_A, G1 pi_1, G1 pi_2, Fr r_I, Fr r_A){
        
        Fr r_r_I;
        Fr r_r_A;
        Fr tau_i[4];
        Fr r_tau_i[4];
        Fr r_delta_i[4];
        Fr delta_i[4];
        r_r_I.setByCSPRNG();
        r_r_A.setByCSPRNG();

        for(int i=0;i<4;i++){

            tau_i[i].setByCSPRNG();
            r_tau_i[i].setByCSPRNG();
            r_delta_i[i].setByCSPRNG();

        }

        delta_i[0]=r_I*tau_i[0];
        delta_i[1]=r_I*tau_i[1];
        delta_i[2]=r_A*tau_i[2];
        delta_i[3]=r_A*tau_i[3];

        this->msg.P_1 = this->setup.g1si[0]*tau_i[0]+this->setup.Frakg*tau_i[1];
        this->msg.P_2 = pi_1+this->setup.Frakg*tau_i[0];
        this->msg.Q_1 = this->setup.g1si[0]*tau_i[2]+this->setup.Frakh*tau_i[3];
        this->msg.Q_2 = pi_2+this->setup.Frakh*tau_i[2];
        this->msg.R_1 = this->setup.g1si[0]*r_tau_i[0]+this->setup.Frakg*r_tau_i[1];
        this->msg.R_2 = this->msg.P_1*r_r_I-this->setup.g1si[0]*r_delta_i[0]-this->setup.Frakg*r_delta_i[1];
        this->msg.R_3 = this->setup.g1si[0]*r_tau_i[2]+this->setup.Frakh*r_tau_i[3];
        this->msg.R_4 = this->msg.Q_1*r_r_A-this->setup.g1si[0]*r_delta_i[2]-this->setup.Frakh*r_delta_i[3];

        GT R_5_1;
        GT R_5_2;
        GT R_5_3;
        GT R_5_4;
        GT R_5_5;
        GT R_5_6;

        pairing(R_5_1, this->setup.Frakg*r_tau_i[0], C_I);
        pairing(R_5_2, this->setup.Frakh*r_tau_i[2], C_A);
        pairing(R_5_3, pi_1+this->setup.Frakg*tau_i[0]*r_r_I, this->setup.h2);
        pairing(R_5_4, pi_2+this->setup.Frakh*tau_i[2]*r_r_A, this->setup.h2);
        pairing(R_5_5, this->setup.Frakg*r_delta_i[0], this->setup.h2);
        pairing(R_5_6, this->setup.Frakh*r_delta_i[2], this->setup.h2);

        this->msg.R_5 = (R_5_1*R_5_2*R_5_3*R_5_4)/(R_5_5*R_5_6);

        string buf=this->msg.P_1.getStr()+this->msg.P_2.getStr()+this->msg.Q_1.getStr()+this->msg.Q_2.getStr()+this->msg.R_1.getStr()+this->msg.R_2.getStr()
        +this->msg.R_3.getStr()+this->msg.R_4.getStr()+this->msg.R_5.getStr();
        Fr c;
        c.setHashOf(buf); 

        this->response.s_r_I=r_r_I+c*r_I;
        this->response.s_r_A=r_r_A+c*r_A;
        for(int i=0;i<4;i++){

            this->response.s_tau_i[i]=r_tau_i[i]+c*tau_i[i];
            this->response.s_delta_i[i]=r_delta_i[i]+c*delta_i[i];

        }

    }
};

bool ZKNMP_verify(ZKNMP pi){
    string buf=pi.msg.P_1.getStr()+pi.msg.P_2.getStr()+pi.msg.Q_1.getStr()+pi.msg.Q_2.getStr()+pi.msg.R_1.getStr()+pi.msg.R_2.getStr()
    +pi.msg.R_3.getStr()+pi.msg.R_4.getStr()+pi.msg.R_5.getStr();
    Fr c;
    c.setHashOf(buf);

    G1 R_1v = pi.msg.P_1*(-c) + pi.setup.g1si[0]*pi.response.s_tau_i[0] + pi.setup.Frakg*pi.response.s_tau_i[1];
    G1 R_2v = pi.msg.P_1*pi.response.s_r_I + pi.setup.g1si[0]*(-pi.response.s_delta_i[0]) + pi.setup.Frakg*(-pi.response.s_delta_i[1]);
    G1 R_3v = pi.msg.Q_1*(-c) + pi.setup.g1si[0]*pi.response.s_tau_i[2] + pi.setup.Frakh*pi.response.s_tau_i[3];
    G1 R_4v = pi.msg.Q_1*pi.response.s_r_A + pi.setup.g1si[0]*(-pi.response.s_delta_i[2]) + pi.setup.Frakh*(-pi.response.s_delta_i[3]);

    GT R_5_1p;
    GT R_5_2p;
    GT R_5_3p;
    GT R_5_4p;
    GT R_5_5p;
    GT R_5_6p;

    pairing(R_5_1p, pi.setup.Frakg*pi.response.s_tau_i[0], pi.C_I);
    pairing(R_5_2p, pi.setup.Frakh*pi.response.s_tau_i[2], pi.C_A);
    pairing(R_5_3p, pi.msg.P_2*pi.response.s_r_I, pi.setup.h2);
    pairing(R_5_4p, pi.msg.Q_2*pi.response.s_r_A, pi.setup.h2);
    pairing(R_5_5p, pi.setup.Frakg*pi.response.s_delta_i[0], pi.setup.h2);
    pairing(R_5_6p, pi.setup.Frakh*pi.response.s_delta_i[2], pi.setup.h2);

    GT R_5p = (R_5_1p*R_5_2p*R_5_3p*R_5_4p)/(R_5_5p*R_5_6p);

    GT R_5_1o;
    GT R_5_2o;
    GT R_5_3o;

    pairing(R_5_1o, pi.msg.P_2, pi.C_I*c);
    pairing(R_5_2o, pi.msg.Q_2, pi.C_A*c);
    pairing(R_5_3o, pi.setup.g1si[0], pi.setup.g2si[0]*c);

    GT R_5o = (R_5_1o*R_5_2o)/R_5_3o;
    R_5o = pi.msg.R_5*R_5o;

    bool flag=false;
    flag = (pi.msg.R_1 == R_1v);
    flag = (pi.msg.R_2 == R_2v);
    flag = (pi.msg.R_3 == R_3v);
    flag = (pi.msg.R_4 == R_4v);
    
    return flag;
}

struct ZKSP {

    zkbpacc_setup setup;
    G2 C_I;
    G2 C_J;
    G2 C_A;

    ZKSP(zkbpacc_setup setup){
        this->setup=setup;
    }

    struct msg{

        G1 P_1;
        G1 P_2;
        G1 P_3;
        G1 P_4;
        G1 R_1;
        G1 R_2;
        G1 R_3;
        G1 R_4;
        GT R_5;
        GT R_6;
        GT R_7;

    }msg;

    struct response{

        Fr s_r_I;
        Fr s_r_J;
        Fr s_r_A;
        Fr s_tau_i[4];
        Fr s_delta_i[4];

    }response;

    void prove(G2 C_I, G2 C_J, G2 C_A, G1 pi_I, G1 pi_J, Fr r_I, Fr r_J, Fr r_A){

        this->C_I=C_I;
        this->C_J=C_J;
        this->C_A=C_A;

        Fr r_r_I;
        Fr r_r_J;
        Fr r_r_A;
        Fr tau_i[4];
        Fr r_tau_i[4];
        Fr r_delta_i[4];
        Fr delta_i[4];
        r_r_I.setByCSPRNG();
        r_r_A.setByCSPRNG();
        r_r_J.setByCSPRNG();

        for(int i=0;i<4;i++){

            tau_i[i].setByCSPRNG();
            r_tau_i[i].setByCSPRNG();
            r_delta_i[i].setByCSPRNG();

        }

        delta_i[0]=r_I*tau_i[0];
        delta_i[1]=r_I*tau_i[1];
        delta_i[2]=r_J*tau_i[2];
        delta_i[3]=r_J*tau_i[3];

        this->msg.P_1 = ((this->setup.g1si[0])*tau_i[0])+((this->setup.Frakg)*tau_i[1]);
        this->msg.P_2 = pi_I+((this->setup.Frakg)*tau_i[0]);
        this->msg.P_3 = ((this->setup.g1si[0])*tau_i[2])+((this->setup.Frakg)*tau_i[3]);
        this->msg.P_4 = pi_J+((this->setup.Frakg)*tau_i[2]);

        this->msg.R_1 = ((this->setup.g1si[0])*r_tau_i[0])+((this->setup.Frakg)*r_tau_i[1]);
        this->msg.R_2 = ((this->msg.P_1)*r_r_I)+((this->setup.g1si[0])*(-r_delta_i[0]))+((this->setup.Frakg)*(-r_delta_i[1]));
        this->msg.R_3 = ((this->setup.g1si[0])*r_tau_i[2])+((this->setup.Frakg)*r_tau_i[3]);
        this->msg.R_4 = ((this->msg.P_3)*r_r_J)+((this->setup.g1si[0])*(-r_delta_i[2]))+((this->setup.Frakg)*(-r_delta_i[3]));


        GT R_5_1;
        GT R_5_2;
        GT R_5_3;
        GT R_5_4;
        GT pre1; pairing(pre1,this->setup.Frakg,this->setup.h2);
        GT pre2; pairing(pre2,this->setup.g1si[0],this->setup.h2);

        //

        pairing(R_5_1, this->setup.Frakg*r_tau_i[0], C_I);

        //
        
        pairing(R_5_2, this->msg.P_2*r_r_I, this->setup.h2);
        //pairing(R_5_3, this->setup.Frakg, this->setup.h2*r_delta_i[0]);
        GT::pow(R_5_3, pre1, r_delta_i[0]);
        //pairing(R_5_4, this->setup.g1si[0]*r_r_A, this->setup.h2);
        GT::pow(R_5_4, pre2, r_r_A);

        this->msg.R_5=(R_5_1*R_5_2)/(R_5_3*R_5_4);


        GT R_6_1;
        GT R_6_2;
        GT R_6_3;
        GT R_6_4;

        pairing(R_6_1, this->setup.Frakg*r_tau_i[2], C_J);
        pairing(R_6_2, this->msg.P_4*r_r_J, this->setup.h2);
        //pairing(R_6_3, this->setup.Frakg, this->setup.h2*r_delta_i[2]);
        GT::pow(R_6_3, pre1, r_delta_i[2]);
        //pairing(R_6_4, this->setup.g1si[0], this->setup.h2*r_r_A);
        R_6_4 = R_5_4;

        this->msg.R_6=(R_6_1*R_6_2)/(R_6_3*R_6_4);

        GT R_7_1;
        GT R_7_2;

        //pairing(R_7_1, this->setup.g1si[0], this->setup.h2*r_r_I);
        GT::pow(R_7_1, pre2, r_r_I);
        pairing(R_7_2, this->setup.Frakg*(-r_tau_i[2]), this->setup.g2si[0]);

        this->msg.R_7 = R_7_1*R_7_2;

        string buf=this->msg.P_1.getStr()+this->msg.P_2.getStr()+this->msg.P_3.getStr()+this->msg.P_4.getStr()+this->msg.R_1.getStr()+this->msg.R_2.getStr()
        +this->msg.R_3.getStr()+this->msg.R_4.getStr()+this->msg.R_5.getStr()+this->msg.R_6.getStr()+this->msg.R_7.getStr();
        Fr c;
        c.setHashOf(buf);


        this->response.s_r_I=r_r_I+c*r_I;
        this->response.s_r_J=r_r_J+c*r_J;
        this->response.s_r_A=r_r_A+c*r_A;
        for(int i=0;i<4;i++){

            this->response.s_tau_i[i]=r_tau_i[i]+c*tau_i[i];
            this->response.s_delta_i[i]=r_delta_i[i]+c*delta_i[i];

        }

    }

};

bool ZKSP_verify(ZKSP pi){

    string buf=pi.msg.P_1.getStr()+pi.msg.P_2.getStr()+pi.msg.P_3.getStr()+pi.msg.P_4.getStr()+pi.msg.R_1.getStr()+pi.msg.R_2.getStr()
        +pi.msg.R_3.getStr()+pi.msg.R_4.getStr()+pi.msg.R_5.getStr()+pi.msg.R_6.getStr()+pi.msg.R_7.getStr();
    Fr c;
    c.setHashOf(buf);
    
    G1 R_1v=pi.msg.P_1*(-c)+pi.setup.g1si[0]*pi.response.s_tau_i[0]+pi.setup.Frakg*pi.response.s_tau_i[1];
    G1 R_2v=pi.msg.P_1*pi.response.s_r_I+pi.setup.g1si[0]*(-pi.response.s_delta_i[0])+pi.setup.Frakg*(-pi.response.s_delta_i[1]);
    G1 R_3v=pi.msg.P_3*(-c)+pi.setup.g1si[0]*pi.response.s_tau_i[2]+pi.setup.Frakg*pi.response.s_tau_i[3];
    G1 R_4v=pi.msg.P_3*pi.response.s_r_J+pi.setup.g1si[0]*(-pi.response.s_delta_i[2])+pi.setup.Frakg*(-pi.response.s_delta_i[3]);

    GT R_5_1p;
    GT R_5_2p;
    GT R_5_3p;
    GT R_5_4p;
    GT pre1;pairing(pre1,pi.setup.Frakg,pi.setup.h2);
    //GT pre2;pairing(pre2,pi.setup.g1si[0],pi.setup.h2);

    pairing(R_5_1p, pi.setup.Frakg*pi.response.s_tau_i[0], pi.C_I);
    pairing(R_5_2p, pi.msg.P_2*pi.response.s_r_I, pi.setup.h2);
    //pairing(R_5_3p, pi.setup.Frakg, pi.setup.h2*pi.response.s_delta_i[0]);
    GT::pow(R_5_3p, pre1, pi.response.s_delta_i[0]);
    pairing(R_5_4p, pi.setup.g1si[0]*pi.response.s_r_A, pi.setup.h2);
    //GT::pow(R_5_4p, pre2, pi.response.s_r_A);

    GT R_5p=R_5_1p*R_5_2p/R_5_3p/R_5_4p;


    GT R_5_1o;
    GT R_5_2o;

    pairing(R_5_1o,pi.msg.P_2*c,pi.C_I);
    pairing(R_5_2o,pi.setup.g1si[0]*c,pi.C_A);

    GT R_5o=pi.msg.R_5*(R_5_1o/R_5_2o);

    GT R_6_1p;
    GT R_6_2p;
    GT R_6_3p;
    GT R_6_4p;

    pairing(R_6_1p, pi.setup.Frakg*pi.response.s_tau_i[2], pi.C_J);
    pairing(R_6_2p, pi.msg.P_4*pi.response.s_r_J, pi.setup.h2);
    //pairing(R_6_3p, pi.setup.Frakg, pi.setup.h2*pi.response.s_delta_i[2]);
    GT::pow(R_6_3p, pre1, pi.response.s_delta_i[2]);
    //pairing(R_6_4p, pi.setup.g1si[0], pi.setup.h2*pi.response.s_r_A);
    R_6_4p = R_5_4p;

    GT R_6p=R_6_1p*R_6_2p/R_6_3p/R_6_4p;


    GT R_6_1o;
    GT R_6_2o;

    pairing(R_6_1o,pi.msg.P_4*c,pi.C_J);
    pairing(R_6_2o,pi.setup.g1si[0]*c,pi.C_A);

    GT R_6o=pi.msg.R_6*(R_6_1o/R_6_2o);

    GT R_7_1p;
    GT R_7_2p;

    pairing(R_7_1p, pi.setup.g1si[0]*pi.response.s_r_I, pi.setup.h2);
    pairing(R_7_2p, pi.setup.Frakg*(-pi.response.s_tau_i[2]), pi.setup.g2si[0]);

    GT R_7p = R_7_1p*R_7_2p;

    GT R_7_1o;
    GT R_7_2o;

    pairing(R_7_1o,pi.setup.g1si[0]*c,pi.C_I);
    pairing(R_7_2o,pi.msg.P_4*c,pi.setup.g2si[0]);

    GT R_7o=pi.msg.R_7*(R_7_1o/R_7_2o);



    bool flag=false;

    flag = (pi.msg.R_1 == R_1v);
    flag = (pi.msg.R_2 == R_2v);
    flag = (pi.msg.R_3 == R_3v);
    flag = (pi.msg.R_4 == R_4v);
    flag = (R_5p == R_5o);
    flag = (R_6p == R_6o);
    flag = (R_7p == R_7o);

    return flag;

}

struct PoK_proof_g2{

    G2 P;
    G2 R;
    Fr s;

    void prove(G2 P, Fr x){

        this->P = P;
        Fr tau;
        tau.setByCSPRNG();
        G2 g2;
        hashAndMapToG2(g2,"g2");
        R = g2*tau;
        Fr c;
        c.setHashOf(R.getStr());
        this->s = tau + c*x;


    }

};

struct PoK_proof_g2_n{

    G2 P;
    G2 R;
    FrVec ss;

    void prove(vector<G2> gg, G2 P, FrVec xx, uint32_t n){

        assert(gg.size()==n);

        FrVec tau; tau.resize(n);
        for(uint32_t i=0;i<n;i++){
            tau[i].setByCSPRNG();
        }

        G2::mulVec(this->R, &gg[0], &tau[0], n);
        Fr c; c.setHashOf(R.getStr());
        this->ss = PolyAdd(tau, PolyMul({c},xx));

    }

};

// int main(){
    
//     initPairing(mcl::BLS12_381);


//     zkbpacc_setup setup1;
//     setup1.init(8);

//     FrVec poly1 = {6, 11, 6, 1};
//     FrVec poly2 = {1, 4};

//     FrvT_3 out=xGCD(poly1, poly2);

//     FrVec alpha = get<0>(out);
//     FrVec beta = get<1>(out);
//     FrVec gcd = get<2>(out);

//     PCS pcs1;
//     pcs1.setup(8);

//     G1 pi_1 = pcs1.commit_G1(alpha);
//     G1 pi_2 = pcs1.commit_G1(beta);

//     G2 C_A = pcs1.commit(poly1);
//     G2 C_I = pcs1.commit(poly2);
//     Fr r_I;
//     r_I.setByCSPRNG();
//     Fr r_A;
//     r_A.setByCSPRNG();
//     C_I = C_I+setup1.h2*r_I;
//     C_A = C_A+setup1.h2*r_A;

//     ZKNMP nmp(setup1);
//     nmp.prove(C_I, C_A, pi_1, pi_2, r_I, r_A);

//     bool flag = ZKNMP_verify(nmp);

//     std::cout<<"===================================="<<endl;

//     zkbpacc_setup setup2;
//     setup2.init(8);

//     FrVec Ax = {24, -50, 35, -10, 1};
//     FrVec Ix = {2, -3, 1};
//     FrVec Jx = {12, -7, 1};

//     PCS pcs2;
//     pcs2.setup(8);

//     C_A = pcs2.commit(Ax);
//     C_I = pcs2.commit(Ix);
//     G2 C_J = pcs2.commit(Jx);

//     G1 pi_I = pcs2.commit_G1(PolyLongDiv(Ax,Ix));
//     G1 pi_J = pcs2.commit_G1(PolyLongDiv(Ax,Jx));

//     r_A;
//     r_A.setByCSPRNG();
//     r_I;
//     r_I.setByCSPRNG();
//     Fr r_J;
//     r_J.setByCSPRNG();

//     C_A = C_A + setup2.h2*r_A;
//     C_I = C_I + setup2.h2*r_I;
//     C_J = C_J + setup2.h2*r_J;

//     ZKSP pi3(setup2);
//     pi3.prove(C_I, C_J, C_A, pi_I, pi_J, r_I, r_J, r_A);
//     ZKSP_verify(pi3);

//     cout<<"===================================="<<endl;

//     ZKMP pi4(setup2);
//     pi4.prove(C_I, C_A, pi_I, r_I, r_A);
//     bool fla = ZKMP_verify(pi4);
//     cout<<fla<<endl;

// }

