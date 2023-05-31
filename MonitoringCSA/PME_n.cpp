#include <mcl/bls12_381.hpp>
#include <iostream>
#include <limits>
#include <cmath>
#include <vector>
#include <malloc.h>
#include <mcl/ntt.hpp>
#include <chrono>
#include <fstream>
#include <map>
#include <cstdlib>
#include <ctime>
#include "primitiveroots.hpp"
#include "polyonfr.hpp"
#include "mainalgorithms.hpp"

using namespace mcl;
using namespace mcl::bn;
using namespace std;
using namespace chrono;

struct PME_PublicInput {
    G1 g1;
    G2 g2;
    G2 A;
    G1 H_1;
    G1 H_2;
    G2 G_1;
    G2 G_2;
};

struct PME_response{

    G1 P_i;
    G1 pi;

    void p(G1 C, Fr d, Fr sigma, G1 H, Fr e, int i){
        this->P_i = C*e*d + H*(pow(1,i))*d;
        this->pi = (C*e*(1/sigma))+(H*(pow(1,i))*(1/sigma));
    }

};

Fr PME_v1(Fr a, Fr x){

    return a*x;

}

G1 PME_v2(Fr a, Fr d1, Fr d2, PME_response P_1, PME_response P_2, Fr e1, Fr e2, PME_PublicInput Input, zkbpacc_setup_pcs setup){

    G1 P = P_1.P_i*(1/d1)+P_2.P_i*(1/d2);

    GT lhs;
    GT rhs;

    GT lhs_1;
    GT lhs_2;

    pairing(lhs_1, P_1.pi+Input.H_2, Input.G_1);
    pairing(lhs_2, P_2.pi+Input.H_1*(-1), Input.G_2);
    lhs = lhs_1*lhs_2;

    pairing(rhs, P, setup.h2);

    GT lhs2;
    GT rhs2;

    FrVec tmp = {d1*d2, d1+d2, 1};
    G2 L;
    L.clear();
    for(int i=0; i<tmp.size(); i++){
        L = L + setup.pcs.pp.g2si[i]*tmp[i];
    }
    FrVec tmp2 = {e1*d2+e2*d1, e1+e2};
    G1 R;
    R.clear();
    for(int i=0; i<tmp2.size(); i++){
        R = R + setup.pcs.pp.g1si[i]*tmp[i];
    }

    pairing(lhs2, P, L);
    pairing(rhs2, R, Input.A);

    if((lhs==rhs)&&(lhs2==rhs2)){
        return P*(1/a);
    }else{
        cout<<"error"<<endl;
    }

}

// int main(){

//     initPairing(mcl::BLS12_381);

//     PCS pcs;
//     pcs.setup(1000);

//     zkbpacc_setup_pcs setup2(pcs);

//     IDset idset;
//     idset.init(100);

//     PME_PublicInput Input;
//     Input.g1 = setup2.pcs.pp.g1si[0];
//     Input.g2 = setup2.pcs.pp.g2si[0];
    

//     G2 A = pcs.commit(idset.idset, 100);
//     FrVec A_p = Polytree(arr2Vec(idset.idset,100));
//     Input.A = A;

//     FrVec I_1 = {idset.idset[3],1};
//     G1 C_1 = pcs.commit_G1(PolyLongDiv(A_p,I_1));
//     FrVec I_2 = {idset.idset[43],1};
//     G1 C_2 = pcs.commit_G1(PolyLongDiv(A_p,I_2));

//     Fr sigma1;sigma1.setByCSPRNG();
//     Fr sigma2;sigma2.setByCSPRNG();

//     Input.H_1 = setup2.Frakh*(sigma1);
//     Input.H_2 = setup2.Frakh*(sigma2);
//     Input.G_1 = setup2.pcs.pp.g2si[0]*(sigma1);
//     Input.G_2 = setup2.pcs.pp.g2si[0]*(sigma2);

//     Fr a;
//     a.setByCSPRNG();
//     Fr x1; x1.setByCSPRNG();
//     Fr x2; x2.setByCSPRNG();
//     Fr e1 = PME_v1(a, x1);
//     Fr e2 = PME_v1(a, x2);

//     PME_response P_1;
//     PME_response P_2;


//     G1 H = setup2.Frakh*(sigma1*sigma2);
//     P_1.p(C_1, idset.idset[3], sigma1, H, e1, 1);
//     P_1.p(C_2, idset.idset[43], sigma2, H, e2, 2);

//     PME_v2(a, idset.idset[3], idset.idset[43], P_1, P_2, e1, e2, Input, setup2);

// }