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

struct PME_public {

    G1 g1;
    G2 g2;
    G2 A;
    G1 H_1;
    G1 H_2;

};

Fr PME_V1(Fr x, Fr alpha){

    Fr e = x*alpha;
    return e;

}

G1 PME_P1(Fr e, G1 w, G1 H, Fr d, int i=1){

    G1 P;

    P = w*(e*d)+H*(pow(-1,i)*d);

    return P;

}

G1 PME_V2(G1 P_1, G1 P_2, G2 A, Fr d1, Fr d2, Fr e1, Fr e2, Fr alpha, PCS pcs){

    G1 P = P_1*(1/d1)+P_2*(1/d2);
    GT L;
    GT R;

    FrVec tmp1 = {d1*d2, d1+d2, 1};
    FrVec tmp2 = PolyAdd(PolyMul({e1},Polytree({d2})),PolyMul({e2},Polytree({d1})));

   G2 L_g2 = pcs.commit(tmp1);
   G1 R_g1 = pcs.commit_G1(tmp2);

   pairing(L, P, L_g2);
   pairing(R, R_g1, A);

   if(L==R){
        return P*(1/alpha);
   }else{
        cout<<"error"<<endl;
   }

}


