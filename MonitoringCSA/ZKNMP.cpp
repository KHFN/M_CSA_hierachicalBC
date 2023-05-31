#include <mcl/bls12_381.hpp>
#include <iostream>
#include <limits>
#include <cmath>
#include <vector>
#include <malloc.h>
#include <mcl/ntt.hpp>
#include <chrono>
#include <fstream>
#include "primitiveroots.hpp"
#include "polyonfr.hpp"
#include "mainalgorithms.hpp"


using namespace mcl;
using namespace mcl::bn;
using namespace std;
using namespace chrono;

struct ZKNMP_n {

    zkbpacc_setup_pcs setup;

    G1 P_1, P_2, Q_1, Q_2, R_1, R_2, R_3, R_4, R_5;
    Fr s_r_I=0, s_r_A=0;
    FrVec s_tau, s_delta;
    

    ZKNMP_n(zkbpacc_setup_pcs set): setup(set) {}

    void prove(G2 C_I, G2 C_A, G1 w_1, G1 w_2, Fr r_I, Fr r_A){

        Fr r_r_I, r_r_A;
        r_r_I.setByCSPRNG(); r_r_A.setByCSPRNG();
        FrVec tau(4), s_tau(4), s_delta(4);

        for(int i=0;i<4;i++){

            tau[i].setByCSPRNG();
            s_tau[i].setByCSPRNG();
            s_delta[i].setByCSPRNG();

        }        

    }

};

int main(){

    initPairing(mcl::BLS12_381);

    PCS pcs;
    pcs.setup(100);

    FrVec p = {1,2,3,4,5,6,7};
    FrVec q = {8,9};

    

    FrVec A = Polytree(p);
    FrVec I = Polytree(q);
    
    G2 W_A = pcs.commit(A);
    G2 W_I = pcs.commit(I);

}