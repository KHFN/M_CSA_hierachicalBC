#include <mcl/bls12_381.hpp>
#include <iostream>
#include <limits>
#include <cmath>
#include <vector>
#include <malloc.h>
#include <mcl/ntt.hpp>
#include <chrono>
#include <fstream>
#include <iomanip> 
#include "primitiveroots.hpp"
#include "polyonfr.hpp"
#include "mainalgorithms.hpp"


using namespace mcl;
using namespace mcl::bn;
using namespace std;
using namespace chrono;

void ZKMP_test(){

    ofstream fout1;
    ofstream fout2;
    fout1.open("ZKMP_p.txt");
    fout2.open("ZKMP_v.txt");

    zkbpacc_setup setup;
    setup.init(8);

    PCS pcs;
    pcs.setup(8);

    FrVec A = {120, 274, 225, 85, 15, 1};
    FrVec I = {1,1};

    G2 C_A = pcs.commit(A);
    G2 C_I = pcs.commit(I);
    Fr r_I;
    r_I.setByCSPRNG();
    Fr r_A;
    r_A.setByCSPRNG();
    C_I = C_I+setup.h2*r_I;
    C_A = C_A+setup.h2*r_A;

    G1 pi_I = pcs.commit_G1(PolyLongDiv(A,I));

    ZKMP zkmp(setup);

    microseconds total_p(0);
    microseconds total_v(0);
    microseconds elaps;

    for(uint32_t i=0;i<10000;i++){

        auto start_1 = chrono::steady_clock::now();

        zkmp.prove(C_I, C_A, pi_I, r_I, r_A);

        auto end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        fout1<< elaps.count()<< ", ";
        total_p = total_p+elaps;

        auto start_2 = chrono::steady_clock::now();

        bool flag = ZKMP_verify(zkmp);

        auto end_2 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_2-start_2);
        fout2<< elaps.count()<< ", ";
        total_v = total_v+elaps;

    }

    std::cout << "ZKMP_prove = : " << total_p.count()/10000 << "(ms)"<<std::endl;
    std::cout << "ZKMP_verifiy = : " << total_v.count()/10000 << "(ms)"<<std::endl;

}

void ZKNMP_test(){

    ofstream fout1;
    ofstream fout2;
    fout1.open("ZKNMP_p.txt");
    fout2.open("ZKNMP_v.txt");

    zkbpacc_setup setup1;
    setup1.init(8);

    FrVec poly1 = {6, 11, 6, 1};
    FrVec poly2 = {1, 4};

    FrvT_3 out=xGCD(poly1, poly2);

    FrVec alpha = get<0>(out);
    FrVec beta = get<1>(out);
    FrVec gcd = get<2>(out);

    PCS pcs1;
    pcs1.setup(8);

    G1 pi_1 = pcs1.commit_G1(alpha);
    G1 pi_2 = pcs1.commit_G1(beta);

    G2 C_A = pcs1.commit(poly1);
    G2 C_I = pcs1.commit(poly2);
    Fr r_I;
    r_I.setByCSPRNG();
    Fr r_A;
    r_A.setByCSPRNG();
    C_I = C_I+setup1.h2*r_I;
    C_A = C_A+setup1.h2*r_A;

    ZKNMP nmp(setup1);

    microseconds total_p(0);
    microseconds total_v(0);
    microseconds elaps;

    for(uint32_t i=0;i<10000;i++){

        auto start_1 = chrono::steady_clock::now();

        nmp.prove(C_I, C_A, pi_1, pi_2, r_I, r_A);

        auto end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        fout1<< elaps.count()<< ", ";
        total_p = total_p+elaps;

        auto start_2 = chrono::steady_clock::now();

        bool flag = ZKNMP_verify(nmp);

        auto end_2 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_2-start_2);
        fout2<< elaps.count()<< ", ";
        total_v = total_v+elaps;

    }

    std::cout << "ZKNMP_prove = : " << total_p.count()/10000 << "(ms)"<<std::endl;
    std::cout << "ZKNMP_verifiy = : " << total_v.count()/10000 << "(ms)"<<std::endl;

}



void ZKSP_test(){

    ofstream fout1;
    ofstream fout2;
    fout1.open("ZKSP_p.txt");
    fout2.open("ZKSP_v.txt");

    zkbpacc_setup setup2;
    setup2.init(8);

    FrVec Ax = {24, -50, 35, -10, 1};
    FrVec Ix = {2, -3, 1};
    FrVec Jx = {12, -7, 1};

    PCS pcs2;
    pcs2.setup(8);

    G2 C_A = pcs2.commit(Ax);
    G2 C_I = pcs2.commit(Ix);
    G2 C_J = pcs2.commit(Jx);

    G1 pi_I = pcs2.commit_G1(PolyLongDiv(Ax,Ix));
    G1 pi_J = pcs2.commit_G1(PolyLongDiv(Ax,Jx));

    Fr r_A;
    r_A.setByCSPRNG();
    Fr r_I;
    r_I.setByCSPRNG();
    Fr r_J;
    r_J.setByCSPRNG();

    C_A = C_A + setup2.h2*r_A;
    C_I = C_I + setup2.h2*r_I;
    C_J = C_J + setup2.h2*r_J;

    ZKSP pi3(setup2);

    microseconds total_p(0);
    microseconds total_v(0);
    microseconds elaps;

    for(uint32_t i=0;i<10000;i++){

        auto start_1 = chrono::steady_clock::now();

        pi3.prove(C_I, C_J, C_A, pi_I, pi_J, r_I, r_J, r_A);

        auto end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        fout1<< elaps.count()<< ", ";
        total_p = total_p+elaps;

        auto start_2 = chrono::steady_clock::now();

        bool flag = ZKSP_verify(pi3);

        auto end_2 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_2-start_2);
        fout2<< elaps.count()<< ", ";
        total_v = total_v+elaps;

    }

    std::cout << "ZKSP_prove = : " << total_p.count()/10000 << "(ms)"<<std::endl;
    std::cout << "ZKSP_verifiy = : " << total_v.count()/10000 << "(ms)"<<std::endl;

}

void zkIPP_test(const std::string& file_name){

    vector<G2> gg;
    vector<G1> ww;

    std::ofstream fout1(file_name, std::ios::app);

    fout1<<"input size, "<<"Prove time"<<"verification time, "<<endl;

    //int values[] = {127, 255, 511, 1023, 2047, 4095, 8191, 16383};

    for(int i=128; i<=16384; i=i*2){

        
        i=i-1;
        microseconds elaps;

        i=nextPowOf2(i);
        fout1<< i << ", ";

        gg.resize(i);
        ww.resize(i);

        for(uint32_t i=0;i<ww.size();i++){

            hashAndMapToG1(ww[i],to_string(i));
            hashAndMapToG2(gg[i],to_string(i));

        }

        GT P =MultiPairing(ww,gg);

        auto start_1 = chrono::steady_clock::now();

        IPPproof ipp_pi = zkIPPprove(gg,P,ww);

        auto end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        fout1<< elaps.count()<< ", ";
        

        auto start_2 = chrono::steady_clock::now();

        bool flag = zkIPPverify(ipp_pi);

        auto end_2 = chrono::steady_clock::now();

        elaps = duration_cast<microseconds>(end_2-start_2);
        fout1<< elaps.count()<< ", ";

        fout1<<endl;
        


    }


}

void setuptest(const std::string& file_name){

    std::ofstream fout1(file_name, std::ios::app);
    microseconds elaps;

        

        for(uint32_t i=2500;i<320001;i=i*2){
            
            for(uint32_t j=0;j<100;j++){

                fout1<<i<<", ";
                

                FrVec IDset;
                IDset.resize(i);
                for(uint32_t k=0;k<i;k++){
                    IDset[k].setByCSPRNG();
                }

                auto start_1 = chrono::steady_clock::now();

                PCS pcs;
                pcs.setup(i+1);
                zkbpacc_setup setup;
                setup.init(i+1);
                setup.g1si = &pcs.pp.g1si[0];
                setup.g2si = &pcs.pp.g2si[0];

                auto end_1 = chrono::steady_clock::now();
                elaps = duration_cast<microseconds>(end_1-start_1);
                fout1<< elaps.count()<< ", ";

                Fr s = pcs.s;

                start_1 = chrono::steady_clock::now();

                // FrVec Ax = Polytree(IDset);
                // G2 C_A = pcs.commit(Ax);

                Fr Ax=1;
                G2 C_A;

                for(uint32_t k=0;k<IDset.size();k++){

                    Ax = Ax*(s+IDset[k]);

                    float progress = (float)k / (IDset.size() - 1);
                    int barWidth = 70;
                    
                    std::cout << "[";
                    int pos = barWidth * progress;
                    for (int i = 0; i < barWidth; ++i) {
                        if (i < pos) std::cout << "=";
                        else if (i == pos) std::cout << ">";
                        else std::cout << " ";
                    }
                    std::cout << "] " << int(progress * 100.0) << " %\r";
                    std::cout.flush();

                }

                C_A = setup.g2si[0]*Ax;

                end_1 = chrono::steady_clock::now();
                elaps = duration_cast<microseconds>(end_1-start_1);
                fout1<< elaps.count()<< ", ";

                start_1 = chrono::steady_clock::now();

                Fr id = IDset[rand()%i];
                Fr Ix = s+id;
                G2 C_id = setup.g2si[0]*(s+id);

                end_1 = chrono::steady_clock::now();
                elaps = duration_cast<microseconds>(end_1-start_1);
                fout1<< elaps.count()<< ", ";
                

                start_1 = chrono::steady_clock::now();

                G1 wit = setup.g1si[0]*(Ax/Ix);

                end_1 = chrono::steady_clock::now();
                elaps = duration_cast<microseconds>(end_1-start_1);
                fout1<< elaps.count()<< ", ";
                fout1<<endl;

            }

        }
    

}

int main(){

    initPairing(mcl::BLS12_381);

    ZKMP_test();
    ZKSP_test();
    ZKNMP_test();
    zkIPP_test("file_name.csv");
    setuptest("setup_new.csv");


    return 0;
}