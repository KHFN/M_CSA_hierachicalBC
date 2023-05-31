#include <mcl/bls12_381.hpp>
#include <mcl/bn_c384.h>
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
#include "PME.hpp"

struct Agg_input {

    PCS pcs;
    vector<G2> C_i;
    FrVec ri;
    FrVec id_i;
    G2 C_A;
    
};

struct Agg_output{

    G2 C_I;
    G2 C_Ip;
    ZKMP zkmp;
    IPPproof PI_ipp;
    G2 pi_open;

};

struct Agg_PME_prover_input{

    G1 H;
    vector<G1> w;

};

Agg_output AggZKMP(Agg_input input){

    assert((input.C_i.size()%2==0));
    assert((input.id_i.size()%2==0));
    assert((input.ri.size()%2==0));
    assert(input.C_i.size()==input.id_i.size());
    assert(input.C_i.size()==input.ri.size());
    assert(input.id_i.size()==input.ri.size());

    Agg_PME_prover_input pme_p_in;
    hashAndMapToG1(pme_p_in.H, "H");

    cout<<"line 52"<<endl;


    FrVec Ix = Polytree(input.id_i);
    FrVec Rx = {0};

    cout<<"line 57"<<endl;

    PCS& pcs=input.pcs;

    cout<<"line 63"<<endl;

    //PME_prover dummy input setting
    for(uint32_t i=0;i<input.id_i.size();i++){

        FrVec tmp = PolyLongDiv(Ix,{input.id_i[i],1});
        pme_p_in.w.push_back(pcs.commit_G1(tmp));

    }
    

    cout<<"line 74"<<endl;

    Fr tmp_2N = 1/(input.id_i.size());
    for(uint32_t i=0;i<input.id_i.size();i++){

        FrVec numer = PolyMul({input.ri[i]},Ix);
        FrVec denom = {input.id_i[i],1};
        Rx = PolyAdd(Rx, PolyLongDiv(numer,denom));
        
    }

    cout<<"line 83"<<endl;

    Fr r_I; r_I.setByCSPRNG();
    G2 C_I = pcs.commit(Ix) + pcs.commit_h2({r_I});
    G2 C_Ip = pcs.commit(Ix) + pcs.commit_h2(Rx);

    cout<<"line 89"<<endl;

    FrVec c;
    for(uint32_t i=0;i<input.id_i.size();i++){

        c.push_back(1/PolyEvaluate(PolyDifferentiate(Ix), -input.id_i[i]));

    }

    cout<<"line 98"<<endl;

    Fr alpha; alpha.setByCSPRNG();
    FrVec e;
    for(uint32_t i=0;i<input.id_i.size();i++){

        e.push_back(PME_V1(c[i],alpha));

    }

    cout<<"line 108"<<endl;

    vector<G1> P;
    for(uint32_t i=0;i<input.id_i.size();i++){

        P.push_back(PME_P1(e[i],pme_p_in.w[i],pme_p_in.H,input.id_i[i],(i+1)));

    }

    cout<<"line 117"<<endl;

    vector<G1> W;
    G1 w_I;w_I.clear();
    for(uint32_t i=0;i<c.size()/2;i++){

        w_I = w_I + PME_V2(P[i],P[i+1],pcs.commit(Ix),input.id_i[i],input.id_i[i+1],e[i],e[i+1],alpha,pcs);

    }

    cout<<"line 127"<<endl;

    
    zkbpacc_setup setup2;
    setup2.init(pcs.pp.g1si.size());
    setup2.g1si = &pcs.pp.g1si[0];
    setup2.g2si = &pcs.pp.g2si[0];

    ZKMP PI_zkmp(setup2);
    PI_zkmp.prove(C_I, input.C_A, w_I, r_I, 0);

    uint32_t temp = nextPowOf2(input.C_i.size());
    input.C_i.resize(temp);
    pme_p_in.w.resize(temp);

    GT E;
    pairing(E,input.pcs.pp.g1si[0]*(input.id_i.size()),C_Ip);
    IPPproof PI_ipp = zkIPPprove(input.C_i, E, pme_p_in.w);

    Fr z;
    z.setHashOf(C_I.getStr()+C_Ip.getStr());
    Fr y = PolyEvaluate(PolySub(Rx,{r_I}),z);
    G2 pi_open = pcs.Eval((C_I-C_Ip),z,y,PolySub(Rx,{r_I}));

    Agg_output out = {C_I, C_Ip, PI_zkmp, PI_ipp, pi_open};

    cout<<"line 157"<<endl;
    return out;


}

void AggZKMP_test(uint32_t n, const std::string& file_name, const std::string& file_name1, const std::string& file_name2, const std::string& file_name3, 
const std::string& file_name4, const std::string& file_name5, const std::string& file_name6){

    //
    PCS pcs;   
    FrVec IDs;
    FrVec r;
    FrVec x;

    for(uint32_t i=0;i<n;i++){

        Fr tmp;
        tmp.setByCSPRNG();
        IDs.push_back(tmp);

    }


    for(uint32_t i=0;i<n;i++){

        Fr tmp;
        tmp.setByCSPRNG();
        r.push_back(tmp);

    }
        
    pcs.setup(nextPowOf2(IDs.size()));
    
    vector<G2> C;
    FrVec tmp;
    for(uint32_t i=0;i<n;i++){
        
        tmp = Polytree({IDs[i]});
        G2 Ctmp=pcs.commit(tmp);
        C.push_back(Ctmp);
    }
    G2 C_A = pcs.commit(Polytree(IDs));
    //

    Agg_input input = {pcs, C, r, IDs, C_A};
    std::ofstream fout1(file_name, std::ios::app);
    fout1<<"N="<<n<< ", "<<endl;
    std::ofstream fout2(file_name1, std::ios::app);
    fout2<<"N="<<n<< ", "<<endl;
    std::ofstream fout3(file_name2, std::ios::app);
    fout3<<"N="<<n<< ", "<<endl;
    std::ofstream fout4(file_name3, std::ios::app);
    fout4<<"N="<<n<< ", "<<endl;
    std::ofstream fout5(file_name4, std::ios::app);
    fout5<<"N="<<n<< ", "<<endl;
    std::ofstream fout6(file_name5, std::ios::app);
    fout6<<"N="<<n<< ", "<<endl;
    std::ofstream fout7(file_name6, std::ios::app);
    fout7<<"N="<<n<< ", "<<endl;


    microseconds total(0);
    microseconds total_1(0);
    microseconds total_2(0);
    microseconds total_3(0);
    microseconds total_4(0);
    microseconds total_5(0);
    microseconds total_6(0);

    microseconds total_t(0);
    microseconds elaps;
    for(int j=0;j<100;j++){

        assert((input.C_i.size()%2==0));
        assert((input.id_i.size()%2==0));
        assert((input.ri.size()%2==0));
        

        Agg_PME_prover_input pme_p_in;
        hashAndMapToG1(pme_p_in.H, "H");


        auto start_1 = chrono::steady_clock::now();
        
        FrVec Ix = Polytree(input.id_i);
        FrVec Rx = {0};


        auto end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        total_1 = total_1 + elaps;

        start_1 = chrono::steady_clock::now();

        PCS& pcs=input.pcs;

        end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        total_1 = total_1 + elaps;
        fout2<< total_1.count()<< ", ";
        if(((j+1)%10)==0){
            fout2<<endl;
        }

        total = total + total_1;


        //PME_prover dummy input setting
        for(uint32_t i=0;i<input.id_i.size();i++){

            FrVec tmp = PolyLongDiv(Ix,{input.id_i[i],1});
            pme_p_in.w.push_back(pcs.commit_G1(tmp));

        }
    

        start_1 = chrono::steady_clock::now();

        Fr tmp_2N = 1/(input.id_i.size());
        for(uint32_t i=0;i<input.id_i.size();i++){

            FrVec numer = PolyMul({input.ri[i]},Ix);
            FrVec denom = {input.id_i[i],1};
            Rx = PolyAdd(Rx, PolyLongDiv(numer,denom));

            
        }

        end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        total_2 = total_2 + elaps;

        start_1 = chrono::steady_clock::now();

        Fr r_I; r_I.setByCSPRNG();
        G2 C_I = pcs.commit(Ix) + pcs.commit_h2({r_I});
        G2 C_Ip = pcs.commit(Ix) + pcs.commit_h2(Rx);


        end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        total_2 = total_2 + elaps;

        start_1 = chrono::steady_clock::now();        

        FrVec c;
        for(uint32_t i=0;i<input.id_i.size();i++){

            c.push_back(1/PolyEvaluate(PolyDifferentiate(Ix), -input.id_i[i]));

        }

        end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        total_2 = total_2 + elaps;

        start_1 = chrono::steady_clock::now();   

        Fr alpha; alpha.setByCSPRNG();
        FrVec e;
        for(uint32_t i=0;i<input.id_i.size();i++){

            e.push_back(PME_V1(c[i],alpha));

        }

        end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        total_2 = total_2 + elaps;


        vector<G1> P;
        for(uint32_t i=0;i<input.id_i.size();i++){

            P.push_back(PME_P1(e[i],pme_p_in.w[i],pme_p_in.H,input.id_i[i],(i+1)));

        }

        start_1 = chrono::steady_clock::now();

        vector<G1> W;
        G1 w_I;w_I.clear();
        for(uint32_t i=0;i<c.size()/2;i++){

            w_I = w_I + PME_V2(P[i],P[i+1],pcs.commit(Ix),input.id_i[i],input.id_i[i+1],e[i],e[i+1],alpha,pcs);

        }

        end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        total_2 = total_2 + elaps;
        fout3<< total_2.count()<< ", ";
        if(((j+1)%10)==0){
            fout3<<endl;
        }
        total = total + total_2;

        zkbpacc_setup setup2;
        setup2.init(pcs.pp.g1si.size());
        setup2.g1si = &pcs.pp.g1si[0];
        setup2.g2si = &pcs.pp.g2si[0];

        ZKMP PI_zkmp(setup2);

        start_1 = chrono::steady_clock::now();

        PI_zkmp.prove(C_I, input.C_A, w_I, r_I, 0);

        end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        total_3 = total_3 + elaps;
        fout4<< total_3.count()<< ", ";
        if(((j+1)%10)==0){
            fout4<<endl;
        }
        total = total + total_3;

        start_1 = chrono::steady_clock::now();

        uint32_t temp = nextPowOf2(input.C_i.size());
        input.C_i.resize(temp);
        pme_p_in.w.resize(temp);
        GT E;
        pairing(E,input.pcs.pp.g1si[0]*(input.id_i.size()),C_Ip);
        IPPproof PI_ipp = zkIPPprove(input.C_i, E, pme_p_in.w);

        end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        total_4 = total_4 + elaps;
        fout5<< total_4.count()<< ", ";
        if(((j+1)%10)==0){
            fout5<<endl;
        }
        total = total + total_4;

        start_1 = chrono::steady_clock::now();

        Fr z;
        z.setHashOf(C_I.getStr()+C_Ip.getStr());
        Fr y = PolyEvaluate(PolySub(Rx,{r_I}),z);
        G2 pi_open = pcs.Eval((C_I-C_Ip),z,y,PolySub(Rx,{r_I}));

        end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        total_5 = total_5 + elaps;
        fout6<< total_5.count()<< ", ";
        if(((j+1)%10)==0){
            fout6<<endl;
        }
        total = total + total_5;

        start_1 = chrono::steady_clock::now();

        Agg_output out = {C_I, C_Ip, PI_zkmp, PI_ipp, pi_open};

        end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        fout7<< elaps.count()<< ", ";
        if(((j+1)%10)==0){
            fout7<<endl;
        }
        total_6 = total_6 + elaps;      
        fout1<< total.count()<< ", ";
        if(((j+1)%10)==0){
            fout1<<endl;
        }
        total_t = total_t + total;
        
        total = microseconds(0);
        total_1 = microseconds(0);
        total_2 = microseconds(0);
        total_3 = microseconds(0);
        total_4 = microseconds(0);
        total_5 = microseconds(0);
        total_6 = microseconds(0);

    }
    double tt = total_t.count();
    fout1<<endl;

}

int main(){

    initPairing(mcl::BLS12_381);
    
    //     struct Agg_input {

    //     PCS pcs;
    //     vector<G2> C_i;
    //     FrVec ri;
    //     FrVec id_i;
    //     G2 C_A;
        
    // };

    // PCS pcs;
    
    // FrVec IDs;
    // FrVec r;
    // FrVec x;

    // for(uint32_t i=0;i<1000;i++){

    //     Fr tmp;
    //     tmp.setByCSPRNG();
    //     IDs.push_back(tmp);

    // }


    // for(uint32_t i=0;i<1000;i++){

    //     Fr tmp;
    //     tmp.setByCSPRNG();
    //     r.push_back(tmp);

    // }

    // pcs.setup(nextPowOf2(IDs.size()));
    // G2 A = pcs.commit(Polytree(IDs));
    
    // vector<G2> C;
    // for(uint32_t i=0;i<IDs.size();i++){
    //     FrVec tmp;
    //     tmp = Polytree({IDs[i]});
    //     G2 Ctmp=pcs.commit(tmp);
    //     C.push_back(Ctmp);
    // }
    
    // G2 C_A = pcs.commit(Polytree(IDs));

    // Agg_input in = {pcs, C, r, IDs, C_A};

    // AggZKMP(in);

    std::string file_name = "Agg_ZKMP.csv";
    std::string file_name1 = "Agg_ZKMP_commit.csv";     
    std::string file_name2 = "Agg_ZKMP_PME.csv";  
    std::string file_name3 = "Agg_ZKMP_ZKMP.csv";  
    std::string file_name4 = "Agg_ZKMP_ZKIPP.csv";  
    std::string file_name5 = "Agg_ZKMP_PCopenP.csv";  
    std::string file_name6 = "Agg_ZKMP_StructGen.csv";  
    
    
    AggZKMP_test(2,file_name, file_name1, file_name2, file_name3, file_name4, file_name5, file_name6);
    AggZKMP_test(256,file_name, file_name1, file_name2, file_name3, file_name4, file_name5, file_name6);
    AggZKMP_test(512,file_name, file_name1, file_name2, file_name3, file_name4, file_name5, file_name6);
    AggZKMP_test(1024,file_name, file_name1, file_name2, file_name3, file_name4, file_name5, file_name6);
    AggZKMP_test(2048,file_name, file_name1, file_name2, file_name3, file_name4, file_name5, file_name6);

}