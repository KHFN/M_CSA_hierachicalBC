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

struct output{
        Fr Adr;
        G2 C;//ID commitment
        G2 A;//Serial Number accumulator
        ZKMP PI_PoA;
};

struct Entry_SN{
        FrVec M;
        PoK_proof_g2 PI_PoK;
};

struct Tx_Entry{

    Fr Tid;

    output output1;
    Entry_SN SN;
   
};

Tx_Entry Tx_Entry_Gen(PCS pcs, ZKMP PI_PoA, Fr id, G1 wit, G2 C_A, FrVec SNs, Fr gamma, Fr delta, const string& file_name){

    std::ofstream fout1(file_name, std::ios::app);

    microseconds ID_comm(0);
    microseconds Acc_SN(0);
    microseconds ID_prove(0);
    microseconds PoKE(0);
    microseconds PoQ(0);
    microseconds Hash(0);
    microseconds elaps;

    G2 C = pcs.commit(Polytree({id}));
    auto start_1 = chrono::steady_clock::now();
    C = C + pcs.pp.h2si[0]*gamma;
    auto end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    ID_comm = ID_comm + elaps;

    start_1 = chrono::steady_clock::now();
    G2 A = pcs.commit(Polytree(SNs)) + pcs.pp.h2si[0]*delta;
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    Acc_SN = Acc_SN + elaps;


    start_1 = chrono::steady_clock::now();
    PI_PoA.prove(C, C_A, wit, gamma, 0);
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    ID_prove = ID_prove + elaps;

    start_1 = chrono::steady_clock::now();
    string buf = PI_PoA.msg.P_1.getStr()+PI_PoA.msg.P_2.getStr()+PI_PoA.msg.R_1.getStr()+PI_PoA.msg.R_2.getStr()+PI_PoA.msg.R_3.getStr()
    +PI_PoA.response.s_delta_1.getStr()+PI_PoA.response.s_delta_2.getStr()+PI_PoA.response.s_r_A.getStr()+PI_PoA.response.s_r_I.getStr()
    +PI_PoA.response.s_tau_1.getStr();
    buf = buf+C.getStr()+A.getStr();
    Fr Adr; Adr.setHashOf(buf);
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    Hash = Hash + elaps;

    
    output output1 = {Adr, C, A, PI_PoA};

    start_1 = chrono::steady_clock::now();
    PoK_proof_g2 PI_PoK;
    PI_PoK.prove(pcs.pp.g2si[0]*delta,delta);
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    PoKE = PoKE + elaps;

    start_1 = chrono::steady_clock::now();
    Fr Tid;
    Tid.setByCSPRNG();
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    Hash = Hash + elaps;

    Entry_SN SN = {SNs, PI_PoK};
    Tx_Entry Entry = {Tid, output1, SN};

    fout1<< ID_comm.count()<< ", ";
    fout1<< Acc_SN.count()<< ", ";
    fout1<< ID_prove.count()<< ", ";
    fout1<< PoKE.count()<< ", ";
    fout1<< PoQ.count()<< ", ";
    fout1<< Hash.count()<<", "<<endl;

    return Entry;

}


struct Transfer_input{

    Fr Adr;
    G2 C_I;
    G2 A_I;
    PoK_proof_g2_n PI_PoKE;

};

struct Tx_Transfer{

    Fr Tid;

    Transfer_input input;
    output output1;
    output output2;
    ZKSP PI_PoQ;

};

Tx_Transfer Tx_Transfer_Gen(PCS pcs, ZKMP PI_PoA_1, ZKMP PI_PoA_2, ZKSP PI_PoQ, Fr id1, Fr id2, G1 wit1, G1 wit2, G2 C_A, FrVec SNs, FrVec S_1, FrVec S_2, output output_old, 
Fr gamma, Fr gamma1, Fr gamma2, Fr delta, Fr delta1, Fr delta2, const std::string& file_name, Fr Tid_old){

    std::ofstream fout1(file_name, std::ios::app);

    microseconds ID_comm(0);
    microseconds Acc_SN(0);
    microseconds ID_prove(0);
    microseconds PoKE(0);
    microseconds PoQ(0);
    microseconds Hash(0);
    microseconds elaps;

    auto start_1 = chrono::steady_clock::now();
    PoK_proof_g2_n PI_PoKE;
    PI_PoKE.prove({pcs.pp.g2si[0],pcs.pp.h2si[0]}, output_old.C-pcs.pp.g2si[1],{id1,gamma},2);
    auto end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    ID_prove = ID_prove + elaps;

    Transfer_input input = {output_old.Adr, output_old.C, output_old.A, PI_PoKE};
    
    G2 C_1 = pcs.commit(Polytree({id1}));
    start_1 = chrono::steady_clock::now();
    C_1 = C_1 + pcs.pp.h2si[0]*gamma1;
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    ID_comm = ID_comm + elaps;

    start_1 = chrono::steady_clock::now();
    G2 A_1 = pcs.commit(Polytree(S_1)) + pcs.pp.h2si[0]*delta1;
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    Acc_SN = Acc_SN + elaps;


    start_1 = chrono::steady_clock::now();
    PI_PoA_1.prove(C_1, C_A, wit1, gamma1, 0);
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    ID_prove = ID_prove + elaps;


    start_1 = chrono::steady_clock::now();
    string buf = PI_PoA_1.msg.P_1.getStr()+PI_PoA_1.msg.P_2.getStr()+PI_PoA_1.msg.R_1.getStr()+PI_PoA_1.msg.R_2.getStr()+PI_PoA_1.msg.R_3.getStr()
    +PI_PoA_1.response.s_delta_1.getStr()+PI_PoA_1.response.s_delta_2.getStr()+PI_PoA_1.response.s_r_A.getStr()+PI_PoA_1.response.s_r_I.getStr()
    +PI_PoA_1.response.s_tau_1.getStr();
    buf = buf + C_1.getStr() + A_1.getStr();
    Fr Adr1; 
    Adr1.setHashOf(buf);
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    Hash = Hash + elaps;

    output output1 = {Adr1,C_1,A_1,PI_PoA_1};

    G2 C_2 = pcs.commit(Polytree({id2}));
    start_1 = chrono::steady_clock::now();
    C_2 = C_2 + pcs.pp.h2si[0]*gamma2;
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    ID_comm = ID_comm + elaps;

    start_1 = chrono::steady_clock::now();
    G2 A_2 = pcs.commit(Polytree(S_2)) + pcs.pp.h2si[0]*delta2;
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    Acc_SN = Acc_SN + elaps;

    start_1 = chrono::steady_clock::now();
    PI_PoA_2.prove(C_2, C_A, wit2, gamma2, 0);
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    ID_prove = ID_prove + elaps;

    start_1 = chrono::steady_clock::now();
    buf = PI_PoA_2.msg.P_1.getStr() + PI_PoA_2.msg.P_2.getStr() + PI_PoA_2.msg.R_1.getStr() + PI_PoA_2.msg.R_2.getStr() + PI_PoA_2.msg.R_3.getStr()
    +PI_PoA_2.response.s_delta_1.getStr()+PI_PoA_2.response.s_delta_2.getStr()+PI_PoA_2.response.s_r_A.getStr()+PI_PoA_2.response.s_r_I.getStr()
    +PI_PoA_2.response.s_tau_1.getStr();
    buf = buf + C_2.getStr() + A_2.getStr();
    Fr Adr2; 
    Adr2.setHashOf(buf);
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    Hash = Hash + elaps;


    output output2 ={Adr2,C_2,A_2,PI_PoA_2};

    start_1 = chrono::steady_clock::now();
    FrVec S;
    S.insert(S.end(),S_1.begin(),S_1.end());
    S.insert(S.end(),S_2.begin(),S_2.end());
    G1 W_1 = pcs.commit_G1(PolyLongDiv(Polytree(S),Polytree(S_1)));
    G1 W_2 = pcs.commit_G1(PolyLongDiv(Polytree(S),Polytree(S_2)));
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    Acc_SN = Acc_SN + elaps;


    start_1 = chrono::steady_clock::now();
    PI_PoQ.prove(A_1, A_2, output_old.A, W_1, W_2, delta, delta1, delta2);
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    PoQ = PoQ + elaps;

    fout1<< ID_comm.count()<< ", ";
    fout1<< Acc_SN.count()<< ", ";
    fout1<< ID_prove.count()<< ", ";
    fout1<< PoKE.count()<< ", ";
    fout1<< PoQ.count()<< ", ";
    fout1<< Hash.count()<<", "<<endl;
   

    Fr Tid;
    Tid.setHashOf(Tid_old.getStr());

    Tx_Transfer Tx = {Tid,input,output1,output2,PI_PoQ};

    return Tx;

}




struct Tx_Exit{

    Transfer_input input;
    output output_;
    FrVec N;
    ZKSP PI_PoQ;

};


Tx_Exit Tx_Exit_Gen(PCS pcs, output output_old, ZKMP PI_PoA, ZKSP PI_PoQ, Fr id, G1 wit, G2 C_A, Fr gamma, Fr gamma1, Fr delta, Fr delta1, FrVec S, FrVec M, FrVec N, const std::string& file_name){

    std::ofstream fout1(file_name, std::ios::app);

    microseconds ID_comm(0);
    microseconds Acc_SN(0);
    microseconds ID_prove(0);
    microseconds PoKE(0);
    microseconds PoQ(0);
    microseconds Hash(0);
    microseconds elaps;

    auto start_1 = chrono::steady_clock::now();
    G2 C_I = output_old.C;
    auto end_1 = chrono::steady_clock::now();

    start_1 = chrono::steady_clock::now();
    G2 A_I = output_old.A;
    end_1 = chrono::steady_clock::now();

    start_1 = chrono::steady_clock::now();
    PoK_proof_g2_n PI_PoK;
    PI_PoK.prove({pcs.pp.g2si[0],pcs.pp.h2si[0]}, C_I-pcs.pp.g2si[1], {id, gamma}, 2);
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    ID_prove = ID_prove + elaps;

    Transfer_input input = {output_old.Adr, C_I, A_I, PI_PoK};

    G2 C_1 = pcs.commit(Polytree({id}));
    start_1 = chrono::steady_clock::now();
    C_1 = C_1 + pcs.pp.h2si[0]*gamma1;
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    ID_comm = ID_comm + elaps;

    start_1 = chrono::steady_clock::now();
    G2 A_1 = pcs.commit(Polytree(N)) + pcs.pp.h2si[0]*delta1;
    G2 A_2 = pcs.commit(Polytree(M));
    G1 W_1 = pcs.commit_G1(PolyLongDiv(Polytree(S),Polytree(M)));
    G1 W_2 = pcs.commit_G1(PolyLongDiv(Polytree(S),Polytree(N)));
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    Acc_SN = Acc_SN + elaps;

    start_1 = chrono::steady_clock::now();
    PI_PoA.prove(C_1,C_A,wit,gamma1,0);
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    ID_prove = ID_prove + elaps;

    start_1 = chrono::steady_clock::now();
    string buf = PI_PoA.msg.P_1.getStr()+PI_PoA.msg.P_2.getStr()+PI_PoA.msg.R_1.getStr()+PI_PoA.msg.R_2.getStr()+PI_PoA.msg.R_3.getStr()
    +PI_PoA.response.s_delta_1.getStr()+PI_PoA.response.s_delta_2.getStr()+PI_PoA.response.s_r_A.getStr()+PI_PoA.response.s_r_I.getStr()
    +PI_PoA.response.s_tau_1.getStr();
    buf = buf+C_1.getStr()+A_1.getStr();
    Fr Adr; Adr.setHashOf(buf);
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    Hash = Hash + elaps;


    output output_ = {Adr, C_1, A_1, PI_PoA};

    start_1 = chrono::steady_clock::now();
    PI_PoQ.prove(A_1,A_2,A_I,W_1,W_2,delta1,0,delta);
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    PoQ = PoQ + elaps;

    

    Tx_Exit Tx = {input, output_, N, PI_PoQ};

    fout1<< ID_comm.count()<< ", ";
    fout1<< Acc_SN.count()<< ", ";
    fout1<< ID_prove.count()<< ", ";
    fout1<< PoKE.count()<< ", ";
    fout1<< PoQ.count()<< ", ";
    fout1<< Hash.count()<<", "<<endl;

    return Tx;

}

struct Transfer_out{

    G2 com_LocalID;
    G2 com_SN_t;

};

struct Tx_Transfer_out{

    Transfer_out IN;
    Transfer_out OUT1;
    Transfer_out OUT2;

};

Tx_Transfer_out Tx_Transfer_out_Gen(PCS pcs, Fr chainID1, Fr chainID2, Tx_Exit Tx1, Tx_Entry Tx2, Fr gamma, Fr t2, Fr u1, Fr mu1, Fr mu2, const std::string& file_name){

    std::ofstream fout1(file_name, std::ios::app);

    microseconds ID_comm(0);
    microseconds Acc_SN(0);
    microseconds elaps;

    G2 C_1 = pcs.commit(Polytree({chainID1}));
    auto start_1 = chrono::steady_clock::now();
    C_1 = C_1+pcs.pp.h2si[0]*gamma;
    auto end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    ID_comm = ID_comm + elaps;

    start_1 = chrono::steady_clock::now();
    G2 ACC_1 = pcs.commit(Polytree(Tx1.N))+pcs.pp.h2si[0]*mu1;
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    Acc_SN = Acc_SN + elaps;

    Transfer_out IN = {C_1, ACC_1};

    G2 C_1_1 = pcs.commit(Polytree({chainID1}));
    start_1 = chrono::steady_clock::now();
    C_1_1 = C_1_1 +pcs.pp.h2si[0]*t2;
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    ID_comm = ID_comm + elaps;

    start_1 = chrono::steady_clock::now();
    G2 ACC_2 = pcs.commit(Polytree(Tx2.SN.M))+pcs.pp.h2si[0]*mu1;
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    Acc_SN = Acc_SN + elaps;

    Transfer_out OUT1 = {C_1_1, ACC_2};

    G2 C_2 = pcs.commit(Polytree({chainID2}));
    start_1 = chrono::steady_clock::now();
    C_2 = C_2+pcs.pp.h2si[0]*u1;
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    ID_comm = ID_comm + elaps;

    start_1 = chrono::steady_clock::now();
    G2 ACC_3 = pcs.commit(PolyLongDiv(Polytree(Tx1.N),Polytree(Tx2.SN.M)))+pcs.pp.h2si[0]*mu2;
    end_1 = chrono::steady_clock::now();
    elaps = duration_cast<microseconds>(end_1-start_1);
    Acc_SN = Acc_SN + elaps;

    Transfer_out OUT2 = {C_2,ACC_3};
    Tx_Transfer_out Tx = {IN, OUT1, OUT2};

    fout1<< ID_comm.count()<< ", ";
    fout1<< Acc_SN.count();
    fout1<<endl;


    return Tx;

}

// void Tx_Entry_Gen_test(uint32_t n, const std::string& file_name){

//     PCS pcs;
//     pcs.setup(100001);
//     zkbpacc_setup setup;
//     setup.init(100001);
//     setup.g1si = &pcs.pp.g1si[0];
//     setup.g2si = &pcs.pp.g2si[0];

//     ZKMP PI_PoA(setup);
//     FrVec IDset;
//     for(uint32_t i=0;i<100000;i++){
//         Fr id;id.setByCSPRNG();
//         IDset.push_back(id);
//     }
//     srand(time(0));
//     uint32_t j = rand()%100000;
//     Fr id = IDset[j];

//     G1 wit = pcs.commit_G1(PolyLongDiv(Polytree(IDset),Polytree({id})));
//     G2 C_A = pcs.commit(Polytree(IDset));
//     FrVec SNs;
//     for(uint32_t i=0;i<n;i++){
//         Fr SN;SN.setByCSPRNG();
//         SNs.push_back(SN);
//     }

//     Fr gamma; gamma.setByCSPRNG();
//     Fr delta; delta.setByCSPRNG();

//     std::ofstream fout1(file_name, std::ios::app);
//     fout1<<"N="<<n<< ", "<<endl;

//     microseconds ID_comm(0);
//     microseconds Acc_SN(0);
//     microseconds ID_prove(0);
//     microseconds PoKE(0);
//     microseconds Hash(0);
//     microseconds elaps;

//     auto start_1 = chrono::steady_clock::now();
//     G2 C = pcs.commit(Polytree({id})) + pcs.pp.h2si[0]*gamma;
//     auto end_1 = chrono::steady_clock::now();
//     elaps = duration_cast<microseconds>(end_1-start_1);
//     ID_comm = ID_comm + elaps;

//     start_1 = chrono::steady_clock::now();
//     G2 A = pcs.commit(SNs) + pcs.pp.h2si[0]*delta;
//     end_1 = chrono::steady_clock::now();
//     elaps = duration_cast<microseconds>(end_1-start_1);
//     Acc_SN = Acc_SN + elaps;

//     start_1 = chrono::steady_clock::now();
//     PI_PoA.prove(C, C_A, wit, gamma, 0);
//     end_1 = chrono::steady_clock::now();
//     elaps = duration_cast<microseconds>(end_1-start_1);
//     ID_prove = ID_prove + elaps;

//     start_1 = chrono::steady_clock::now();
//     string buf = PI_PoA.msg.P_1.getStr()+PI_PoA.msg.P_2.getStr()+PI_PoA.msg.R_1.getStr()+PI_PoA.msg.R_2.getStr()+PI_PoA.msg.R_3.getStr()
//     +PI_PoA.response.s_delta_1.getStr()+PI_PoA.response.s_delta_2.getStr()+PI_PoA.response.s_r_A.getStr()+PI_PoA.response.s_r_I.getStr()
//     +PI_PoA.response.s_tau_1.getStr();
//     buf = buf+C.getStr()+A.getStr();
//     Fr Adr; Adr.setHashOf(buf);
//     end_1 = chrono::steady_clock::now();
//     elaps = duration_cast<microseconds>(end_1-start_1);
//     Hash = Hash + elaps;

//     output output1 = {Adr, C, A, PI_PoA};

//     start_1 = chrono::steady_clock::now();
//     PoK_proof_g2 PI_PoK;
//     PI_PoK.prove(pcs.pp.g2si[0]*delta,delta);
//     end_1 = chrono::steady_clock::now();
//     elaps = duration_cast<microseconds>(end_1-start_1);
//     PoKE = PoKE + elaps;

//     fout1<< ID_comm.count()<< ", ";
//     fout1<< Acc_SN.count()<< ", ";
//     fout1<< ID_prove.count()<< ", ";
//     fout1<< PoKE.count()<< ", ";
//     fout1<<Hash.count()<<endl;
//     fout1<<endl;

//     Entry_SN SN = {SNs, PI_PoK};
//     Tx_Entry Entry = {output1, SN};

//     fout1<<"Tx_Transfer"<<endl;
//     fout1<<"N="<<n<< ", "<<endl;

//     ZKMP PI_PoA_1(setup);
//     ZKMP PI_PoA_2(setup);
//     ZKSP PI_PoQ(setup);
//     uint32_t m = rand()%10;
//     Fr id2 = IDset[m];

//     G1 wit2 = pcs.commit_G1(PolyLongDiv(Polytree(IDset),Polytree({id2})));

//     vector<FrVec> outVec;
//     uint32_t ntmp = SNs.size() / 5;
//     FrVec S_2(SNs.begin(), SNs.begin() + ntmp);
//     FrVec S_1(SNs.begin() + ntmp, SNs.end());
//     Fr gamma1;gamma1.setByCSPRNG();
//     Fr gamma2;gamma2.setByCSPRNG();

//     Fr delta1;delta1.setByCSPRNG();
//     Fr delta2;delta2.setByCSPRNG();

//     Tx_Transfer Tx1 = Tx_Transfer_Gen(pcs, PI_PoA_1, PI_PoA_2, PI_PoQ, id, id2, wit, wit2, C_A, SNs, S_1, S_2, Entry.output1, 
//     gamma, gamma1, gamma2, delta, delta1, delta2, file_name);

//     fout1<<"Tx_Exit"<<endl;
//     fout1<<"N="<<n<< ", "<<endl;

//     ZKMP PI_PoA_3(setup);
//     ZKSP PI_PoQ_2(setup);
//     gamma1.setByCSPRNG();
//     delta1.setByCSPRNG();

//     ntmp = S_2.size() / 5;
//     FrVec M(S_2.begin(), S_2.begin() + ntmp);
//     FrVec N(S_2.begin() + ntmp, S_2.end());

//     Tx_Exit Tx2 = Tx_Exit_Gen(pcs, Tx1.output2, PI_PoA_3, PI_PoQ_2, id, wit, C_A, gamma, gamma1, delta, delta1, SNs, S_1, S_2, file_name);

//     // fout1<<"Tx_inter"<<endl;
//     // fout1<<"N="<<n<< ", "<<endl;

//     // Fr chainID1; chainID1.setByCSPRNG();
//     // Fr chainID2; chainID2.setByCSPRNG();
//     // gamma.setByCSPRNG();
//     // delta.setByCSPRNG();

//     // ZKMP PI_PoA6(setup);

//     // Tx_Entry Tx3 = Tx_Entry_Gen(pcs, PI_PoA6, id, wit, C_A, S_1, gamma, delta);

//     // Fr t2;t2.setByCSPRNG();
//     // Fr u1;u1.setByCSPRNG(); 
//     // Fr mu1;mu1.setByCSPRNG();  
//     // Fr mu2;mu2.setByCSPRNG();

//     // Tx_Transfer_out_Gen(pcs, chainID1, chainID2, Tx2, Entry, gamma, t2, u1, mu1, mu2, file_name);

// }



int main(){
    
    initPairing(mcl::BLS12_381);

    PCS pcs;
    pcs.setup(2001);
    zkbpacc_setup setup;
    setup.init(2001);
    setup.g1si = &pcs.pp.g1si[0];
    setup.g2si = &pcs.pp.g2si[0];

    ZKMP PI_PoA(setup);
    FrVec IDset;
    for(uint32_t i=0;i<10;i++){
        Fr id;id.setByCSPRNG();
        IDset.push_back(id);
    }
    srand(time(0));
    uint32_t j = rand()%10;
    Fr id = IDset[j];

    const std::string& file_name = "EntryTxGen.csv"; 
    ZKMP PI_PoA1(setup);

    G1 wit = pcs.commit_G1(PolyLongDiv(Polytree(IDset),Polytree({id})));
    G2 C_A = pcs.commit(Polytree(IDset));

    FrVec SNs;
    for(uint32_t i=0;i<125;i++){
        Fr SN;SN.setByCSPRNG();
        SNs.push_back(SN);
    }

    FrVec SNs2;
    for(uint32_t i=0;i<250;i++){
        Fr SN;SN.setByCSPRNG();
        SNs2.push_back(SN);
    }

    FrVec SNs3;
    for(uint32_t i=0;i<500;i++){
        Fr SN;SN.setByCSPRNG();
        SNs3.push_back(SN);
    }

    FrVec SNs4;
    for(uint32_t i=0;i<1000;i++){
        Fr SN;SN.setByCSPRNG();
        SNs4.push_back(SN);
    }

    FrVec SNs5;
    for(uint32_t i=0;i<2000;i++){
        Fr SN;SN.setByCSPRNG();
        SNs5.push_back(SN);
    }

    Fr gamma; gamma.setByCSPRNG();
    Fr delta; delta.setByCSPRNG();

    // for(uint8_t i=0;i<100;i++){
    
    //     Tx_Entry Tx1 = Tx_Entry_Gen(pcs, PI_PoA1, id, wit, C_A, SNs, gamma, delta, "Tx_Entry_625.csv");
    
    // }

    // for(uint8_t i=0;i<100;i++){

    //     Tx_Entry Tx2 = Tx_Entry_Gen(pcs, PI_PoA1, id, wit, C_A, SNs2, gamma, delta, "Tx_Entry_1250.csv");

    // }

    // for(uint8_t i=0;i<100;i++){

    //     Tx_Entry Tx3 = Tx_Entry_Gen(pcs, PI_PoA1, id, wit, C_A, SNs3, gamma, delta, "Tx_Entry_2500.csv");

    // }

    // for(uint8_t i=0;i<100;i++){

    //     Tx_Entry Tx4 = Tx_Entry_Gen(pcs, PI_PoA1, id, wit, C_A, SNs4, gamma, delta, "Tx_Entry_5000.csv");

    // }


        Tx_Entry Tx1 = Tx_Entry_Gen(pcs, PI_PoA1, id, wit, C_A, SNs, gamma, delta, "Tx_Entry_125.csv");
        Tx_Entry Tx2 = Tx_Entry_Gen(pcs, PI_PoA1, id, wit, C_A, SNs2, gamma, delta, "Tx_Entry_250.csv");
        Tx_Entry Tx3 = Tx_Entry_Gen(pcs, PI_PoA1, id, wit, C_A, SNs3, gamma, delta, "Tx_Entry_500.csv");
        Tx_Entry Tx4 = Tx_Entry_Gen(pcs, PI_PoA1, id, wit, C_A, SNs4, gamma, delta, "Tx_Entry_1000_tmp.csv");   
        Tx_Entry Tx5 = Tx_Entry_Gen(pcs, PI_PoA1, id, wit, C_A, SNs5, gamma, delta, "Tx_Entry_20000_tmp.csv");

 
    // int split = SNs5.size()*(1/5);

    

    ZKMP PI_PoA_1(setup);
    ZKMP PI_PoA_2(setup);
    ZKSP PI_PoQ(setup);

    uint32_t k = rand()%10;
    Fr id2 = IDset[k];
    G1 wit2 = pcs.commit_G1(PolyLongDiv(Polytree(IDset),Polytree({id2})));
    Fr gamma1; gamma1.setByCSPRNG();
    Fr gamma2; gamma2.setByCSPRNG();
    Fr delta1; delta1.setByCSPRNG();
    Fr delta2; delta2.setByCSPRNG();

    FrVec S_2(SNs.begin(),SNs.begin()+SNs.size()*(1/5));
    FrVec S_1(SNs.begin()+SNs.size()*(1/5),SNs.end());
    
    FrVec S_4(SNs2.begin(),SNs2.begin()+SNs2.size()*(1/5));
    FrVec S_3(SNs2.begin()+SNs2.size()*(1/5),SNs2.end());

    FrVec S_6(SNs3.begin(),SNs3.begin()+SNs3.size()*(1/5));
    FrVec S_5(SNs3.begin()+SNs3.size()*(1/5),SNs3.end());

    FrVec S_8(SNs4.begin(),SNs4.begin()+SNs4.size()*(1/5));
    FrVec S_7(SNs4.begin()+SNs4.size()*(1/5),SNs4.end());

    FrVec S_10(SNs5.begin(),SNs5.begin()+SNs5.size()*(1/5));
    FrVec S_9(SNs5.begin()+SNs5.size()*(1/5),SNs5.end());




        // Tx_Transfer Tx6 = Tx_Transfer_Gen(pcs, PI_PoA_1, PI_PoA_2, PI_PoQ, id, id2, wit, wit2, C_A, SNs, 
        // S_1, S_2, Tx1.output1, gamma, gamma1, gamma2, delta, delta1, delta2, "Transfer_125.csv", Tx1.Tid);

        // Tx_Transfer Tx7 = Tx_Transfer_Gen(pcs, PI_PoA_1, PI_PoA_2, PI_PoQ, id, id2, wit, wit2, C_A, SNs2, 
        // S_3, S_4, Tx2.output1, gamma, gamma1, gamma2, delta, delta1, delta2, "Transfer_250.csv", Tx2.Tid);

        // Tx_Transfer Tx8 = Tx_Transfer_Gen(pcs, PI_PoA_1, PI_PoA_2, PI_PoQ, id, id2, wit, wit2, C_A, SNs3, 
        // S_5, S_6, Tx3.output1, gamma, gamma1, gamma2, delta, delta1, delta2, "Transfer_500.csv", Tx3.Tid);

        // Tx_Transfer Tx9 = Tx_Transfer_Gen(pcs, PI_PoA_1, PI_PoA_2, PI_PoQ, id, id2, wit, wit2, C_A, SNs4, 
        // S_7, S_8, Tx4.output1, gamma, gamma1, gamma2, delta, delta1, delta2, "Transfer_1000.csv", Tx4.Tid);

        // Tx_Transfer Tx10 = Tx_Transfer_Gen(pcs, PI_PoA_1, PI_PoA_2, PI_PoQ, id, id2, wit, wit2, C_A, SNs5, 
        // S_9, S_10, Tx5.output1, gamma, gamma1, gamma2, delta, delta1, delta2, "Transfer_2000.csv", Tx5.Tid);



    // // FrVec S_2_2(SNs.begin(),SNs.begin()+SNs.size()*(2/5));
    // // FrVec S_1_2(SNs.begin()+SNs.size()*(2/5),SNs.end());
    // // FrVec S_2_3(SNs.begin(),SNs.begin()+SNs.size()*(3/5));
    // // FrVec S_1_3(SNs.begin()+SNs.size()*(3/5),SNs.end());

        Tx_Transfer Tx6 = Tx_Transfer_Gen(pcs, PI_PoA_1, PI_PoA_2, PI_PoQ, id, id2, wit, wit2, C_A, SNs, 
        {}, SNs, Tx1.output1, gamma, gamma1, gamma2, delta, delta1, delta2, "Transfer_625.csv", Tx1.Tid);

        Tx_Transfer Tx7 = Tx_Transfer_Gen(pcs, PI_PoA_1, PI_PoA_2, PI_PoQ, id, id2, wit, wit2, C_A, SNs2, 
        {}, SNs2, Tx2.output1, gamma, gamma1, gamma2, delta, delta1, delta2, "Transfer_1250.csv", Tx2.Tid);

        Tx_Transfer Tx8 = Tx_Transfer_Gen(pcs, PI_PoA_1, PI_PoA_2, PI_PoQ, id, id2, wit, wit2, C_A, SNs3, 
        {}, SNs3, Tx3.output1, gamma, gamma1, gamma2, delta, delta1, delta2, "Transfer_2500.csv", Tx3.Tid);

        Tx_Transfer Tx9 = Tx_Transfer_Gen(pcs, PI_PoA_1, PI_PoA_2, PI_PoQ, id, id2, wit, wit2, C_A, SNs4, 
        {}, SNs4, Tx4.output1, gamma, gamma1, gamma2, delta, delta1, delta2, "Transfer_5000.csv", Tx4.Tid);

        Tx_Transfer Tx10 = Tx_Transfer_Gen(pcs, PI_PoA_1, PI_PoA_2, PI_PoQ, id, id2, wit, wit2, C_A, SNs5, 
        {}, SNs5, Tx5.output1, gamma, gamma1, gamma2, delta, delta1, delta2, "Transfer_10000.csv", Tx5.Tid);

    for(int i=0;i<100;i++){

        Tx_Exit Tx11 = Tx_Exit_Gen(pcs,Tx6.output1, PI_PoA, PI_PoQ, id2, wit2, C_A, gamma1, gamma2, delta1, delta2, SNs, {}, SNs, "Exit_125.csv");
        Tx_Exit Tx12 = Tx_Exit_Gen(pcs,Tx7.output1, PI_PoA, PI_PoQ, id2, wit2, C_A, gamma1, gamma2, delta1, delta2, SNs2, {}, SNs2, "Exit_250.csv");
        Tx_Exit Tx13 = Tx_Exit_Gen(pcs,Tx8.output1, PI_PoA, PI_PoQ, id2, wit2, C_A, gamma1, gamma2, delta1, delta2, SNs3, {}, SNs3, "Exit_500.csv");
        Tx_Exit Tx14 = Tx_Exit_Gen(pcs,Tx9.output1, PI_PoA, PI_PoQ, id2, wit2, C_A, gamma1, gamma2, delta1, delta2, SNs4, {}, SNs4, "Exit_1000.csv");
        Tx_Exit Tx15 = Tx_Exit_Gen(pcs,Tx10.output1, PI_PoA, PI_PoQ, id2, wit2, C_A, gamma1, gamma2, delta1, delta2, SNs5, {}, SNs5, "Exit_2000.csv");

    }
    // ZKMP PI_PoA2(setup);
    // Fr gamma6; gamma6.setByCSPRNG();
    // Fr delta6; delta6.setByCSPRNG();

    // gamma.setByCSPRNG();
    // Fr t2;t2.setByCSPRNG();
    // Fr u1;u1.setByCSPRNG();
    // Fr mu1;mu1.setByCSPRNG();
    // Fr mu2;mu2.setByCSPRNG();

    // Tx_Entry Tx16 = Tx_Entry_Gen(pcs,PI_PoA2, id, wit, C_A, SNs, gamma6, delta6, "Temp.csv");
    // Tx_Entry Tx17 = Tx_Entry_Gen(pcs,PI_PoA2, id, wit, C_A, SNs, gamma6, delta6, "Temp.csv");
    // Tx_Entry Tx18 = Tx_Entry_Gen(pcs,PI_PoA2, id, wit, C_A, SNs, gamma6, delta6, "Temp.csv");
    // Tx_Entry Tx19 = Tx_Entry_Gen(pcs,PI_PoA2, id, wit, C_A, SNs, gamma6, delta6, "Temp.csv");
    // Tx_Entry Tx20 = Tx_Entry_Gen(pcs,PI_PoA2, id, wit, C_A, SNs, gamma6, delta6, "Temp.csv");
    // Fr chainID1;Fr chainID2;
    // chainID1.setByCSPRNG();chainID2.setByCSPRNG();

    // for(uint8_t i=0;i<100;i++){

    //     Tx_Transfer_out Tx16_1 = Tx_Transfer_out_Gen(pcs,chainID1,chainID2,Tx11,Tx16,gamma,t2,u1,mu1,mu2,"Tout_625.csv");
    //     Tx_Transfer_out Tx17_1 = Tx_Transfer_out_Gen(pcs,chainID1,chainID2,Tx12,Tx17,gamma,t2,u1,mu1,mu2,"Tout_1250.csv");
    //     Tx_Transfer_out Tx18_1 = Tx_Transfer_out_Gen(pcs,chainID1,chainID2,Tx13,Tx18,gamma,t2,u1,mu1,mu2,"Tout_2500.csv");
    //     Tx_Transfer_out Tx19_1 = Tx_Transfer_out_Gen(pcs,chainID1,chainID2,Tx14,Tx19,gamma,t2,u1,mu1,mu2,"Tout_5000.csv");
    //     Tx_Transfer_out Tx20_1 = Tx_Transfer_out_Gen(pcs,chainID1,chainID2,Tx15,Tx20,gamma,t2,u1,mu1,mu2,"Tout_10000.csv");

    // }
    

}