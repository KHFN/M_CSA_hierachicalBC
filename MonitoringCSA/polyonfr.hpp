#ifndef POLYONFR_H__
#define POLYONFR_H__

#include <mcl/bls12_381.hpp>
#include <iostream>
#include <limits>
#include <climits>
#include <bitset>
#include <cmath>
#include "primitiveroots.hpp"
#include "ntt.hpp"
#include "polyonfr.cpp"

typedef tuple<FrVec,FrVec> FrvT_2;
typedef tuple<FrVec,FrVec,FrVec> FrvT_3;

using namespace mcl;
using namespace mcl::bn;
using namespace std;

bool IsPolyZero(FrVec a);

FrVec PolyCondense(FrVec a);


bool IsPolyEqual(FrVec a, FrVec b);


FrVec PolyAdd(FrVec a, FrVec b);

//a-b 
FrVec PolySub(FrVec a, FrVec b);


FrVec PolyMul(const FrVec a, const FrVec b);


//A(x)/B(x)
FrVec PolyLongDiv(FrVec A, FrVec B);

Fr PolyEvaluate(const FrVec &a, const Fr &x);


tuple<FrVec,FrVec> PolyDiv(FrVec A, FrVec B);


tuple<FrVec,FrVec> PolyDiv_na(FrVec A, FrVec B);

FrVec PolyDifferentiate(FrVec A);




/*tuple<FrVec,FrVec,FrVec> xGCD2(FrVec a, FrVec b)
{
    if(b.size() > a.size()){
        tuple<FrVec,FrVec,FrVec> out = xGCD2(b,a);
        cout<<"line 298"<<endl;
        return out;
    }else{
        FrVec s={0};
        FrVec old_s={1};

        FrVec r=b;
        FrVec old_r=a;
        cout<<"line 306"<<endl;

        for(;IsPolyZero(r)==0;){
            FrVec quotient; FrVec remainder;
            tuple<FrVec,FrVec>(quotient, remainder)=PolyDiv(old_r,r);
            old_r=r; r=remainder;
            old_s=s; s=PolySub(old_s, PolyMul(quotient,s));
        }
        cout<<"line 314"<<endl;

        FrVec bezout_t;
        if (IsPolyZero(b)==0){
            bezout_t=PolyLongDiv(PolySub(old_r,PolyMul(old_s,a)),b);
            cout<<"line 319"<<endl;
        }else{
            bezout_t={0};
            cout<<"line 322"<<endl;
        }
        cout<<"line 325"<<endl;
        return tuple<FrVec,FrVec,FrVec>(old_r, old_s, bezout_t);
    }
}*/

tuple<FrVec,FrVec,FrVec> xGCD(FrVec a, FrVec b);


FrVec Polytree(FrVec a); 

vector<vector<FrVec>> SubProductTree(FrVec a);


#endif