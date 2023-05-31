/*this polynomial operation code is originated from https://github.com/accumulators-agg/go-poly/blob/master/fft/poly.go */


#include <mcl/bls12_381.hpp>
#include <mcl/ntt.hpp>
#include <iostream>
#include <limits>
#include <climits>
#include <bitset>
#include <cmath>
#include "primitiveroots.hpp"
#include "ntt.hpp"

typedef tuple<FrVec,FrVec> FrvT_2;
typedef tuple<FrVec,FrVec,FrVec> FrvT_3;

using namespace mcl;
using namespace mcl::bn;
using namespace std;

bool IsPolyZero(FrVec a)
{
    uint32_t d=a.size();
    if (d == 0) {
            throw runtime_error("IsPolyZeroError");
        }

    bool flag;
    flag = true;
    for(uint32_t i=0;i<d && flag == true; i++)
    {
        flag = flag && a[i].isZero();
    }
    return flag;
}

FrVec PolyCondense(FrVec a){

    uint32_t d=a.size();
    if (d == 0) {
            throw runtime_error("IsPolyZeroError");
        }


    uint32_t tmpd=d;
    for(uint32_t i=d-1;i>0;i--)
    {

        if( a[i].isZero() != true){
            break;
        }
        tmpd--;

    }
    FrVec res(a.begin(),a.begin()+tmpd);
    return res;

}


bool IsPolyEqual(FrVec a, FrVec b)
{
    FrVec p = PolyCondense(a);
    FrVec q = PolyCondense(b);

    if(p.size() != q.size())
    {
        return false;
    }

    bool flag=true;

    for(uint32_t i=0;i<p.size() && flag == true;i++)
    {
        flag = flag && p[i]==q[i];
    }
    return flag;

}

FrVec PolyAdd(FrVec a, FrVec b)
{
    if (IsPolyZero(a)){
		return PolyCondense(b);
	}

	if (IsPolyZero(b)){
		return PolyCondense(a);
	}

    uint32_t ad = a.size();
    uint32_t bd = b.size();

    uint32_t n =Max(ad, bd);

    FrVec c;
    c.resize(n,0);

    for(uint32_t i=0;i<n;i++)
    {
        if(i<ad)
        {
            c[i]=c[i]+a[i];
        }
        if(i<bd)
        {
            c[i]=c[i]+b[i];
        }
    }

    c=PolyCondense(c);
    return c;

}

//a-b 
FrVec PolySub(FrVec a, FrVec b)
{
    
	if (IsPolyZero(b)){
		return PolyCondense(a);
	}

    uint32_t ad = a.size();
    uint32_t bd = b.size();

    uint32_t n =Max(ad, bd);

    FrVec c;
    c.resize(n,0);

    for(uint32_t i=0;i<n;i++)
    {
        if(i<ad)
        {
            c[i]=c[i]+a[i];
        }
        if(i<bd)
        {
            c[i]=c[i]-b[i];
        }
    }

    c=PolyCondense(c);
    return c;
}

FrVec PolyMul(FrVec a, FrVec b){

    uint32_t N = a.size()+b.size()-1;
    N=nextPowOf2(N);
    a.resize(N,0);
    b.resize(N,0);

    FrVec res(N);

    Ntt<Fr> ntt;

    ntt.init(N);

    Fr* pt_a=&a[0];
    Fr* pt_b=&b[0];
    Fr* pt_res=&res[0];

    ntt.ntt(pt_a);
    ntt.ntt(pt_b);


    for(uint32_t i=0;i<res.size();i++){
        res[i]=a[i]*b[i];
    }

    ntt.intt(pt_res);

    return PolyCondense(res);

}

/*FrVec PolyMul(const FrVec a, const FrVec b)
{
    FrVec res = poly_mul(a, b);
    res=PolyCondense(res);
    return res;
}*/

//A(x)/B(x)
FrVec PolyLongDiv(FrVec A, FrVec B)
{
    if (IsPolyZero(B)) {
            throw runtime_error("PolyDiv: Cannot divide by zero polynomial.");
        }
    
    FrVec a;
    a.resize(A.size(),0);
    for(uint32_t i=0;i<a.size();i++){

        a[i]=A[i];

    }

    int64_t aPos=a.size()-1;
    int64_t bPos=B.size()-1;
    int64_t diff=aPos-bPos;

    FrVec out;
    out.resize(diff+1,0);
    Fr quot;
    for(;diff >= 0;)
    {   

        quot=a[aPos]/B[bPos];
        out[diff]=quot;
        Fr tmp, tmp2;
        for(int64_t i=bPos;i>=0;i--){
            tmp=quot*B[i];
            tmp2=a[diff+i]-tmp;
            a[diff+i]=tmp2;
        }
        aPos-=1;
        diff-=1;
        
    }
   
    return out;
}

Fr PolyEvaluate(const FrVec &a, const Fr &x)
{
    if (IsPolyZero(a)) {
        throw runtime_error("PolyEvaluateError: Zero polynomial");
    }

    FrVec condensed_a = PolyCondense(a);
    uint32_t d = condensed_a.size();

    Fr result;
    Fr power_of_x;
    result.clear();
    power_of_x = 1;

    for (uint32_t i = 0; i < d; i++) {
        Fr term = condensed_a[i] * power_of_x;
        result += term;
        power_of_x *= x;
    }

    return result;
}

tuple<FrVec,FrVec> PolyDiv(FrVec A, FrVec B)
{
    /*if (IsPolyZero(B)) {
            throw runtime_error("PolyDiv: Cannot divide by zero polynomial.");
        }*/

    if (B.size() > A.size()) {
		throw runtime_error("PolyDiv: Deg(B) should be <= Ded(A)");
	}

    FrVec a;
    a.resize(A.size(),0);
    for(uint32_t i=0;i<a.size();i++){

        a[i]=A[i];

    }

    int64_t aPos=a.size()-1;
    int64_t bPos=B.size()-1;
    int64_t diff=aPos-bPos;

    FrVec out;
    out.resize(diff+1,0);

    for(;diff>=0;){
        Fr& quot=out[diff];
        quot=a[aPos]/B[bPos];
        Fr tmp, tmp2;
        for(int64_t i=bPos;i>=0;i--){
            tmp=quot*B[i];
            tmp2=a[diff+i]-tmp;
            a[diff+i]=tmp2;
        }
        aPos-=1;
        diff-=1;
    }
    out=PolyCondense(out);
    a=PolyCondense(a);

    return tuple<FrVec,FrVec>(out,a);
    
}

tuple<FrVec,FrVec> PolyDiv_na(FrVec A, FrVec B)
{
    FrVec quotient=PolyLongDiv(A,B);
    FrVec tmp=PolyMul(quotient,B);
    FrVec remainder=PolySub(A,tmp);

    return tuple<FrVec,FrVec>(quotient,remainder);
}
FrVec PolyDifferentiate(FrVec A)
{
    uint32_t n=A.size();

    if(n==0){
        throw runtime_error("PolyDifferentiate: Input is empty");
    }

    if(n==1){
        FrVec out;
        out.resize(1,0);
        return out;
    }

    FrVec c;
    c.resize(n,0);
    Fr temp;
    for(uint32_t i=1;i<n;i++){
        temp=i;
        c[i]=A[i]*temp;
    }
    FrVec out(c.begin()+1,c.end());
    return out;
}



tuple<FrVec,FrVec,FrVec> xGCD2(FrVec a, FrVec b)
{ 
    if(b.size()>a.size()){
        FrvT_3 out=xGCD2(b,a);
        FrVec s=get<0>(out);
        FrVec t=get<1>(out);
        FrVec g=get<2>(out);
        return FrvT_3(s,t,g);
    }else{

        FrVec s={0};
        FrVec old_s={1};

        FrVec r=b;
        FrVec old_r=a;

        while(!IsPolyZero(r)){
            
            tuple<FrVec,FrVec> qr = PolyDiv(old_r, r);
            FrVec quotient=get<0>(qr);
            FrVec remainder=get<1>(qr);
            old_r = r;
            r = remainder;
            old_s = s;
            s = PolySub(old_s, PolyMul(quotient, s));
            
        }

        FrVec bezout_t;

        if(!IsPolyZero(b)){
            bezout_t =PolyLongDiv(PolySub(old_r,PolyMul(old_s, a)),b);
        }else{
            bezout_t={0};
        }
        return tuple<FrVec,FrVec,FrVec>(old_r,old_s,bezout_t);
    }
}

tuple<FrVec,FrVec,FrVec> xGCD(FrVec a, FrVec b)
{
    if(a.size()<b.size()){
        
        tuple<FrVec,FrVec,FrVec> out = xGCD(b,a);
        FrVec s=get<0>(out);
        FrVec t=get<1>(out);
        FrVec gcd=get<2>(out);

        return tuple<FrVec,FrVec,FrVec>(s,t,gcd);

    }

    FrVec r=a; FrVec r_n=b;

    FrVec s_o={1}; FrVec s={0};
    FrVec t_o={0}; FrVec t={1};

    FrVec qu;
    FrvT_2 qr;
    FrVec s_n;
    FrVec t_n;

    while(IsPolyZero(r_n)==false){

        qr=PolyDiv(r,r_n);
        qu=get<0>(qr);

        s_n=PolySub(s_o,PolyMul(s,qu));
        t_n=PolySub(t_o,PolyMul(t,qu));

        s_o=s;
        t_o=t;

        s=s_n;
        t=t_n;

        r=r_n;
        r_n=get<1>(qr);

    }

    s=s_o;
    t=t_o;
    
    // s=PolyLongDiv(s_o,r);
    // t=PolyLongDiv(t_o,r);
    // r=PolyLongDiv(r,r);

    s=PolyCondense(s);
    t=PolyCondense(t);
    r=PolyCondense(r);

    return FrvT_3(s,t,r);
    
}

FrVec Polytree_n(FrVec a) {

    uint32_t n = a.size();
    uint32_t aLen =n;
    n=nextPowOf2(n);
    uint8_t l = uint8_t(std::log2(n));
    a.resize(n,0);

    

    vector<FrVec> M(n, FrVec(n));
    for(uint32_t j=0;j<n;j++){
        if(j < aLen){
           M[j].resize(2);
           M[j][0]= (-a[j]);
           M[j][1]= 1; 
        }else{
            M[j].resize(1);
            M[j][0]=1;
        }   
    }

    FrVec x;
    FrVec y;


    int64_t index;
    index=0;
    uint32_t L;
    vector<FrVec> m(L, FrVec(L));
    for(uint8_t i=1;i<=l;i++){
        L = (1<<(l-i));
        m.resize(L,FrVec(L));
        for(uint64_t j=0;j<L;j++){
            x=M[index];
            index++;
            y=M[index];
            index++;
            m[j]=PolyMul(x,y);

        }
        index = 0;
        M.clear();
        M=m;
        m.clear();
        x.clear();
        y.clear();
    }
    return M[0];

}

FrVec Polytree(FrVec a) {

    FrVec tmp;
    FrVec A={1};

    if(a.size()==0){
        A={1};
    }else{for(int i=0;i<a.size();i++){
        tmp = {a[i],1};
        A = PolyMul(A, tmp);
    }}

    return A;

}


vector<vector<FrVec>> SubProductTree(FrVec a)
{
    uint32_t n=a.size();

    if(IsPowerOfTwo(n)==false){
        throw runtime_error("SubProductTree inputs needs to be power of two");
    }

    uint8_t l= uint8_t(log2(n));
    vector<vector<FrVec>> M(l+1,vector<FrVec>());
    M[0] = vector<FrVec>(n, FrVec());

    for(uint32_t j=0;j<n;j++){
        M[0][j] = FrVec(2);
        M[0][j][0] = (-a[j]);
        M[0][j][1] = 1;

    }

    FrVec x;
    FrVec y;
    int64_t index=0;
    for(uint8_t i=1;i<=l;i++){
        uint32_t L(1<<(l-i));
        M[i] = vector<FrVec>(L,FrVec(n,0));
        for(uint32_t j=0;j<L;j++){
            x=M[i-1][index];
            index++;
            y=M[i-1][index];
            index++;
            M[i][j]=PolyMul(x,y);
        }
        index=0;
    }
    return M;
}

class Poly {
public:
    vector<Fr> coeffs;

    Poly() {}
    Poly(vector<Fr> c) : coeffs(c) {}

    Poly& operator=(const std::initializer_list<Fr>& il) {
        this->coeffs = vector<Fr>(il.begin(), il.end());
        return *this;
    }

    Poly operator+(const Poly &rhs) const {
        FrVec res_coeffs = PolyAdd(this->coeffs, rhs.coeffs);
        return Poly(res_coeffs);
    }

    Poly operator-(const Poly &rhs) const {
        FrVec res_coeffs = PolySub(this->coeffs, rhs.coeffs);
        return Poly(res_coeffs);
    }

    Poly operator*(const Poly &rhs) const {
        FrVec res_coeffs = PolyMul(this->coeffs, rhs.coeffs);
        return Poly(res_coeffs);
    }

    tuple<Poly, Poly> operator/(const Poly &rhs) const {
        FrvT_2 res_coeffs = PolyDiv(this->coeffs, rhs.coeffs);
        return make_tuple(Poly(get<0>(res_coeffs)), Poly(get<1>(res_coeffs)));
    }

    Poly operator%(const Poly &rhs) const {
        FrvT_2 res_coeffs = PolyDiv(this->coeffs, rhs.coeffs);
        return Poly(get<1>(res_coeffs));
    }

    bool operator==(const Poly &rhs) const {
        return IsPolyEqual(this->coeffs, rhs.coeffs);
    }

    bool operator!=(const Poly &rhs) const {
        return !IsPolyEqual(this->coeffs, rhs.coeffs);
    }

    void print() {
        for (size_t i = 0; i < coeffs.size(); i++) {
            cout << coeffs[i] << "x^" << i;
            if (i + 1 < coeffs.size()) {
                cout << " + ";
            }
        }
        cout << endl;
    }

};



/*int main(){

    initPairing(mcl::BLS12_381);

    cout<<"================================================="<<endl;

    FrVec aa={1,3,3,1}; FrVec bb={6,5,1};

    

    FrVec s0={1}; FrVec s1={0};
    FrVec t0={0}; FrVec t1={1};

    FrVec r0=PolyAdd(PolyMul(s0,aa),PolyMul(t0,bb));
    put("r0",r0);
    FrVec r1=PolyAdd(PolyMul(s1,aa),PolyMul(t1,bb));
    put("r1",r1);

    r0=aa;
    r1=bb;

    FrvT_2 temp;

    cout<<"++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

    
    temp=PolyDiv(r0,r1);
    FrVec quo2=get<0>(temp);
    
    FrVec s2=PolySub(s0,PolyMul(s1,quo2));
    
    put("s2",s2);
    FrVec t2=PolySub(t0,PolyMul(t1,quo2));
    
    put("t2",t2);
    FrVec r2=PolyAdd(PolyMul(aa,s2),PolyMul(bb,t2));
    FrVec inspection=get<1>(temp);

    
    put("r2",inspection);
    
    

    cout<<"++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    temp=PolyDiv(r1,r2);
    FrVec quo3=get<0>(temp);

    FrVec s3=PolySub(s1,PolyMul(s2,quo3));
    put("s3",s3);

    FrVec t3=PolySub(t1,PolyMul(t2,quo3));
    put("t3",t3);

    FrVec r3=PolyAdd(PolyMul(aa,s3),PolyMul(bb,t3));
    inspection=get<1>(temp);

    
    put("r3",inspection);
    
    


    cout<<"++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

    temp=PolyDiv(r2,r3);
    FrVec quo4=get<0>(temp);

    FrVec s4=PolySub(s2,PolyMul(s3,quo4));
    put("s4",s4);

    FrVec t4=PolySub(t2,PolyMul(t3,quo4));
    put("t4",t4);

    FrVec r4=PolyAdd(PolyMul(aa,s4),PolyMul(bb,t4));
    inspection=get<1>(temp);

    
    put("r4",inspection);
    cout<<"++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

    FrvT_3 out2=xGCD(aa,bb);

    FrVec ss=get<0>(out2);
    FrVec tt=get<1>(out2);
    FrVec gg=get<2>(out2);

    put("s",ss);
    put("t",tt);
    put("gcd",gg);

    inspection=PolyAdd(PolyMul(ss,aa),PolyMul(tt,bb));

    put("inspection",inspection);

    cout<<"++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

    FrVec uyt;

    uyt=Polytree(aa);
    
    
    put("uyt",uyt);
    
    vector<vector<FrVec>> dgg;

    dgg=SubProductTree(aa);

    put("1",dgg[0][0]);
    put("2",dgg[1][0]);
    put("3",dgg[2][0]);

}*/
