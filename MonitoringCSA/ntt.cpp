#include <mcl/bls12_381.hpp>
#include <iostream>
#include <limits>
#include <cmath>
#include "primitiveroots.hpp"

using namespace mcl;
using namespace mcl::bn;
using namespace std;

typedef std::vector<Fr> FrVec;
inline void put(const char *msg, const Fr& x)
{
	printf("%s=%s\n", msg, x.getStr(10).c_str());
}

inline void put(const char *msg, const FrVec& x)
{
	printf("%s=", msg);
	for (size_t i = 0; i < x.size(); i++) {
		printf(" %s", x[i].getStr(10).c_str());
	}
	printf("\n");
}

inline uint64_t nextPowOf2(const uint64_t v) {
    if (v == 0) {
        return 1;
    }
    return uint64_t(1) << static_cast<uint64_t>(log2(v-1)+1);
}

struct NTT {
	int log2N;
	int N;
	Fr g;
	Fr q;
	Fr w;
	Fr invN;
	Fr invW;
	NTT(int log2N)
		: log2N(log2N)
		, N(1 << log2N)
	{
		assert(log2N < 32);
		{
			Vint t;
			Vint::invMod(t, N, Fr::getOp().mp);
			invN.setMpz(t);
		}
		q = -1;
		q /= N;
		uint32_t root = 5;
		w.setStr(scale2roots[log2N]);
		Fr::inv(invW, w);
		//put("invN", invN);
		//put("w", w);
	}
	template<class T>
	void _fft(std::vector<T>& out, const std::vector<T>& in, const Fr& g) const
	{
		out.resize(in.size());
		for (int i = 0; i < N; i++) {
			T& v = out[i];
			v.clear();
			Fr t0;
			Fr::pow(t0, g, i);
			Fr t = 1;
			for (int j = 0; j < N; j++) {
				v += in[j] * t;
				t *= t0;
			}
		}
	}
	template<class T>
	void fft(std::vector<T>& out, const std::vector<T>& in) const
	{
		_fft(out, in, w);
	}
	template<class T>
	void ifft(std::vector<T>& out, const std::vector<T>& in) const
	{
		_fft(out, in, invW);
		for (int i = 0; i < N; i++) {
			out[i] *= invN;
		}
	}
};


FrVec poly_mul(const FrVec f, const FrVec g)
{
	uint32_t l=f.size()+g.size()-1;
	l=nextPowOf2(l);
	NTT ntt(log2(l));
	NTT* data;
	data = new NTT(ntt);
	
	FrVec a, f_n = f;
	f_n.resize(ntt.N,0);
	FrVec b, g_n = g;
	g_n.resize(ntt.N,0);

	ntt.fft(a, f_n);
	ntt.fft(b, g_n);

	FrVec c;
	c.resize(ntt.N,0);

	for(uint32_t i=0;i<ntt.N;i++)
	{
		c[i]=a[i]*b[i];
	}

	FrVec d;

	ntt.ifft(d,c);

	FrVec res;
	res.resize(l);
	for(uint32_t i=0;i<l;i++)
	{
		res[i]=d[i];
	}

	delete data;

	return res;
}

/*int main()
	try
{
	initPairing(mcl::BLS12_381);
	NTT ntt(4);
	FrVec xs(ntt.N), ys, zs, as, bs;
	for (int i = 0; i < ntt.N; i++) {
		xs[i] = 1;
	}
	

	put("xs", xs);
	ntt.fft(ys, xs);
	put("ys", ys);
	ntt.ifft(zs, ys);
	put("zs", zs);

	NTT ntt_2(5);
	xs.resize(ntt_2.N,0);
	ntt_2.fft(as,xs);
	
	

	cout<<ntt.N*ntt.invN<<endl;
	cout<<Fr(480)*ntt.invN<<endl;

	xs={1,2};

	NTT ntt2(6);
	xs.resize(ntt2.N,0);

	put("xs",xs);
	ntt2.fft(as,xs);
	ntt2.fft(bs,xs);
	FrVec cs, ds, fs, gs;
	cs.resize(ntt2.N);

	for(uint32_t i=0;i<ntt2.N;i++){

		cs[i]=as[i]*bs[i]*bs[i]*bs[i];

	}

	ntt2.ifft(ds, cs);

	put("ds",ds);

	FrVec cs, ds;

	cs=poly_mul(xs,xs);
	
	put("cs",cs);

	ds=poly_mul(cs,xs);

	put("ds",ds);

	Fr s;
	s.setByCSPRNG();
	Fr t;
	Fr u;

	for(uint32_t i=0;i<cs.size();i++)
	{
		Fr s_tmp=s;
		Fr::pow(s_tmp,s_tmp,i);
		t=s_tmp*cs[i];
	}

	for(uint32_t i=0;i<xs.size();i++)
	{
		Fr s_tmp=s;
		Fr::pow(s_tmp,s_tmp,i);
		u=s_tmp*xs[i];
	}
	Fr invu;
	Fr::inv(invu,u);

	Fr pi=t*invu;

	G1 g1;
	hashAndMapToG1(g1,"g1");
	G2 g2;
	hashAndMapToG2(g2,"g2");
	G1 g11=g1*pi;
	G2 g22=g2*u;
	G1 g222=g1*t;

	GT gt;
	GT gtt;

	pairing(gt,g11,g22);
	pairing(gtt,g222,g2);
	cout<<g1.getStr(16).c_str()<<endl;
	cout<<gt.getStr(16).c_str()<<endl;
	cout<<gtt.getStr(16).c_str()<<endl;
	

} catch (std::exception& e) {
	printf("err %s\n", e.what());
	return 1;
}*/
//export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../mcl/lib
