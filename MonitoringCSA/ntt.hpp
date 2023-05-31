#ifndef NTT_H__
#define NTT_H__

#include <mcl/bls12_381.hpp>
#include <iostream>
#include <limits>
#include <cmath>
#include "ntt.cpp"
#include "primitiveroots.hpp"

using namespace mcl;
using namespace mcl::bn;
using namespace std;

typedef std::vector<Fr> FrVec;

inline void put(const char *msg, const Fr& x);

inline void put(const char *msg, const FrVec& x);

inline uint64_t nextPowOf2(const uint64_t v);

struct NTT; /*{
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
		put("invN", invN);
		put("w", w);
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
};*/

FrVec poly_mul(const FrVec f, const FrVec g);
/*{
	uint32_t l=f.size()+g.size()-1;
	l=nextPowOf2(l);
	NTT ntt(log2(l));
	cout<<l<<""<<log2(l)<<endl;

	
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
	return res;

}*/

#endif