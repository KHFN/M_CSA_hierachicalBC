#ifndef MAINALGORITHMS_H__
#define MAINALGORITHMS_H__


#include <mcl/bls12_381.hpp>
#include <iostream>
#include <limits>
#include <cmath>
#include <vector>
#include <malloc.h>
#include "primitiveroots.hpp"
#include "ntt.hpp"
#include "polyonfr.hpp"
#include "mainalgorithms.cpp"

using namespace mcl;
using namespace mcl::bn;
using namespace std;

template <typename T>
void putG(vector<T> a);

struct IDset;


struct zkbpacc_setup;

template <typename T>
T* Vec2arr(vector<T> a);

template <typename T>
vector<T> arr2Vec(T* a, uint32_t n);
template <typename T>
vector<T> convG(vector<T> a, vector<T> b);

/*
a in GG1^n, b in GG2
return \prod_{i=0}^{n-1} e(a[i],b[i]);
*/
GT MultiPairing(vector<G1> a, vector<G2> b);

struct IPPproof;

bool IPPverify(IPPproof pi);

struct zkIPPproof;

IPPproof zkIPPprove(vector<G2> gg, GT P, vector<G1> ww);

bool zkIPPverify(IPPproof pi);

struct PCS;
struct PCS_h;
struct zkbpacc_setup;
struct zkbpacc_setup_pcs;

struct ZKMP;

bool ZKMP_verify(ZKMP pi);

struct ZKNMP;

bool ZKNMP_verify(ZKNMP pi);

struct ZKSP;

bool ZKSP_verify(ZKSP pi);

#endif