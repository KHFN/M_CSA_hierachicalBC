#include <mcl/bls12_381.hpp>
#include <iostream>

using namespace mcl::bn;
using namespace std;

bool IsPowerOfTwo(uint64_t v) {
    return (v & (v - 1)) == 0;
}

uint32_t Max(uint32_t x, uint32_t y) {
    return (x < y) ? y : x;
}

uint32_t Min(uint32_t x, uint32_t y) {
    return (x < y) ? x : y;
}
