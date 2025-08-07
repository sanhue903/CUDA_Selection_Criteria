#ifndef CRITERIA_SMH_CUH
#define CRITERIA_SMH_CUH

#ifdef __CUDACC__
#include <cuda_runtime.h>
#endif
#include <stdint.h>


#ifdef __CUDACC__
__device__ bool CB(double tau, double card_A, double card_B) {
    double gamma = card_A / card_B;
    return (gamma >= tau);
}

__device__ bool smh_a(const uint64_t* v1, const uint64_t* v2, uint n_rows, uint n_bands) {
    for(uint band_id = 0; band_id < n_bands; band_id++) {
        bool check = true;
        for(uint j = 0; j < n_rows; j++) {
            if(v1[band_id * n_rows + j] != v2[band_id * n_rows + j]) {
                check = false;
                break;
            }
        }
        if (check) return true;
    }
    return false;
}

extern __constant__ float d_pow2neg[64];


__device__ double hll_union_card(const uint8_t* a,
                           const uint8_t* b)
{
    const int p = 14;
    const int m = 1 << p;
    uint32_t zeros = 0;
    long double sum = 0.0;

    for (int j = 0; j < m; ++j) {
        uint8_t r = std::max(a[j], b[j]);
        if (r == 0) {
            ++zeros;
        } else {
            /*  ldexp(1.0, -r)   ==  2^(-r)  */
            sum += std::ldexp(1.0L, -int(r));
        }
    }

    const long double alpha =
        (m == 16) ? 0.673L :
        (m == 32) ? 0.697L :
        (m == 64) ? 0.709L :
        0.7213L / (1 + 1.079L / m);

    long double raw = alpha * m * m / (zeros + sum);

    /* 1) low-range correction (linear counting) */
    if (raw < 2.5 * m && zeros) {
        raw = m * std::log(static_cast<long double>(m) / zeros);
    }
    /* 2) high-range correction  (Flajolet et al. §4.3) */
    else if (raw > (1ULL<<32) / 30.0) {
        raw = -(1ULL<<32) * std::log1p(-raw / (1ULL<<32));
    }
    return static_cast<double>(raw);
}

#endif

#endif
