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

__device__ double hll_union_card(const uint8_t* a, const uint8_t* b, int m)
{
    double Z = 0.f;
    for (int j = 0; j < m; ++j) {
        uint8_t r = max(a[j], b[j]);
        // Equivalent to pow2neg: 2^-r
        if (Z) Z += ldexpf(1.0f, -int(r));
    }
    if (Z == 0.0) Z = 1e-9; // Clamp to small positive value
    const double alpha = 0.7213f / (1.f + 1.079f / m);
    return alpha * m * m / Z;
}

#endif

#endif
