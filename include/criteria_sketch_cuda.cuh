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

__device__ double hll_union_card(const uint8_t* a,
                           const uint8_t* b)
{
    const int p = 14;
    const int m = 1 << p;
    uint32_t zeros = 0;
    double sum = 0.0;

    for (int j = 0; j < m; ++j) {
        uint8_t r = max(a[j], b[j]);
        if (r == 0) {
            ++zeros;
        } else {
            /*  ldexp(1.0, -r)   ==  2^(-r)  */
            sum += ldexp(1.0, -int(r));
        }
    }

    const double alpha =
        (m == 16) ? 0.673 :
        (m == 32) ? 0.697 :
        (m == 64) ? 0.709 :
        0.7213 / (1 + 1.079 / m);

    double raw = alpha * m * m / (zeros + sum);

    /* 1) low-range correction (linear counting) */
    if (raw < 2.5 * m && zeros) {
        raw = m * log(double(m) / zeros);
    }
    /* 2) high-range correction  (Flajolet et al. §4.3) */
    else if (raw > (1ULL<<32) / 30.0) {
        raw = -(1ULL<<32) * log1p(-raw / (1ULL<<32));
    }
    return raw;
}

#endif

#endif
