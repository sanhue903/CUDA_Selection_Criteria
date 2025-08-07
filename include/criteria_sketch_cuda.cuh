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
                                 const uint8_t* b,
                                 int            m)
{
    double Z = 0.0;
    int    zeros = 0;

    // 1) fusión por max() y acumulación harmónica -------------------------
    for (int j = 0; j < m; ++j) {
        uint8_t r = max(a[j], b[j]);

        // *** NO añadir 2^0 cuando r == 0 ***
        if (r) Z += ldexpf(1.0f, -int(r));

        zeros += (r == 0);
    }

    const double alpha = 0.7213 / (1.0 + 1.079 / m);
    double est = alpha * m * m / Z;

    // 2) small-range correction (linear-counting) --------------------------
    if (est < 2.5 * m && zeros) {
        est = m * log(static_cast<double>(m) / zeros);
    }

    // 3) large-range correction -------------------------------------------
    constexpr double THRESH = double(1ULL << 32);   // 4 294 967 296
    if (est > THRESH) {
        est = -THRESH * log1p(-est / THRESH);
    }
    return est;
}

#endif

#endif
