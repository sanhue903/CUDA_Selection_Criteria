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

__device__ bool CB_smh_a(double tau, double card_A, double card_B,
                         const uint64_t* v1, const uint64_t* v2, uint n_rows, uint n_bands) {
    return CB(tau, card_A, card_B) || smh_a(v1, v2, n_rows, n_bands);
}
#endif

#endif
