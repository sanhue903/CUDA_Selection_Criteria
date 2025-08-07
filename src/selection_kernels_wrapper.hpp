#pragma once
#include <cstdint>
#include <cuda_runtime.h>


void upload_pow2neg(cudaStream_t stream = 0);
struct Result {
    int x, y;
    float sim;
};

void launch_kernel_smh(
    const uint8_t* main_sketches,
    const uint64_t* aux_sketches,

    const double* cards,
    const int2* pairs,

    int total_pairs,
    double tau,

    int m_aux, int m_hll,
    int n_rows, int n_bands,

    Result* out,
    int* out_count,
    int blockSize
);

void launch_kernel_CBsmh(
    const uint8_t* main_sketches,
    const uint64_t* aux_sketches,

    const double* cards,
    const int2* pairs,

    int total_pairs,
    double tau,

    int m_aux, int m_hll,
    int n_rows, int n_bands,

    Result* out,
    int* out_count,
    int blockSize
);
