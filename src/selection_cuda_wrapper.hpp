#pragma once
#include <cstdint>
#include <cuda_runtime.h>

void launch_kernel_smh(
    const uint64_t* d_sketches,
    const double* d_cards,
    const uint2* d_pairs,
    int m,
    int n_rows, int n_bands,
    int total_pairs,
    int* d_out,
    int blockSize,
    cudaStream_t stream);

void launch_kernel_CBsmh(
    const uint64_t* d_sketches,
    const double* d_cards,
    const uint2* d_pairs,
    int m,
    int n_rows, int n_bands,
    double tau,
    int total_pairs,
    int* d_out,
    int blockSize,
    cudaStream_t stream);
