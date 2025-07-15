#include <cuda_runtime.h>
#include <cstdint>
#include <cstdio>
#include "include/criteria_sketch_cuda.cuh"

// ======= DEVICE FUNCTIONS =====================================
// (Assumes smh_a and CB are defined in criteria_sketch_cuda.cuh)

// ======= KERNELS ==============================================

// kernel 1: solo smh_a, now uses precomputed pairs
__global__ void kernel_smh(
    const uint64_t* sketches,
    const double* cards,
    const uint2* pairs,   // NEW: precomputed (i, k) table
    int m,
    int n_rows, int n_bands,
    int total_pairs,
    int* out
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total_pairs) return;

    uint2 p = pairs[idx];
    int i = p.x;
    int k = p.y;

    const uint64_t* v1 = sketches + i * m;
    const uint64_t* v2 = sketches + k * m;

    out[idx] = smh_a(v1, v2, n_rows, n_bands);
}

// kernel 2: CB + smh_a, now uses precomputed pairs
__global__ void kernel_CBsmh(
    const uint64_t* sketches,
    const double* cards,
    const uint2* pairs,   // NEW: precomputed (i, k) table
    int m,
    int n_rows, int n_bands,
    double tau,
    int total_pairs,
    int* out
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total_pairs) return;

    uint2 p = pairs[idx];
    int i = p.x;
    int k = p.y;

    double e1 = cards[i];
    double e2 = cards[k];
    const uint64_t* v1 = sketches + i * m;
    const uint64_t* v2 = sketches + k * m;

    out[idx] = CB(tau, e1, e2) || smh_a(v1, v2, n_rows, n_bands);
}

// ====== WRAPPERS WITH STREAM SUPPORT ===========================

void launch_kernel_smh(
    const uint64_t* d_sketches,
    const double* d_cards,
    const uint2* d_pairs,  // device (i, k) pairs
    int m,
    int n_rows, int n_bands,
    int total_pairs,
    int* d_out,
    int blockSize,
    cudaStream_t stream // NEW: pass stream
) {
    int gridSize = (total_pairs + blockSize - 1) / blockSize;
    kernel_smh<<<gridSize, blockSize, 0, stream>>>(
        d_sketches, d_cards, d_pairs, m, n_rows, n_bands, total_pairs, d_out
    );
    // Optionally: check for errors without sync
    #ifndef NDEBUG
    cudaError_t err = cudaPeekAtLastError();
    if (err != cudaSuccess) {
        printf("kernel_smh launch error: %s\n", cudaGetErrorString(err));
    }
    #endif
}

void launch_kernel_CBsmh(
    const uint64_t* d_sketches,
    const double* d_cards,
    const uint2* d_pairs,  // device (i, k) pairs
    int m,
    int n_rows, int n_bands,
    double tau,
    int total_pairs,
    int* d_out,
    int blockSize,
    cudaStream_t stream // NEW: pass stream
) {
    int gridSize = (total_pairs + blockSize - 1) / blockSize;
    kernel_CBsmh<<<gridSize, blockSize, 0, stream>>>(
        d_sketches, d_cards, d_pairs, m, n_rows, n_bands, tau, total_pairs, d_out
    );
    #ifndef NDEBUG
    cudaError_t err = cudaPeekAtLastError();
    if (err != cudaSuccess) {
        printf("kernel_CBsmh launch error: %s\n", cudaGetErrorString(err));
    }
    #endif
}
