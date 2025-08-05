#include <cuda_runtime.h>
#include <cstdint>
#include <cstdio>
#include <cmath>
#include "include/criteria_sketch_cuda.cuh"

struct Result {
    int x, y;
    float sim;
};



__global__ void kernel_smh(
    const uint8_t* main_sketches,
    const uint64_t* aux_sketches,

    const double* cards,
    const int2* pairs,

    int total_pairs,
    double tau,

    int m_hll, int m_aux,
    int n_rows, int n_bands,

    Result* out,
    int* out_count
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total_pairs) return;

    int2 p = pairs[idx];
    int i = p.x;
    int k = p.y;

    const uint64_t* v1 = aux_sketches + i * m_aux;
    const uint64_t* v2 = aux_sketches + k * m_aux;

    if (!smh_a(v1, v2, n_rows, n_bands))
        return;

    double c1 = cards[i];
    double c2 = cards[k];

    const uint8_t* main1 = main_sketches + i * m_hll;
    const uint8_t* main2 = main_sketches + k * m_hll;

    double union_card = hll_union_card(main1, main2, m_hll);
    double jacc14 = (c1 + c2 - union_card) / union_card;
    if (jacc14 < tau)
        return;

    int out_idx = atomicAdd(out_count, 1);
    out[out_idx] = {i, k, (float)jacc14};
}

// kernel 2: CB + smh_a, now uses precomputed pairs
__global__ void kernel_CBsmh(
    const uint8_t* main_sketches,
    const uint64_t* aux_sketches,

    const double* cards,
    const int2* pairs,   

    int total_pairs,
    double tau,

    int m_aux, int m_hll,
    int n_rows, int n_bands,
    Result* out,
    int* out_count
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total_pairs) return;

    int2 p = pairs[idx];
    int i = p.x;
    int k = p.y;

    double c1 = cards[i];
    double c2 = cards[k];

    if (!CB(tau, c1, c2))
        return;
    
    const uint64_t* aux1 = aux_sketches + i * m_aux;
    const uint64_t* aux2 = aux_sketches + k * m_aux;
   
    if (!smh_a (aux1, aux2, n_rows, n_bands))
        return;

    const uint8_t* main1 = main_sketches + i * m_hll;
    const uint8_t* main2 = main_sketches + k * m_hll;
    
    double union_card = hll_union_card(main1, main2, m_hll); 
    double jacc14 = (c1 + c2 - union_card) / union_card;
    if (jacc14 < tau)
        return;
    
    // Only write valid candidates
    int out_idx = atomicAdd(out_count, 1);
    out[out_idx] = {i, k, (float)jacc14};
}

void upload_pow2neg(cudaStream_t stream = 0){
    float h_pow2neg[64];
    for (int k = 0; k < 64; ++k)
        h_pow2neg[k] = ldexpf(1.0f, -k);
    
    cudaMemcpyToSymbolAsync(d_pow2neg,
                            h_pow2neg,
                            sizeof(h_pow2neg),
                            0,
                            cudaMemcpyHostToDevice,
                            stream);

}

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
) {
    int gridSize = (total_pairs + blockSize - 1) / blockSize;
    cudaMemset(out_count, 0, sizeof(int));
    kernel_smh<<<gridSize, blockSize>>>(
        main_sketches, aux_sketches, cards, pairs, total_pairs, tau, m_hll, m_aux, n_rows, n_bands, out, out_count
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
) {
    int gridSize = (total_pairs + blockSize - 1) / blockSize;
    // Reset output counter
    cudaMemset(out_count, 0, sizeof(int));
    kernel_CBsmh<<<gridSize, blockSize>>>(
        main_sketches, aux_sketches, cards, pairs, total_pairs, tau, m_aux, m_hll, n_rows, n_bands, out, out_count
    );
    #ifndef NDEBUG
    cudaError_t err = cudaPeekAtLastError();
    if (err != cudaSuccess) {
        printf("kernel_CBsmh launch error: %s\n", cudaGetErrorString(err));
    }
    #endif
}
