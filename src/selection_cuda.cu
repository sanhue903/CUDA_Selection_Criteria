// selection_cuda.cu

#include <cuda_runtime.h>
#include <cstdint>
#include "include/criteria_sketch_cuda.cuh"

// === kernel 1: solo smh_a ===============================================
__global__ void kernel_smh(const uint64_t* sketches,
                           const double* cards,
                           int N, int sketch_size,
                           int n_rows, int n_bands,
                           double tau,
                           int* out)
{
    int total_pairs = N*(N-1)/2;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total_pairs) return;

    // map idx -> (i, k)
    int i = N - 2 - int(sqrtf(-8*idx + 4*N*(N-1)-7)*0.5f - 0.5f);
    int k = idx + i + 1 - N*(N-i)/2 + (N-i)*((N-i)-1)/2;

    const uint64_t* v1 = sketches + i*sketch_size;
    const uint64_t* v2 = sketches + k*sketch_size;

    out[idx] = smh_a(v1, v2, n_rows, n_bands, sketch_size);
}

// === kernel 2: CB + smh_a  ==============================================
__global__ void kernel_CBsmh(const uint64_t* sketches,
                             const double* cards,
                             int N, int sketch_size,
                             int n_rows, int n_bands,
                             double tau,
                             int* out)
{
    int total_pairs = N*(N-1)/2;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total_pairs) return;

    int i = N - 2 - int(sqrtf(-8*idx + 4*N*(N-1)-7)*0.5f - 0.5f);
    int k = idx + i + 1 - N*(N-i)/2 + (N-i)*((N-i)-1)/2;

    double e1 = cards[i];
    double e2 = cards[k];
    if (!CB(tau, e1, e2)) { out[idx] = -1; return; }

    const uint64_t* v1 = sketches + i*sketch_size;
    const uint64_t* v2 = sketches + k*sketch_size;
    out[idx] = smh_a(v1, v2, n_rows, n_bands, sketch_size);
}

// === WRAPPERS ===========================================================

extern "C" void launch_kernel_smh(const uint64_t* d_sketches,
                                  const double* d_cards,
                                  int N, int sketch_size,
                                  int n_rows, int n_bands,
                                  double tau,
                                  int* d_out,
                                  int blockSize)
{
    int total_pairs = N * (N - 1) / 2;
    int gridSize = (total_pairs + blockSize - 1) / blockSize;
    kernel_smh<<<gridSize, blockSize>>>(d_sketches, d_cards, N, sketch_size, n_rows, n_bands, tau, d_out);
    cudaDeviceSynchronize(); // Optionally check errors!
}

extern "C" void launch_kernel_CBsmh(const uint64_t* d_sketches,
                                    const double* d_cards,
                                    int N, int sketch_size,
                                    int n_rows, int n_bands,
                                    double tau,
                                    int* d_out,
                                    int blockSize)
{
    int total_pairs = N * (N - 1) / 2;
    int gridSize = (total_pairs + blockSize - 1) / blockSize;
    kernel_CBsmh<<<gridSize, blockSize>>>(d_sketches, d_cards, N, sketch_size, n_rows, n_bands, tau, d_out);
    cudaDeviceSynchronize(); // Optionally check errors!
}
