#include "criteria_sketch_cuda.cuh"
#include <cuda_runtime.h>
#include "selection_cuda_wrapper.hpp"

// CUDA kernel for comparing all pairs (only smh_a)
__global__ void comparar_pares_smh(
    const uint64_t* sketches,    // [num_files * sketch_size]
    const double* cardinalidades, // [num_files]
    int num_files,
    int sketch_size,
    int n_rows, int n_bands,
    double threshold,
    int* resultado // [total_pairs]: 1 if passes, 0 if not
) {
    int total_pairs = num_files * (num_files - 1) / 2;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total_pairs) return;

    int i = num_files - 2 - int(sqrtf(-8 * idx + 4 * num_files * (num_files - 1) - 7) / 2.0f - 0.5f);
    int k = idx + i + 1 - num_files * (num_files - 1) / 2 + (num_files - i) * ((num_files - i) - 1) / 2;

    double e1 = cardinalidades[i];
    double e2 = cardinalidades[k];

    const uint64_t* v1 = sketches + i * sketch_size;
    const uint64_t* v2 = sketches + k * sketch_size;

    resultado[idx] = CB_smh_a(threshold, e1, e2, v1, v2, n_rows, n_bands);
}

// Wrapper function for C++
void launch_comparar_pares_smh(const uint64_t* sketches, const double* cardinalidades, int num_files, int sketch_size, int n_rows, int n_bands, double threshold, int* resultado) {
    int total_pairs = num_files * (num_files - 1) / 2;
    uint64_t* d_sketches;
    double* d_cards;
    int* d_result;
    cudaMalloc(&d_sketches, num_files * sketch_size * sizeof(uint64_t));
    cudaMalloc(&d_cards, num_files * sizeof(double));
    cudaMalloc(&d_result, total_pairs * sizeof(int));
    cudaMemcpy(d_sketches, sketches, num_files * sketch_size * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_cards, cardinalidades, num_files * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemset(d_result, 0, total_pairs * sizeof(int));

    int blockSize = 128;
    int gridSize = (total_pairs + blockSize - 1) / blockSize;
    comparar_pares_smh<<<gridSize, blockSize>>>(d_sketches, d_cards, num_files, sketch_size, n_rows, n_bands, threshold, d_result);
    cudaDeviceSynchronize();

    cudaMemcpy(resultado, d_result, total_pairs * sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(d_sketches);
    cudaFree(d_cards);
    cudaFree(d_result);
}
