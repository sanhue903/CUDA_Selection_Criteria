// ...existing code...
#include "sketch/sketch.h"
#include <seqan/seq_io.h>
#include <string>
#include <stdint.h>
#include <cuda_runtime.h>
#include "build_sketch_cuda_wrapper.hpp"


__device__ uint64_t canonical_kmer_cuda(uint64_t kmer, uint k) {
    uint64_t reverse = 0;
    uint64_t b_kmer = kmer;
    kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
    kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
    kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
    kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
    kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
    reverse = (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - (k << 1));
    return (b_kmer < reverse) ? b_kmer : reverse;
}

__device__ uint64_t hash_kmer_cuda(uint64_t kmer) {
    kmer ^= (kmer >> 33);
    kmer *= 0xff51afd7ed558ccdULL;
    kmer ^= (kmer >> 33);
    kmer *= 0xc4ceb9fe1a85ec53ULL;
    kmer ^= (kmer >> 33);
    return kmer;
}

__global__ void canonicalize_and_hash_kmers_batch(const uint64_t* kmers, uint64_t* out, uint k, int do_hash, int total_kmers) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total_kmers) return;
    uint64_t canon = canonical_kmer_cuda(kmers[idx], k);
    out[idx] = do_hash ? hash_kmer_cuda(canon) : canon;
}




// Host-side batch implementation
void launch_sketch_files_cuda(const std::vector<std::vector<uint64_t>> &all_kmers, uint k, int sketch_type, std::vector<void*> &sketches) {
    int do_hash = 1; // Always hash after canonicalization
    int total_kmers = 0;
    for (const auto& v : all_kmers) total_kmers += v.size();
    // Flatten all k-mers into one buffer
    std::vector<uint64_t> flat_kmers;
    std::vector<int> offsets(all_kmers.size()+1, 0);
    for (size_t i = 0; i < all_kmers.size(); ++i) {
        offsets[i] = flat_kmers.size();
        flat_kmers.insert(flat_kmers.end(), all_kmers[i].begin(), all_kmers[i].end());
    }
    offsets[all_kmers.size()] = flat_kmers.size();
    // Allocate device memory
    uint64_t *d_kmers, *d_out;
    cudaMalloc(&d_kmers, flat_kmers.size() * sizeof(uint64_t));
    cudaMalloc(&d_out, flat_kmers.size() * sizeof(uint64_t));
    cudaMemcpy(d_kmers, flat_kmers.data(), flat_kmers.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
    // Launch kernel
    int blockSize = 128;
    int gridSize = (flat_kmers.size() + blockSize - 1) / blockSize;
    canonicalize_and_hash_kmers_batch<<<gridSize, blockSize>>>(d_kmers, d_out, k, do_hash, flat_kmers.size());
    cudaDeviceSynchronize();
    // Copy results back
    std::vector<uint64_t> flat_out(flat_kmers.size());
    cudaMemcpy(flat_out.data(), d_out, flat_kmers.size() * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    // Split results back into per-file sketches
    for (size_t i = 0; i < all_kmers.size(); ++i) {
        std::vector<uint64_t> result(flat_out.begin() + offsets[i], flat_out.begin() + offsets[i+1]);
        if (sketch_type == 0) {
            auto hll = std::make_shared<sketch::hll_t>(14); // or your chosen precision
            for (auto val : result) hll->addh(val);
            sketches[i] = hll;
        } else {
            auto smh = std::make_shared<sketch::SuperMinHash<>>(result.size());
            for (auto val : result) smh->addh(val);
            sketches[i] = smh;
        }
    }
    cudaFree(d_kmers);
    cudaFree(d_out);
}
