#pragma once
#include <cstdint>

void launch_kernel_smh(
    const uint64_t* d_sketches, const double* d_cards,
    int N, int m, int n_rows, int n_bands, double tau, int* d_out,
    int block, int grid);

void launch_kernel_CBsmh(
    const uint64_t* d_sketches, const double* d_cards,
    int N, int m, int n_rows, int n_bands, double tau, int* d_out,
    int blockSize, int gridSize);
