#pragma once
#include <cstdint>

void launch_kernel_smh(
    const uint64_t* d_sk, const double* d_cd,
    int N, int m, int n_rows, int n_bands, double threshold, int* d_out,
    int block, int grid);

void launch_kernel_CBsmh(
    const uint64_t* d_sk, const double* d_cd,
    int N, int m, int n_rows, int n_bands, double threshold, int* d_out,
    int block, int grid);
