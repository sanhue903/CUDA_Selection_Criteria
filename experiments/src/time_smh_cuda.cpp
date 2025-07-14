// experiments/src/time_smh_cuda.cpp

#include "sketch/sketch.h"
#include <fstream>
#include <omp.h>
#include <iostream>
#include <seqan/seq_io.h>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <getopt.h>
#include <cuda_runtime.h>
#include "include/metrictime2.hpp"
#include "include/criteria_sketch.hpp"
#include "src/selection_cuda_wrapper.hpp"  // This is your wrapper!

// Function to load file list (reuse)
void load_file_list(std::vector<std::string> & files, std::string & list_file, std::string path = "") {
    std::string line;
    if (list_file.empty ()) {
        std::cerr << "No input file provided\n";
        exit (-1);
    }
    std::ifstream file(list_file);
    if (!file.is_open ()) {
        std::cerr << "No valid input file provided\n";
        exit (-1);
    }
    while (getline (file, line)) files.push_back (path + line);
    file.close();
}

// Read SuperMinHash sketch from file (.smhN format)
std::vector<uint64_t> read_smh(std::string path) {
    gzFile fp(gzopen(path.data(), "rb"));
    if(fp == nullptr) throw std::runtime_error("Could not open file at '" + path + "' for reading");

    uint32_t smh_size;
    if(static_cast<uint64_t>(gzread(fp, &smh_size, sizeof(smh_size))) != sizeof(smh_size)) {
        throw std::runtime_error("Error reading from file (header)\n");
    }
    std::vector<uint64_t> smh_vector(smh_size);
    if(static_cast<uint64_t>(gzread(fp, smh_vector.data(), smh_vector.size() * sizeof(smh_vector[0]))) != smh_vector.size() * sizeof(smh_vector[0])) {
        throw std::runtime_error("Error reading from file (vector)\n");
    }
    gzclose(fp);
    return smh_vector;
}

int main(int argc, char *argv[])
{
    std::vector<std::string> files;
    std::string list_file = "";
    uint threads = 8;
    float threshold = 0.9;
    int mh_size = 8;
    int total_rep = 1;

    char c;
    while ((c = getopt(argc, argv, "xl:t:h:m:R:")) != -1) {
        switch (c) {
            case 'x':
                std::cout << "Usage: -l -t -h -m\n";
                return 0;
            case 'l':
                list_file = std::string(optarg);
                break;
            case 't':
                threads = std::stoi(optarg);
                break;
            case 'h':
                threshold = std::stof(optarg);
                break;
            case 'm':
                mh_size = std::stoi(optarg);
                break;
            case 'R':
                total_rep = std::stoi(optarg);
                break;
            default:
                break;
        }
    }

    // Step 1: Load files
    omp_set_num_threads(threads);
    load_file_list(files, list_file);

    std::cout << list_file << ";build_smh;" << threshold << ";";
    TIMERSTART(construccion)

    std::vector<std::pair<std::string, double>> card_name(files.size());
    std::map<std::string, std::shared_ptr<sketch::hll_t>> name2hll;
    std::map<std::string, std::vector<uint64_t>> name2mhv;

    for (size_t i_processed = 0; i_processed < files.size(); ++i_processed) {
        std::string filename = files.at(i_processed);
        name2hll[filename] = std::make_shared<sketch::hll_t>(14);
        name2mhv[filename] = read_smh(filename + ".smh" + std::to_string(mh_size));
    }

    // Collect cardinalities
    for (size_t i = 0; i < files.size(); ++i) {
        auto c = name2hll[files[i]]->report();
        card_name.at(i) = std::make_pair(files[i], c);
    }

    TIMERSTOP(construccion)
    std::cout << ";m:" << mh_size << "\n";

    // Sort by cardinality
    std::sort(card_name.begin(), card_name.end(),
        [](const std::pair<std::string, double> &x, const std::pair<std::string, double> &y) {
            return x.second < y.second;
        });

    // Calculate bands/rows
    int n_rows = 1, n_bands = 1;
    for (int band = 1; band <= mh_size; band++) {
        if (mh_size % band != 0) continue;
        float P_r = 1.0 - pow(1.0 - pow(threshold, (float)mh_size / band), (float)band);
        if (P_r >= 0.95) {
            n_bands = band;
            n_rows = mh_size / band;
            break;
        }
    }

    // Flatten data for GPU
    int N = files.size();
    std::vector<uint64_t> sketches_flat(N * mh_size);
    std::vector<double> cards_sorted(N);

    for (int i = 0; i < N; ++i) {
        const std::string& fname = card_name[i].first;
        std::copy(name2mhv[fname].begin(), name2mhv[fname].end(), sketches_flat.begin() + i * mh_size);
        cards_sorted[i] = card_name[i].second;
    }

    // Allocate device arrays
    uint64_t* d_sk = nullptr;
    double* d_cd = nullptr;
    int* d_out1 = nullptr;
    int* d_out2 = nullptr;
    int total_pairs = N * (N - 1) / 2;
    cudaMalloc(&d_sk, sketches_flat.size() * sizeof(uint64_t));
    cudaMalloc(&d_cd, cards_sorted.size() * sizeof(double));
    cudaMalloc(&d_out1, total_pairs * sizeof(int));
    cudaMalloc(&d_out2, total_pairs * sizeof(int));
    cudaMemcpy(d_sk, sketches_flat.data(), sketches_flat.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_cd, cards_sorted.data(), cards_sorted.size() * sizeof(double), cudaMemcpyHostToDevice);

    int block = 256;
    int grid = (total_pairs + block - 1) / block;

    for (int rep = 0; rep < total_rep; ++rep) {
        // ---- SMH CUDA
        std::cout << list_file << ";smh_a_cuda;" << threshold << ";";
        TIMERSTART(criterio_smh_cuda)
        launch_kernel_smh(d_sk, d_cd, N, mh_size, n_rows, n_bands, threshold, d_out1, block, grid);
        TIMERSTOP(criterio_smh_cuda)
        std::cout << ";r:" << n_rows << "_" << "b:" << n_bands << "\n";

        // ---- CB+SMH CUDA
        std::cout << list_file << ";CB+smh_a_cuda;" << threshold << ";";
        TIMERSTART(criterio_CBsmh_cuda)
        launch_kernel_CBsmh(d_sk, d_cd, N, mh_size, n_rows, n_bands, threshold, d_out2, block, grid);
        TIMERSTOP(criterio_CBsmh_cuda)
        std::cout << ";r:" << n_rows << "_" << "b:" << n_bands << "\n";
    }

    cudaFree(d_sk); cudaFree(d_cd); cudaFree(d_out1); cudaFree(d_out2);
    return 0;
}
