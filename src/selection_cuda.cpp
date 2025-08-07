// selection_main.cpp

#include "sketch/sketch.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>
#include <memory>
#include <getopt.h>
#include <cuda_runtime.h>
#include "include/criteria_sketch.hpp"
#include "include/criteria_sketch_cuda.cuh"
#include "selection_kernels_wrapper.hpp"

// Helper to read SMH sketch from .smhN file
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

void load_file_list(std::vector<std::string>& files, std::string& list_file, std::string path = "") {
    std::string line;
    if (list_file.empty()) {
        std::cerr << "No input file provided\n";
        exit(-1);
    }
    std::ifstream file(list_file);
    if (!file.is_open()) {
        std::cerr << "No valid input file provided\n";
        exit(-1);
    }
    while (getline(file, line)) {
        // Remove leading/trailing whitespace and carriage returns
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (!line.empty()) {
            files.push_back(path + line);
        }
    }
    file.close();
}

int main(int argc, char *argv[]) {
    std::vector<std::string> files;
    std::string list_file = "";
    float threshold = 0.9;
    int aux_bytes = 256;  // will use as m = aux_bytes/8
    std::string criterion = "smh_a";
    int block_size = 256;

    char c;
    while ((c = getopt(argc, argv, "xl:b:a:h:c:")) != -1) {
        switch (c) {
            case 'x':
                std::cout << "Usage: -l -t -a -h -c\n";
                return 0;
            case 'l':
                list_file = std::string(optarg);
                break;
            case 'b':
                block_size = std::stoi(optarg);
                break;
            case 'a':
                aux_bytes = std::stoi(optarg);
                break;
            case 'h':
                threshold = std::stof(optarg);
                break;
            default:
                break;
        }
    }

    load_file_list(files, list_file);

    int m = aux_bytes / 8;
    int N = files.size();

    // Read sketches and cardinalities, build mappings
    std::vector<std::pair<std::string, double>> card_name(N);
    std::map<std::string, std::shared_ptr<sketch::hll_t>> name2hll14;
    std::map<std::string, std::vector<uint64_t>> name2aux;
    for (int i = 0; i < N; ++i) {
        std::string filename = files.at(i);
        name2hll14[filename] = std::make_shared<sketch::hll_t>(filename + ".hll");
        name2aux[filename] = read_smh(filename + ".smh" + std::to_string(m));
    }

    // Calculate cardinalities
    for (int i = 0; i < N; ++i) {
        auto c = name2hll14[files[i]]->report();
        card_name.at(i) = std::make_pair(files[i], c);
    }

    // Sort by cardinality
    std::sort(card_name.begin(), card_name.end(),
              [](const std::pair<std::string, double>& x,
                 const std::pair<std::string, double>& y) {
                  return x.second < y.second;
              });

    upload_pow2neg();
    // Compute n_rows and n_bands
    int n_rows = 1, n_bands = 1;
    for (int band = 1; band <= m; ++band) {
        if (m % band != 0) continue;
        float P_r = 1.0 - pow(1.0 - pow(threshold, (float)m / band), (float)band);
        if (P_r >= 0.95) {
            n_bands = band;
            n_rows = m / band;
            break;
        }
    }

     // Flatten data for GPU
    const size_t m_hll = 1 << 14;
    const size_t buckets_aux = aux_bytes/8;

    std::vector<uint8_t> hll_flat(N * m_hll);
    std::vector<uint64_t> aux_smh_flat(N * buckets_aux);
    std::vector<double> cards_sorted(N);

    for (int i = 0; i < N; ++i) {
        const std::string& fname = card_name[i].first;
        std::memcpy(aux_smh_flat.data() + i * m, name2aux[fname].data(), aux_bytes);
        std::memcpy(hll_flat.data() + i *m_hll, name2hll14[fname]->data(), m_hll);
        cards_sorted[i] = card_name[i].second;
    }

   
    std::vector<int2> pairs;
    for (int i = 0; i < N - 1; ++i)
        for (int k = i + 1; k < N; ++k)
            pairs.push_back({i, k});
    size_t total_pairs = pairs.size();

    // --- Allocate device arrays
    uint8_t* d_main = nullptr;
    uint64_t* d_aux = nullptr;
    double* d_cd = nullptr;
    int2* d_pairs = nullptr;
    Result* d_out = nullptr;
    int* d_out_count = nullptr;

    cudaMalloc(&d_main, hll_flat.size() * sizeof(uint8_t));
    cudaMalloc(&d_aux, aux_smh_flat.size() * sizeof(uint64_t));
    cudaMalloc(&d_cd, cards_sorted.size() * sizeof(double));
    cudaMalloc(&d_pairs, total_pairs * sizeof(int2));
    cudaMalloc(&d_out, total_pairs * sizeof(Result));
    cudaMalloc(&d_out_count, sizeof(int));

    cudaMemcpy(d_main, hll_flat.data(), hll_flat.size() * sizeof(uint8_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_aux, aux_smh_flat.data(), aux_smh_flat.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_cd, cards_sorted.data(), cards_sorted.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_pairs, pairs.data(), total_pairs * sizeof(int2), cudaMemcpyHostToDevice);

    launch_kernel_CBsmh(d_main, d_aux, d_cd, d_pairs, total_pairs, threshold,
                    aux_bytes, m_hll, n_rows, 
                    n_bands, d_out, d_out_count, block_size);

    // Get number of valid results
    int h_out_count = 0;
    cudaMemcpy(&h_out_count, d_out_count, sizeof(int), cudaMemcpyDeviceToHost);
    std::vector<Result> result(h_out_count);
    cudaMemcpy(result.data(), d_out, h_out_count * sizeof(Result), cudaMemcpyDeviceToHost);

    cudaFree(d_aux); cudaFree(d_cd); cudaFree(d_out); cudaFree(d_main); cudaFree(d_pairs); cudaFree(d_out_count);

    for (const auto &pair : result){
        std::cout << card_name[pair.x].first << " " << card_name[pair.y].first << " " << pair.sim << "\n";
    }

    return 0;
}
