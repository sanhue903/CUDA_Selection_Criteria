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
#include "selection_cuda_wrapper.hpp"


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

    char c;
    while ((c = getopt(argc, argv, "xl:t:a:h:c:")) != -1) {
        switch (c) {
            case 'x':
                std::cout << "Usage: -l -t -a -h -c\n";
                return 0;
            case 'l':
                list_file = std::string(optarg);
                break;
            case 'a':
                aux_bytes = std::stoi(optarg);
                break;
            case 'h':
                threshold = std::stof(optarg);
                break;
            case 'c':
                criterion = std::string(optarg);
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

    // Flatten sketches and cards to arrays
    std::vector<uint64_t> sketches_flat(N * m);
    std::vector<double> cards_sorted(N);
    for (int i = 0; i < N; ++i) {
        std::copy(name2aux[card_name[i].first].begin(),
                  name2aux[card_name[i].first].end(),
                  sketches_flat.begin() + i * m);
        cards_sorted[i] = card_name[i].second;
    }

    // --- Allocate/copy on GPU
    uint64_t* d_sketches;
    double* d_cards;
    int* d_out1;
    int* d_out2;
    int total_pairs = N * (N - 1) / 2;

    cudaMalloc(&d_sketches,  sketches_flat.size() * sizeof(uint64_t));
    cudaMalloc(&d_cards,  cards_sorted.size() * sizeof(double));
    cudaMalloc(&d_out1, total_pairs * sizeof(int));
    cudaMalloc(&d_out2, total_pairs * sizeof(int));

    cudaMemcpy(d_sketches, sketches_flat.data(), sketches_flat.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_cards, cards_sorted.data(), cards_sorted.size() * sizeof(double), cudaMemcpyHostToDevice);

    int block = 256;
    int grid = (total_pairs + block - 1) / block;

    // --- Launch CUDA kernels via wrappers!
    launch_kernel_smh(d_sketches, d_cards, N, m, n_rows, n_bands, threshold, d_out1, block, grid);
    launch_kernel_CBsmh(d_sketches, d_cards, N, m, n_rows, n_bands, threshold, d_out2, block, grid);

    // --- Copy results back
    std::vector<int> smh_result(total_pairs), cbsmh_result(total_pairs);
    cudaMemcpy(smh_result.data(), d_out1, total_pairs * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(cbsmh_result.data(), d_out2, total_pairs * sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(d_sketches); cudaFree(d_cards); cudaFree(d_out1); cudaFree(d_out2);

    // --- Output results with final Jaccard computation (CPU)
    for (int idx = 0, i = 0; i < N - 1; ++i) {
        std::string fn1 = card_name[i].first;
        double e1 = card_name[i].second;
        for (int k = i + 1; k < N; ++k, ++idx) {
            std::string fn2 = card_name[k].first;
            double e2 = card_name[k].second;
            if (cbsmh_result[idx] == 1) {
                double t = name2hll14[fn1]->union_size(*name2hll14[fn2]);
                double jacc14 = (e1 + e2 - t) / t;
                if (jacc14 >= threshold) {
                    std::cout << fn1 << " " << fn2 << " " << jacc14 << "\n";
                }
            }
        }
    }

    return 0;
}
