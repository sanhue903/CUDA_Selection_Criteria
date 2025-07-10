#include "sketch/sketch.h"
#include <fstream>
#include <iostream>
#include <seqan/seq_io.h>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <zlib.h>
#include <cuda_runtime.h>
#include <include/metrictime2.hpp>
#include <include/criteria_smh.cuh>  // Tu header device

// -------------------------------------------
// Utilidad: lectura sketch SMH desde disco
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
    while (getline(file, line)) files.push_back(path + line);
    file.close();
}

// -------------------------------------------
// Kernel para comparar todos los pares (solo smh_a)
__global__ void comparar_pares_smh(
    const uint64_t* sketches,
    const double* cardinalidades,
    int num_files,
    int sketch_size,
    int n_rows, int n_bands,
    double threshold,
    int* resultado // [total_pairs]: 1 si pasa, 0 si no
) {
    int total_pairs = num_files * (num_files - 1) / 2;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total_pairs) return;

    // idx → (i, k)
    int i = num_files - 2 - int(sqrtf(-8 * idx + 4 * num_files * (num_files - 1) - 7) / 2.0f - 0.5f);
    int k = idx + i + 1 - num_files * (num_files - 1) / 2 + (num_files - i) * ((num_files - i) - 1) / 2;

    double e1 = cardinalidades[i];
    double e2 = cardinalidades[k];

    const uint64_t* v1 = sketches + i * sketch_size;
    const uint64_t* v2 = sketches + k * sketch_size;

    resultado[idx] = smh_a(v1, v2, n_rows, n_bands, sketch_size);
}

// -------------------------------------------
// Kernel para comparar todos los pares (CB + smh_a + lógica break)
__global__ void comparar_pares_CBsmh(
    const uint64_t* sketches,
    const double* cardinalidades,
    int num_files,
    int sketch_size,
    int n_rows, int n_bands,
    double threshold,
    int* resultado // [total_pairs]: 1 si pasa, 0 si no
) {
    int total_pairs = num_files * (num_files - 1) / 2;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total_pairs) return;

    // idx → (i, k)
    int i = num_files - 2 - int(sqrtf(-8 * idx + 4 * num_files * (num_files - 1) - 7) / 2.0f - 0.5f);
    int k = idx + i + 1 - num_files * (num_files - i) / 2 + (num_files - i) * ((num_files - i) - 1) / 2;

    double e1 = cardinalidades[i];
    double e2 = cardinalidades[k];

    const uint64_t* v1 = sketches + i * sketch_size;
    const uint64_t* v2 = sketches + k * sketch_size;

    // Lógica break por fila: si CB falla, marca resultado como -1 y “rompe” comparaciones posteriores para ese i
    __shared__ int break_flag[1];
    if (threadIdx.x == 0) break_flag[0] = 0;
    __syncthreads();

    if (break_flag[0]) {
        resultado[idx] = -1; // “anulado por break”
        return;
    }
    if (!CB(threshold, e1, e2)) {
        break_flag[0] = 1; // activa flag break para todos los hilos siguientes del mismo i
        resultado[idx] = -1;
        return;
    }

    int select = smh_a(v1, v2, n_rows, n_bands, sketch_size);
    resultado[idx] = select;
}

// -------------------------------------------
// MAIN
int main(int argc, char *argv[]) {
    std::vector<std::string> files;
    std::string list_file = "";
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

    load_file_list(files, list_file);

    std::cout << list_file << ";build_smh;" << threshold << ";";
    TIMERSTART(construccion)
    std::vector<std::pair<std::string, double>> card_name(files.size());
    std::map<std::string, std::shared_ptr<sketch::hll_t>> name2hll;
    std::map<std::string, std::vector<uint64_t>> name2mhv;

    for (size_t i = 0; i < files.size(); ++i) {
        std::string filename = files.at(i);
        name2mhv[filename] = read_smh(filename + ".smh" + std::to_string(mh_size));
        name2hll[filename] = std::make_shared<sketch::hll_t>(filename + ".hll");
    }

    for (size_t i = 0; i < files.size(); ++i) {
        auto c = name2hll[files[i]]->report();
        card_name.at(i) = std::make_pair(files[i], c);
    }

    TIMERSTOP(construccion)
    std::cout << ";m:" << mh_size << "\n";

    // Ordena
    std::sort(card_name.begin(), card_name.end(),
        [](const std::pair<std::string, double>& x, const std::pair<std::string, double>& y) {
            return x.second < y.second;
        });

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

    // Serializa los sketches en un buffer plano para GPU
    std::vector<uint64_t> sketches_flat(files.size() * mh_size);
    std::vector<double> cards_sorted(files.size());
    for (size_t i = 0; i < card_name.size(); ++i) {
        const std::string& fname = card_name[i].first;
        std::copy(name2mhv[fname].begin(), name2mhv[fname].end(), sketches_flat.begin() + i * mh_size);
        cards_sorted[i] = card_name[i].second;
    }

    int num_files = files.size();
    int sketch_size = mh_size;
    int total_pairs = num_files * (num_files - 1) / 2;

    uint64_t* d_sketches;
    double* d_cards;
    int* d_result_smh;
    int* d_result_CBsmh;
    cudaMalloc(&d_sketches, sketches_flat.size() * sizeof(uint64_t));
    cudaMalloc(&d_cards, cards_sorted.size() * sizeof(double));
    cudaMalloc(&d_result_smh, total_pairs * sizeof(int));
    cudaMalloc(&d_result_CBsmh, total_pairs * sizeof(int));
    cudaMemcpy(d_sketches, sketches_flat.data(), sketches_flat.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_cards, cards_sorted.data(), cards_sorted.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemset(d_result_smh, 0, total_pairs * sizeof(int));
    cudaMemset(d_result_CBsmh, 0, total_pairs * sizeof(int));

    for (int rep = 0; rep < total_rep; rep++) {
        // -------- CRITERIO SMHa (CUDA) --------
        std::cout << list_file << ";smh_a_cuda;" << threshold << ";";
        TIMERSTART(criterio_smh_cuda)

        int blockSize = 128;
        int gridSize = (total_pairs + blockSize - 1) / blockSize;
        comparar_pares_smh<<<gridSize, blockSize>>>(
            d_sketches, d_cards, num_files, sketch_size, n_rows, n_bands, threshold, d_result_smh
        );
        cudaDeviceSynchronize();
        TIMERSTOP(criterio_smh_cuda)
        std::cout << ";r:" << n_rows << "_" << "b:" << n_bands << "\n";

        // -------- CRITERIO CB+SMHa (CUDA) --------
        std::cout << list_file << ";CB+smh_a_cuda;" << threshold << ";";
        TIMERSTART(criterio_CBsmh_cuda)
        comparar_pares_CBsmh<<<gridSize, blockSize>>>(
            d_sketches, d_cards, num_files, sketch_size, n_rows, n_bands, threshold, d_result_CBsmh
        );
        cudaDeviceSynchronize();
        TIMERSTOP(criterio_CBsmh_cuda)
        std::cout << ";r:" << n_rows << "_" << "b:" << n_bands << "\n";
    }

    cudaFree(d_sketches);
    cudaFree(d_cards);
    cudaFree(d_result_smh);
    cudaFree(d_result_CBsmh);

    return 0;
}
