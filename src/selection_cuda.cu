#include "sketch/bbmh.h"
#include <fstream>
#include <iostream>
#include <seqan/seq_io.h>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <zlib.h>
#include "criteria_smh.cuh" // <-- Usa el header device que definiste
#include <cuda_runtime.h>

// -------------------------------------------------------------------
// FUNCIONES DE UTILIDAD (
// -------------------------------------------------------------------

std::vector<uint64_t> read_smh(std::string path) {
    gzFile fp(gzopen(path.data(), "rb"));
    if(fp == nullptr) throw std::runtime_error(std::string("Could not open file at '") + path + "' for reading");

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

void load_file_list (std::vector<std::string> & files, std::string & list_file, std::string path = "") {
    std::string line;
    if (list_file.empty ()) {
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

// -------------------------------------------------------------------
// KERNEL CUDA PARA COMPARAR TODOS LOS PARES (solo smh_a)
// -------------------------------------------------------------------

__global__ void comparar_pares_smh(
    const uint64_t* sketches,    // [num_files * sketch_size]
    const double* cardinalidades, // [num_files]
    int num_files,
    int sketch_size,
    int n_rows, int n_bands,
    double threshold,
    int* resultado // [total_pairs]: 1 si pasa, 0 si no
) {
    int total_pairs = num_files * (num_files - 1) / 2;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total_pairs) return;

    // idx -> (i, k)
    int i = num_files - 2 - int(sqrtf(-8 * idx + 4 * num_files * (num_files - 1) - 7) / 2.0f - 0.5f);
    int k = idx + i + 1 - num_files * (num_files - 1) / 2 + (num_files - i) * ((num_files - i) - 1) / 2;

    double e1 = cardinalidades[i];
    double e2 = cardinalidades[k];

    const uint64_t* v1 = sketches + i * sketch_size;
    const uint64_t* v2 = sketches + k * sketch_size;

    resultado[idx] = CB_smh_a(threshold, e1, e2, v1, v2, n_rows, n_bands);
}

// -------------------------------------------------------------------
// MAIN
// -------------------------------------------------------------------

int main(int argc, char *argv[])
{
    // --------- PARÁMETROS Y ARGUMENTOS --------------
    std::vector<std::string> files;
    float z_score = 1.96;
    int order_n = 1;
    std::string list_file = "";
    uint threads = 8;
    uint aux_bytes = 256;
    float threshold = 0.9;
    std::string criterion = "";
    char c;

    while ((c = getopt(argc, argv, "xl:t:a:h:c:")) != -1) {
        switch (c) {
            case 'x':
                std::cout << "Usage: -l -t -a -h -c\n";
                return 0;
            case 'l':
                list_file = std::string(optarg);
                break;
            case 't':
                threads = std::stoi(optarg);
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

    std::vector<std::pair<std::string, double>> card_name(files.size());
    std::map<std::string, std::vector<uint64_t>> name2aux;
    uint m = aux_bytes / 8;

    // ------- LEE Y GUARDA LOS SKETCHES Y CARDINALIDADES ---------
    for (size_t i = 0; i < files.size(); ++i) {
        std::string filename = files.at(i);
        name2aux[filename] = read_smh(filename + ".smh" + std::to_string(m));
    }

    std::vector<double> cardinalidades(files.size());
    for (size_t i = 0; i < files.size(); ++i) {
        // Aquí deberías calcular la cardinalidad real de cada archivo:
        // (por ahora usamos el tamaño del sketch como dummy)
        cardinalidades[i] = static_cast<double>(m);  // <--- Modifica si tienes la real
        card_name.at(i) = std::make_pair(files[i], cardinalidades[i]);
    }

    // ---- ORDENAR POR CARDINALIDAD
    std::sort(card_name.begin(), card_name.end(),
        [](const std::pair<std::string, double> &x, const std::pair<std::string, double> &y) {
            return x.second < y.second;
        });

    // ---- REORDENAR SKETCHES Y CARDINALIDADES EN ORDEN
    std::vector<uint64_t> sketches_flat(files.size() * m);
    std::vector<double> cards_sorted(files.size());
    for (size_t i = 0; i < card_name.size(); ++i) {
        const std::string& fname = card_name[i].first;
        std::copy(name2aux[fname].begin(), name2aux[fname].end(), sketches_flat.begin() + i * m);
        cards_sorted[i] = card_name[i].second;
    }

    // ---- PARÁMETROS DE LSH ----
    int n_rows = 1, n_bands = 1;
    for (uint band = 1; band <= m; band++) {
        if (m % band != 0) continue;
        n_bands = band;
        n_rows = m / n_bands;
        float P_r = 1.0 - pow(1.0 - pow(threshold, (float)m / band), (float)band);
        if (P_r >= 0.95) {
            break;
        }
    }

    // ---- COPIA A GPU ----
    int num_files = files.size();
    int sketch_size = m;
    int total_pairs = num_files * (num_files - 1) / 2;

    uint64_t* d_sketches;
    double* d_cards;
    int* d_result;
    cudaMalloc(&d_sketches, sketches_flat.size() * sizeof(uint64_t));
    cudaMalloc(&d_cards, cards_sorted.size() * sizeof(double));
    cudaMalloc(&d_result, total_pairs * sizeof(int));
    cudaMemcpy(d_sketches, sketches_flat.data(), sketches_flat.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_cards, cards_sorted.data(), cards_sorted.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemset(d_result, 0, total_pairs * sizeof(int));

    // ---- LANZAR KERNEL ----
    int blockSize = 128;
    int gridSize = (total_pairs + blockSize - 1) / blockSize;
    comparar_pares_smh<<<gridSize, blockSize>>>(
        d_sketches, d_cards, num_files, sketch_size, n_rows, n_bands, threshold, d_result
    );
    cudaDeviceSynchronize();

    // ---- RECUPERAR RESULTADOS ----
    std::vector<int> h_result(total_pairs);
    cudaMemcpy(h_result.data(), d_result, total_pairs * sizeof(int), cudaMemcpyDeviceToHost);

    // ---- IMPRIMIR PARES QUE PASAN ----
    for (int idx = 0; idx < total_pairs; ++idx) {
        if (h_result[idx]) {
            // idx -> (i, k)
            int i = num_files - 2 - int(std::sqrt(-8 * idx + 4 * num_files * (num_files - 1) - 7) / 2.0 - 0.5);
            int k = idx + i + 1 - num_files * (num_files - 1) / 2 + (num_files - i) * ((num_files - i) - 1) / 2;
            std::cout << card_name[i].first << " " << card_name[k].first << " PASA smh_a\n";
        }
    }

    // ---- LIBERA ----
    cudaFree(d_sketches);
    cudaFree(d_cards);
    cudaFree(d_result);

    return 0;
}
