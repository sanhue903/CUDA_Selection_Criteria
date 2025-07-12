#include "sketch/sketch.h"
#include <fstream>
#include <iostream>
#include <seqan/seq_io.h>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <zlib.h>
#include "criteria_sketch_cuda.cuh"
#include "selection_cuda_wrapper.hpp"

// -------------------------------------------------------------------
// FUNCIONES DE UTILIDAD 
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

    int num_files = files.size();
    int sketch_size = m;
    int total_pairs = num_files * (num_files - 1) / 2;

    // ---- LLAMAR CUDA DESDE WRAPPER ----
    std::vector<int> h_result(total_pairs);
    launch_comparar_pares_smh(sketches_flat.data(), cards_sorted.data(), num_files, sketch_size, n_rows, n_bands, threshold, h_result.data());

    // ---- IMPRIMIR PARES QUE PASAN ----
    for (int idx = 0; idx < total_pairs; ++idx) {
        if (h_result[idx]) {
            int i = num_files - 2 - int(std::sqrt(-8 * idx + 4 * num_files * (num_files - 1) - 7) / 2.0 - 0.5);
            int k = idx + i + 1 - num_files * (num_files - 1) / 2 + (num_files - i) * ((num_files - i) - 1) / 2;
            std::cout << card_name[i].first << " " << card_name[k].first << " PASA smh_a\n";
        }
    }

    return 0;
}
