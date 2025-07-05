// selection_cuda.cu

#include <cuda_runtime.h>
#include <vector>
#include <string>
#include <utility>
#include <map>
#include <memory>
#include "criteria_sketch.hpp" 
#include <zlib.h>
#include <stdexcept>

void load_sketches_and_report_cardinality_cuda(const std::vector<std::string>& files, std::vector<std::pair<std::string, double>>& card_name) {
    int n = files.size();
    card_name.resize(n);
    for (int i = 0; i < n; ++i) {
        std::string filename = files[i];
        std::shared_ptr<sketch::hll_t> hll = std::make_shared<sketch::hll_t>(filename + ".hll");
        double c = hll->report();
        card_name[i] = std::make_pair(filename, c);
    }
}

// CUDA kernel for pairwise comparison (example for hll_a and hll_an criteria)
__global__ void pairwise_compare_kernel(const double* card, int n, double threshold, int* out_matrix) {
    int i = blockIdx.x;
    int j = threadIdx.x + i + 1;
    if (i < n - 1 && j < n) {
        double e1 = card[i];
        double e2 = card[j];
        if (e2 == 0) return;
        // Example: simple threshold comparison (replace with real logic)
        if (fabs(e1 - e2) < threshold) {
            out_matrix[i * n + j] = 1; // Mark as selected
        } else {
            out_matrix[i * n + j] = 0;
        }
    }
}

// Host-side function to perform pairwise comparison using CUDA
void pairwise_compare_cuda(const std::vector<std::pair<std::string, double>>& card_name, double threshold, std::vector<std::vector<int>>& out_matrix) {
    int n = card_name.size();
    std::vector<double> card(n);
    for (int i = 0; i < n; ++i) card[i] = card_name[i].second;
    out_matrix.assign(n, std::vector<int>(n, 0));

    // Device memory
    double* d_card;
    int* d_out_matrix;
    cudaMalloc(&d_card, n * sizeof(double));
    cudaMalloc(&d_out_matrix, n * n * sizeof(int));
    cudaMemcpy(d_card, card.data(), n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemset(d_out_matrix, 0, n * n * sizeof(int));

    // Launch kernel: one block per i, threads for j
    int blocks = n - 1;
    int threads = n;
    pairwise_compare_kernel<<<blocks, threads>>>(d_card, n, threshold, d_out_matrix);
    cudaDeviceSynchronize();

    // Copy results back
    std::vector<int> out_flat(n * n);
    cudaMemcpy(out_flat.data(), d_out_matrix, n * n * sizeof(int), cudaMemcpyDeviceToHost);
    cudaFree(d_card);
    cudaFree(d_out_matrix);

    // Fill out_matrix
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            out_matrix[i][j] = out_flat[i * n + j];
}

std::vector<uint64_t> read_smh(std::string path){
    gzFile fp(gzopen(path.data(), "rb"));
    if(fp == nullptr) throw std::runtime_error(std::string("Could not open file at '") + path + "' for reading");
    uint32_t smh_size;
    if(static_cast<uint64_t>(gzread(fp, &smh_size, sizeof(smh_size))) != sizeof(smh_size)) {
        throw sketch::exception::ZlibError(std::string("Error reading from file\n"));
    }
    std::vector<uint64_t> smh_vector(smh_size);
    if(static_cast<uint64_t>(gzread(fp, smh_vector.data(), smh_vector.size() * sizeof(smh_vector[0]))) != smh_vector.size() * sizeof(smh_vector[0])) {
        throw sketch::exception::ZlibError(std::string("Error reading from file\n"));
    }
    gzclose(fp);
    return smh_vector;
}

void load_file_list (std::vector<std::string> & files, std::string & list_file, std::string path = "")
{
    std::string line;
    if (list_file.empty ()) {
        std::cerr << "No input file provided\n";
        exit (-1);
    }
    std::ifstream file (list_file);
    if (!file.is_open ()) {
        std::cerr << "No valid input file provided\n";
        exit (-1);
    }
    while (getline (file, line)) files.push_back (path + line);
    file.close();
}

// ----------------------------------------------------------
// ------------------------- MAIN ---------------------------
// ----------------------------------------------------------

int main(int argc, char *argv[])
{
    std::vector<std::string> files;
    float z_score = 1.96;
    int order_n = 1;
    std::string list_file = "";
    uint threads = 8;
    uint aux_bytes = 256;
    float threshold = 0.9;
    std::string criterion = "";
    char c;
    while ((c = getopt(argc, argv, "xl:t:a:h:c:")) != -1)
    {
        switch (c) {
            case 'x':
                std::cout << "Usage: -l -t -a -h -c\n";
                return 0;
            case 'l':
                list_file = std::string (optarg);
                break;
            case 't':
                threads = std::stoi (optarg);
                break;
            case 'a':
                aux_bytes = std::stoi (optarg);
                break;
            case 'h':
                threshold = std::stof (optarg);
                break;
            case 'c':
                criterion = std::string (optarg);
                break;
            default:
                break;
        }
    }
    load_file_list (files, list_file);
    std::vector<std::pair<std::string, double>> card_name (files.size ());
    std::map<std::string, std::shared_ptr<sketch::hll_t>> name2hll14;
    std::string out[files.size()];
    if (criterion == "hll_a" || criterion == "hll_an"){
        std::map<std::string, std::shared_ptr<sketch::hll_t>> name2aux;
        uint p = __builtin_ctz (aux_bytes);
        for (size_t i_processed = 0; i_processed < files.size (); ++i_processed) {
            std::string filename = files.at (i_processed);
            name2hll14[filename] = std::make_shared<sketch::hll_t> (14);
            name2aux[filename] = std::make_shared<sketch::hll_t> (p);
        }
        for (size_t i_processed = 0; i_processed < files.size (); ++i_processed) {
            std::string filename = files.at (i_processed);
            name2hll14[filename] = std::make_shared<sketch::hll_t>(filename + ".hll");
            name2aux[filename] = std::make_shared<sketch::hll_t>(filename + ".hll_" + std::to_string(p));
            auto c = name2hll14[filename]->report ();
            card_name.at (i_processed) = std::make_pair (filename, c);
        }
        std::sort (card_name.begin (), card_name.end (),
            [](const std::pair<std::string, double> &x,
                const std::pair<std::string, double> &y)
            { return x.second < y.second; });
        // CUDA pairwise comparison
        std::vector<std::vector<int>> out_matrix;
        pairwise_compare_cuda(card_name, threshold, out_matrix);
        for (size_t i_processed = 0; i_processed < card_name.size () - 1; ++i_processed) {
            std::string out_str;
            std::string fn1 = card_name[i_processed].first;
            size_t e1 = card_name[i_processed].second;
            for (size_t k = i_processed + 1; k < card_name.size (); ++k) {
                if (out_matrix[i_processed][k]) {
                    std::string fn2 = card_name[k].first;
                    double t = name2hll14[fn1]->union_size(*name2hll14[fn2]);
                    double jacc14 = ((double)e1+(double)card_name[k].second-t) / t;
                    if (jacc14 >= threshold) out_str += fn1 + " " + fn2 + " " + std::to_string(jacc14) + "\n";
                }
            }
            out[i_processed] = out_str;
        }
    } else if (criterion == "smh_a") {
        std::map<std::string, std::vector<uint64_t>> name2aux;
        uint m = aux_bytes/8;
        for (size_t i_processed = 0; i_processed < files.size (); ++i_processed) {
            std::string filename = files.at (i_processed);
            std::vector<uint64_t> v(m);
            name2hll14[filename] = std::make_shared<sketch::hll_t> (14);
            name2aux[filename] = v;
        }
        for (size_t i_processed = 0; i_processed < files.size (); ++i_processed) {
            std::string filename = files.at (i_processed);
            name2hll14[filename] = std::make_shared<sketch::hll_t>(filename + ".hll");
            name2aux[filename] = read_smh(filename + ".smh"+std::to_string(m));
            auto c = name2hll14[filename]->report ();
            card_name.at (i_processed) = std::make_pair (filename, c);
        }
        std::sort (card_name.begin (), card_name.end (),
            [](const std::pair<std::string, double> &x,
                const std::pair<std::string, double> &y)
            { return x.second < y.second; });
        int n_rows = 1, n_bands = 1;
        for (uint band=1; band <=m; band++){
            if (m % band != 0) continue;
            n_bands = band;
            n_rows = m / n_bands;
            float P_r = 1.0-pow(1.0-pow(threshold, (float)m/band), (float)band);
            if (P_r >= 0.95){
                break;
            }
        }
        // CUDA pairwise comparison (reuse kernel or write new one for SMH if needed)
        std::vector<std::vector<int>> out_matrix;
        pairwise_compare_cuda(card_name, threshold, out_matrix);
        for (size_t i_processed = 0; i_processed < card_name.size () - 1; ++i_processed) {
            std::string out_str;
            std::string fn1 = card_name[i_processed].first;
            size_t e1 = card_name[i_processed].second;
            for (size_t k = i_processed + 1; k < card_name.size (); ++k) {
                if (out_matrix[i_processed][k]) {
                    std::string fn2 = card_name[k].first;
                    double t = name2hll14[fn1]->union_size(*name2hll14[fn2]);
                    double jacc14 = ((double)e1+(double)card_name[k].second-t) / t;
                    if (jacc14 >= threshold) out_str += fn1 + " " + fn2 + " " + std::to_string(jacc14) + "\n";
                }
            }
            out[i_processed] = out_str;
        }
    } else {
        std::cout << "Option -c invalid. The accepted criteria are hll_a, hll_an and smh_a.\n";
    }
    for (size_t i_processed = 0; i_processed < card_name.size (); ++i_processed) {
        std::cout << out[i_processed];
    }
    return 0;
}

