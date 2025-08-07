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
#include "src/selection_kernels_wrapper.hpp"




uint64_t canonical_kmer (uint64_t kmer, uint k = 31)
{
    uint64_t reverse = 0;
    uint64_t b_kmer = kmer;

    kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
    kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
    kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
    kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
    kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
    reverse = (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - (k << 1));

    return (b_kmer < reverse) ? b_kmer : reverse;
}

void sketch_file (std::vector<uint64_t> &mh_vector, std::string filename, uint k)
{
    sketch::SuperMinHash<> smh(mh_vector.size()-1);
       seqan::SeqFileIn seqFileIn;
       if (!open(seqFileIn, filename.c_str ()))
       {
               // Try opening as gzipped file
               std::string gz_filename = filename + ".gz";
               if (!open(seqFileIn, gz_filename.c_str ()))
               {
                       std::cerr << "ERROR: Could not open the file " << filename << " or " << gz_filename << ".\n";
                       return;
               }
       }

    seqan::CharString id;
    seqan::IupacString seq;

    while (!atEnd (seqFileIn))
    {
        try {
            seqan::readRecord(id, seq, seqFileIn);
        }
        catch (seqan::ParseError &a) {
            break;
        }

        uint64_t kmer = 0;
        uint bases = 0;
        for (size_t i = 0; i < length(seq); ++i)
        {
            uint8_t two_bit = 0;//(char (seq[i]) >> 1) & 0x03;
            bases++;

            switch (char (seq[i]))
            {
                case 'A': two_bit = 0; break;
                case 'C': two_bit = 1; break;
                case 'G': two_bit = 2; break;
                case 'T': two_bit = 3; break;
                case 'a': two_bit = 0; break;
                case 'c': two_bit = 1; break;
                case 'g': two_bit = 2; break;
                case 't': two_bit = 3; break;
                      // Ignore kmer
                default: two_bit = 0; bases = 0; kmer = 0; break;
            }

            kmer = (kmer << 2) | two_bit;
            kmer = kmer & ((1ULL << (k << 1)) - 1);

            if (bases == k)
            {
                //s->add (XXH3_64bits ((const void *) &kmer, sizeof (uint64_t)));
                smh.addh(canonical_kmer (kmer));
                bases--;
            }
        }
    }
    for(size_t i=0;i<mh_vector.size();i++){
        mh_vector[i] =smh.h_[i];
    }
    close (seqFileIn);
}

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
    while (getline (file, line)) {
        // Remove leading/trailing whitespace and carriage returns
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (!line.empty()) {
            files.push_back(path + line);
        }
    }
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
    const uint k = 31;
    const uint sketch_bits = 14;
    float threshold = 0.9;
    int mh_size = 8;
    int block_size = 256;
    int total_rep = 1;

    char c;
    while ((c = getopt(argc, argv, "xl:h:m:b:")) != -1) {
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
            case 'b':
                block_size = std::stoi(optarg);
                break;
            default:
                break;
        }
    }

    // Step 1: Load files
    omp_set_num_threads(threads);
    load_file_list(files, list_file);


    TIMERSTART(construccion)
    std::vector<std::pair<std::string, double>> card_name (files.size ());
    std::map<std::string, std::shared_ptr<sketch::hll_t>> name2hll;
    std::map<std::string, std::vector<uint64_t>> name2mhv;
    //int add_mem = 8*mh_size;

    for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
    {
        std::string filename = files.at (i_processed);
        std::vector<uint64_t> v(mh_size);
        name2hll[filename] = std::make_shared<sketch::hll_t> (sketch_bits);
        name2mhv[filename] = v;
    }


    #pragma omp parallel for schedule(dynamic)
    for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
    {
        std::string filename = files.at (i_processed);
        // Cargar .hll de tamaño 2^14 ya antes creados
        name2hll[filename] = std::make_shared<sketch::hll_t>(filename + ".hll");
        // Cosntruir nuevos sketches auxiliares (hll y mh)
        sketch_file (name2mhv[filename], filename, k);

        auto c = name2hll[filename]->report ();
        card_name.at (i_processed) = std::make_pair (filename, c);
    }

    std::cout << list_file << ";build_smh;" << threshold << ";";
    TIMERSTOP(construccion);
    std::cout << std::endl;


    // Sort by cardinality
    std::sort(card_name.begin(), card_name.end(),
        [](const std::pair<std::string, double> &x, const std::pair<std::string, double> &y) {
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

    // Flatten data for GPU
    const size_t m_hll = 1 << 14;
    const size_t m_aux = mh_size*8;

    int N = files.size();
    std::vector<uint8_t> hll_flat(N * m_hll);
    std::vector<uint64_t> aux_smh_flat(N * mh_size);
    std::vector<double> cards_sorted(N);

    for (int i = 0; i < N; ++i) {
        const std::string& fname = card_name[i].first;
        std::memcpy(aux_smh_flat.data() + i * m_aux, name2mhv[fname].data(), m_aux);
        std::memcpy(hll_flat.data() + i *m_hll, name2hll[fname]->data(), m_hll);
        cards_sorted[i] = card_name[i].second;
    }

   
    std::vector<int2> pairs;
    for (int i = 0; i < N - 1; ++i)
        for (int k = i + 1; k < N; ++k)
            pairs.push_back({   i, k});
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

    for (int rep = 0; rep < total_rep; ++rep) {
        // ---- SMH CUDA
        std::cout << list_file << ";smh_a;" << threshold << ";";
        TIMERSTART(criterio_smh_cuda)
        launch_kernel_smh(d_main, d_aux, d_cd, d_pairs, total_pairs, threshold,
                        mh_size, m_hll, n_rows, 
                        n_bands, d_out, d_out_count, block_size);
        TIMERSTOP(criterio_smh_cuda);
        std::cout << std::endl;

        
        // Get matches count but don't print
        int h_out_count = 0;
        cudaMemcpy(&h_out_count, d_out_count, sizeof(int), cudaMemcpyDeviceToHost);
        std::vector<Result> result(h_out_count);
        cudaMemcpy(result.data(), d_out, h_out_count * sizeof(Result), cudaMemcpyDeviceToHost);

        // ---- CB+SMH CUDA
        std::cout << list_file << ";CB+smh_a;" << threshold << ";";
        TIMERSTART(criterio_CBsmh_cuda)
        launch_kernel_CBsmh(d_main, d_aux, d_cd, d_pairs, total_pairs, threshold,   
                        mh_size, m_hll, n_rows, 
                        n_bands, d_out, d_out_count, block_size);
        TIMERSTOP(criterio_CBsmh_cuda);
        std::cout << std::endl;


        // Get matches count but don't print
        cudaMemcpy(&h_out_count, d_out_count, sizeof(int), cudaMemcpyDeviceToHost);
        result.resize(h_out_count);
        cudaMemcpy(result.data(), d_out, h_out_count * sizeof(Result), cudaMemcpyDeviceToHost);
    }

    cudaFree(d_aux); cudaFree(d_cd); cudaFree(d_out); cudaFree(d_main); cudaFree(d_pairs); cudaFree(d_out_count);
    return 0;
}
