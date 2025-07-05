#include "sketch/sketch.h"
#include <fstream>
#include <iostream>
#include <seqan/seq_io.h>
#include <include/metrictime2.hpp>
#include <cmath>
#include <cuda_runtime.h>
#include <omp.h>
#include <zlib.h>

void write_smh(std::shared_ptr<sketch::SuperMinHash<>> smh, std::string path){
    std::vector<uint64_t> smh_vector = smh->h_;
    uint32_t smh_size = smh->h_.size();
    gzFile fp(gzopen(path.data(), "wb"));
    if(!fp) throw sketch::exception::ZlibError(Z_ERRNO, std::string("Could not open file at '") + path + "' for writing");
    if(gzwrite(fp, &(smh_size), sizeof(smh_size)) == 0) throw std::runtime_error("Error writing to file.");
    if(gzwrite(fp, smh_vector.data(), smh_vector.size()*sizeof(smh_vector[0])) == 0) throw std::runtime_error("Error writing to file.");
    gzclose(fp);
}

void write_hll(std::shared_ptr<sketch::hll_t> hll, std::string path){
    hll->write(path);
}

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

void sketch_hll (std::shared_ptr<sketch::hll_t> s, std::string filename, uint k)
{
    seqan::SeqFileIn seqFileIn;
    if (!open(seqFileIn, filename.c_str ()))
    {
        std::cerr << "ERROR: Could not open the file " << filename << ".\n";
        return;
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
            uint8_t two_bit = 0;
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
                default: two_bit = 0; bases = 0; kmer = 0; break;
            }
            kmer = (kmer << 2) | two_bit;
            kmer = kmer & ((1ULL << (k << 1)) - 1);
            if (bases == k)
            {
                s->addh(canonical_kmer (kmer));
                bases--;
            }
        }
    }
    close (seqFileIn);
}

void sketch_smh (std::shared_ptr<sketch::SuperMinHash<>> smh, std::string filename, uint k)
{
    seqan::SeqFileIn seqFileIn;
    if (!open(seqFileIn, filename.c_str ()))
    {
        std::cerr << "ERROR: Could not open the file " << filename << ".\n";
        return;
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
            uint8_t two_bit = 0;
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
                default: two_bit = 0; bases = 0; kmer = 0; break;
            }
            kmer = (kmer << 2) | two_bit;
            kmer = kmer & ((1ULL << (k << 1)) - 1);
            if (bases == k)
            {
                smh->addh(canonical_kmer (kmer));
                bases--;
            }
        }
    }
    close (seqFileIn);
}

void load_file_list (std::vector<std::string> & files, std::string & list_file, std::string path = "")
{
    std::string line;
    if (list_file.empty ())
    {
        std::cerr << "No input file provided\n";
        exit (-1);
    }
    std::ifstream file (list_file);
    if (!file.is_open ())
    {
        std::cerr << "No valid input file provided\n";
        exit (-1);
    }
    while (getline (file, line)) files.push_back (path + line);
    file.close();
}

void parallel_sketch_hll_cuda(const std::vector<std::string>& files, std::map<std::string, std::shared_ptr<sketch::hll_t>>& name2hll, int k) {
    int n = files.size();
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; ++i) {
        std::string filename = files[i];
        sketch_hll(name2hll[filename], filename, k);
        write_hll(name2hll[filename], filename + ".hll");
    }
}

void parallel_aux_hll_cuda(const std::vector<std::string>& files, std::map<std::string, std::shared_ptr<sketch::hll_t>>& name2hll, int k, int p) {
    int n = files.size();
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; ++i) {
        std::string filename = files[i];
        sketch_hll(name2hll[filename], filename, k);
        write_hll(name2hll[filename], filename + ".hll_" + std::to_string(p));
    }
}

void parallel_smh_cuda(const std::vector<std::string>& files, std::map<std::string, std::shared_ptr<sketch::SuperMinHash<>>>& name2smh, int k, int m) {
    int n = files.size();
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; ++i) {
        std::string filename = files[i];
        sketch_smh(name2smh[filename], filename, k);
        write_smh(name2smh[filename], filename + ".smh" + std::to_string(m));
    }
}

int main(int argc, char *argv[])
{
    std::vector<std::string> files;
    const uint k = 31;
    std::string list_file = "";
    uint threads = 8;
    uint aux_bytes = 256;
    std::string criterion = "";
    char c;
    while ((c = getopt(argc, argv, "l:t:a:c:")) != -1)
    {
        switch (c) {
            case 'l':
                list_file = std::string (optarg);
                break;
            case 't':
                threads = std::stoi (optarg);
                break;
            case 'a':
                aux_bytes = std::stoi (optarg);
                break;
            case 'c':
                criterion = std::string (optarg);
                break;
            default:
                break;
        }
    }
    omp_set_num_threads (threads);
    load_file_list (files, list_file);
    std::map<std::string, std::shared_ptr<sketch::hll_t>> name2hll14;
    for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
    {
        std::string filename = files.at (i_processed);
        name2hll14[filename] = std::make_shared<sketch::hll_t> (14);
    }
    parallel_sketch_hll_cuda(files, name2hll14, k);
    if (criterion == "hll_a" || criterion == "hll_an"){
        uint p = __builtin_ctz (aux_bytes);
        std::map<std::string, std::shared_ptr<sketch::hll_t>> name2hll;
        for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
        {
            std::string filename = files.at (i_processed);
            name2hll[filename] = std::make_shared<sketch::hll_t> (p);
        }
        parallel_aux_hll_cuda(files, name2hll, k, p);
    }else if (criterion == "smh_a"){
        uint m = aux_bytes/8;
        std::map<std::string, std::shared_ptr<sketch::SuperMinHash<>>> name2smh;
        for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
        {
            std::string filename = files.at (i_processed);
            name2smh[filename] = std::make_shared<sketch::SuperMinHash<>> (m);
        }
        parallel_smh_cuda(files, name2smh, k, m);
    }else{
        printf("Option -c invalid. The accepted criteria are hll_a, hll_an and smh_a.\n");
    }
    return 0;
}
