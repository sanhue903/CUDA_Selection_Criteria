#include "sketch/sketch.h"
#include <fstream>
#include <iostream>
#include <seqan/seq_io.h>
#include <include/metrictime2.hpp>
#include <include/criteria_sketch.hpp>
#include <cmath>
#include <cuda_runtime.h>
#include <omp.h>


// CUDA-style parallel region for pairwise comparisons (host-side, as in previous translations)
void cuda_pairwise_time_hll(
    const std::vector<std::pair<std::string, double>>& card_name,
    const std::map<std::string, std::shared_ptr<sketch::hll_t>>& name2hll,
    const std::map<std::string, std::shared_ptr<sketch::hll_t>>& name2hll_aux,
    uint hll_aux_bits, float threshold, float z_score, int order_n,
    std::string* out, int mode)
{
    int n = card_name.size();
    #pragma omp parallel for schedule(dynamic)
    for (int i_processed = 0; i_processed < n - 1; ++i_processed) {
        std::string out_str;
        std::string fn1 = card_name[i_processed].first;
        size_t e1 = card_name[i_processed].second;
        for (int k = i_processed + 1; k < n; ++k) {
            std::string fn2 = card_name[k].first;
            size_t e2 = card_name[k].second;
            if(e2 == 0) continue;
            bool select = false;
            if (mode == 0) select = hll_a(threshold, e1, e2, name2hll_aux.at(fn1), name2hll_aux.at(fn2), hll_aux_bits, z_score);
            else if (mode == 1) select = hll_an(threshold, e1, e2, name2hll_aux.at(fn1), name2hll_aux.at(fn2), hll_aux_bits, z_score, order_n);
            else if (mode == 2) { if (!CB(threshold, e1, e2)) break; select = hll_a(threshold, e1, e2, name2hll_aux.at(fn1), name2hll_aux.at(fn2), hll_aux_bits, z_score); }
            else if (mode == 3) { if (!CB(threshold, e1, e2)) break; select = hll_an(threshold, e1, e2, name2hll_aux.at(fn1), name2hll_aux.at(fn2), hll_aux_bits, z_score, order_n); }
            if (!select) continue;
            double t = name2hll.at(fn1)->union_size(*name2hll.at(fn2));
            double jacc14 = ((double)e1+(double)e2-t) / t;
            if (jacc14 >= threshold) out_str += fn1 + " " + fn2 + " " + std::to_string(jacc14) + "\n";
        }
        out[i_processed] = out_str;
    }
}

// Main function (adapted for CUDA-style workflow)
int main(int argc, char *argv[])
{
    // ...existing code...
    for (int rep=0; rep<total_rep; rep++)
    {
        // HLLa
        std::cout << list_file << ";hll_a;" << threshold << ";";
        TIMERSTART(criterio_hllp_directo)
        cuda_pairwise_time_hll(card_name, name2hll, name2hll_aux, hll_aux_bits, threshold, z_score, order_n, out, 0);
        TIMERSTOP(criterio_hllp_directo)
        std::cout << ";z:" << z_score << "_" << "p:" << hll_aux_bits << "\n";
        // HLLan
        std::cout << list_file << ";hll_an;" << threshold << ";";
        TIMERSTART(criterio_hllp_ordenn)
        cuda_pairwise_time_hll(card_name, name2hll, name2hll_aux, hll_aux_bits, threshold, z_score, order_n, out, 1);
        TIMERSTOP(criterio_hllp_ordenn)
        std::cout << ";z:" << z_score << "_" << "p:" << hll_aux_bits << "_" << "n:" << order_n << "\n";
        // CB+HLLa
        std::cout << list_file << ";CB+hll_a;" << threshold << ";";
        TIMERSTART(criterio_CBhllp_directo)
        cuda_pairwise_time_hll(card_name, name2hll, name2hll_aux, hll_aux_bits, threshold, z_score, order_n, out, 2);
        TIMERSTOP(criterio_CBhllp_directo)
        std::cout << ";z:" << z_score << "_" << "p:" << hll_aux_bits << "\n";
        // CB+HLLan
        std::cout << list_file << ";CB+hll_an;" << threshold << ";";
        TIMERSTART(criterio_CBhllp_ordenn)
        cuda_pairwise_time_hll(card_name, name2hll, name2hll_aux, hll_aux_bits, threshold, z_score, order_n, out, 3);
        TIMERSTOP(criterio_CBhllp_ordenn)
        std::cout << ";z:" << z_score << "_" << "p:" << hll_aux_bits << "_" << "n:" << order_n << "\n";
    }
    return 0;
}
