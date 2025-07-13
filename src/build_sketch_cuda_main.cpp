#include <seqan/seq_io.h>
#include "sketch/sketch.h"
#include <fstream>
#include "build_sketch_cuda_wrapper.hpp"
#include <iostream>
#include <seqan/seq_io.h>
#include <include/metrictime2.hpp>
#include <cmath>
#include <map>
#include <vector>
#include <string>

// Host-side k-mer extraction
void extract_kmers_from_file(const std::string& filename, uint k, std::vector<uint64_t>& kmers) {
	seqan::SeqFileIn seqFileIn;
	if (!open(seqFileIn, filename.c_str())) return;
	seqan::CharString id;
	seqan::IupacString seq;
	while (!atEnd(seqFileIn)) {
		try {
			seqan::readRecord(id, seq, seqFileIn);
		} catch (seqan::ParseError&) {
			break;
		}
		uint64_t kmer = 0;
		uint bases = 0;
		for (size_t i = 0; i < length(seq); ++i) {
			uint8_t two_bit = 0;
			bases++;
			switch (char(seq[i])) {
				case 'A': case 'a': two_bit = 0; break;
				case 'C': case 'c': two_bit = 1; break;
				case 'G': case 'g': two_bit = 2; break;
				case 'T': case 't': two_bit = 3; break;
				default: two_bit = 0; bases = 0; kmer = 0; break;
			}
			kmer = (kmer << 2) | two_bit;
			kmer = kmer & ((1ULL << (k << 1)) - 1);
			if (bases == k) {
				kmers.push_back(kmer);
				bases--;
			}
		}
	}
	close(seqFileIn);
}

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

void load_file_list (std::vector<std::string> & files, std::string & list_file, std::string path = "") {
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

	load_file_list (files, list_file);

	// Use CUDA for parallel sketching
	std::vector<void*> sketches(files.size());
	int sketch_type = 0;
	if (criterion == "smh_a") sketch_type = 1;
	launch_sketch_files_cuda(files, k, sketch_type, sketches);
	// Write results from sketches to disk
	for (size_t i = 0; i < files.size(); ++i) {
		std::string out_path = files[i] + (sketch_type == 0 ? ".hll" : ".smh");
		if (sketch_type == 0) {
			write_hll(std::static_pointer_cast<sketch::hll_t>(sketches[i]), out_path);
		} else {
			write_smh(std::static_pointer_cast<sketch::SuperMinHash<>>(sketches[i]), out_path);
		}
	}

	return 0;
}
