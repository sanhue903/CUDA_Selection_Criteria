#include "sketch/sketch.h"
#include <fstream>
#include <omp.h>
#include <iostream>
#include <seqan/seq_io.h>
#include <sketch/include/metrictime2.hpp>
#include <sketch/include/criteria_sketch.hpp>
#include <cmath>


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

void sketch_file (std::shared_ptr<sketch::SuperMinHash<>> smh, std::string filename, int k)
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
		catch (seqan::ParseError a) {
			break;
		}

		uint64_t kmer = 0;
		size_t bases = 0;
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


// ----------------------------------------------------------
// ------------------------- MAIN ---------------------------
// ----------------------------------------------------------

int main(int argc, char *argv[])
{
  // Paso 0: Cargar parámetros ##############################
	std::vector<std::string> files;

	std::string list_file = "";
	uint threads = 8;
	const uint k = 31;
	const uint sketch_bits = 14;
	float threshold = 0.9;
	int mh_size = 8;

	char c;

	while ((c = getopt(argc, argv, "xl:t:h:m:")) != -1)
	{
		switch (c) {
			case 'x':
				std::cout << "l:t:h:m:\n";
				return 0;
			case 'l':
				list_file = std::string (optarg);
			break;
			case 't':
				threads = std::stoi (optarg);
			break;
			case 'h':
				threshold = std::stof (optarg);
			break;
			case 'm':
				mh_size = std::stoi (optarg);
			break;
			default:
			break;
		}
	}

	// Inicializar variables:
	omp_set_num_threads (threads);
	load_file_list (files, list_file);

	// Paso 1: Cargar los Hyperloglog de tamaño 2^14 ################
	std::map<std::string, std::shared_ptr<sketch::SuperMinHash<>>> name2smh;
	//std::map<std::string, std::vector<uint64_t>> name2mhv;

	for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
	{
		std::string filename = files.at (i_processed);
		//std::vector<uint64_t> v(mh_size);
		name2smh[filename] = std::make_shared<sketch::SuperMinHash<>> (mh_size);
	}

	// Cargar sketches desde archivos.
	// card_hll: hll's de tamaño p=14
	// sktwo_hll: hll's auxiliares de tamaño p=r

	// card_name: Pares (nombre, estimado) para card_hll

	#pragma omp parallel for schedule(dynamic)
	for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
	{
		std::string filename = files.at (i_processed);
		// Cargar .hll de tamaño 2^14 ya antes creados
		// Cosntruir nuevos sketches auxiliares (hll y mh)
		sketch_file (name2smh[filename], filename, k);

		std::string sketch_aux_size = filename + ".smh"+std::to_string(mh_size);

		name2smh[filename]->write_unfinalized (sketch_aux_size);
	}



	return 0;
}
