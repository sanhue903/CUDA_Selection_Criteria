#include "sketch/hll.h"
#include <fstream>
#include <omp.h>
#include <iostream>
#include <seqan/seq_io.h>

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

double jaccard_index (uint u_card, uint i_card)
{
	if (i_card > u_card) return 0;
	return i_card / double (u_card);
}

double compute_distance (std::shared_ptr<sketch::hll_t> s_1, std::shared_ptr<sketch::hll_t> s_2, size_t c_1, size_t c_2)
{
	double dist = 0;
	size_t j = s_1->jaccard_index (*s_2);
	size_t inter = (c_1 + c_2) - j;
	dist = jaccard_index (j, inter);
	std::cout<<" dist "<<dist<<std::endl;

	return s_1->jaccard_index (*s_2);
}

void estimate_file (std::shared_ptr<sketch::hll_t> s, std::string filename, int k)
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
				s->addh(canonical_kmer (kmer));
				bases--;
			}
		}
	}
	close (seqFileIn);
}

void load_file_list (std::vector<std::string> & files, std::string & list_file, std::string path = "./")
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

struct data_s {
	std::shared_ptr<sketch::hll_t> sketch;
	std::string filename;
	double card;
};

int main(int argc, char *argv[])
{
	std::vector<std::string> files;

	std::string list_file = "";
	uint threads = 8;
	const uint k = 31;
	const uint sketch_bits = 14;
	float threshold = 0.9;
	std::string filename_matrix;
	std::ofstream o_matrix;

	bool build_sketches = false;

	char c;

	while ((c = getopt(argc, argv, "b:l:c:p:t:o:w")) != -1)
	{
		switch (c) {
			case 'l':
				list_file = std::string (optarg);
			break;
			case 'p':
				threads = std::stoi (optarg);
			break;
			case 't':
				threshold = std::stof (optarg);
			break;
			case 'o':
				filename_matrix = optarg;
			break;
			case 'w':
				build_sketches = true;
			break;
			default:
			break;
		}
	}

	omp_set_num_threads (threads);
	load_file_list (files, list_file);

	if (!filename_matrix.empty ()) o_matrix.open (filename_matrix, std::ios::out);
	std::ostream & output_matrix = (filename_matrix.empty ()) ? std::cout : o_matrix;

	std::vector<std::pair<std::string, double>> card_name (files.size ());
	std::map<std::string, std::shared_ptr<sketch::hll_t>> card_hll;
	std::string out[files.size ()];

	auto es0_time = std::chrono::steady_clock::now ();

	for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
	{
		std::string filename = files.at (i_processed);
		card_hll[filename] = std::make_shared<sketch::hll_t> (sketch_bits);
	}

	std::cout << "Filling sketches for " << files.size () << " files\n";
	// Load sketches from file
	if (build_sketches)
	{
		#pragma omp parallel for schedule(dynamic)
		 for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
		 {
			std::string filename = files.at (i_processed);
			estimate_file (card_hll[filename], filename, k);
			card_hll[filename]->write (filename + ".hll");

			auto c = card_hll[filename]->report ();
			card_name.at (i_processed) = std::make_pair (filename, c);
		}
	}
	else
	{
		#pragma omp parallel for schedule(dynamic)
		for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
		{
			std::string filename = files.at (i_processed);
			card_hll[filename] = std::make_shared<sketch::hll_t>(filename + ".hll");

			auto c = card_hll[filename]->report ();
			card_name.at (i_processed) = std::make_pair (filename, c);
		}
	}

	auto es1_time = std::chrono::steady_clock::now ();
	auto time_hll = std::chrono::duration_cast<std::chrono::microseconds>(es1_time - es0_time).count();

	std::cout << "time_hll " << (time_hll / 1000) << "[ms]\n";

	std::sort (card_name.begin (), card_name.end (),
	           [](const std::pair<std::string, double> &x,
	              const std::pair<std::string, double> &y)
	{
		return x.second < y.second;
	});

	std::cout << "Sorting " << card_name.size () << " cardinalities\n";

	for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
	{
		std::cout << card_name.at (i_processed).first << " " << size_t(card_name.at (i_processed).second) << "\n";
	}
	//return 0;

	auto es2_time = std::chrono::steady_clock::now ();
	auto time_sort = std::chrono::duration_cast<std::chrono::microseconds>(es2_time - es1_time).count();

	std::cout << "time_sort " << (time_sort / 1000) << "[ms]\n";

	std::cout << "Computing distances\n";

	std::cout << "Threshold " << threshold << "\n";

	#pragma omp parallel for schedule(dynamic)
	for (size_t i_processed = 0; i_processed < card_name.size () - 1; ++i_processed)
	{
		std::string fn1 = card_name[i_processed].first;
		size_t e1 = card_name[i_processed].second;
		std::string out_str;

		size_t k = i_processed + 1;
		for (; k < card_name.size (); ++k)
		{
			double fraction = e1 / card_name[k].second;
			std::string fn2 = card_name[k].first;
			if (fraction < threshold) break;
		}

		for (size_t i = i_processed + 1; i < k; ++i)
		{
			std::string fn2 = card_name[i].first;
			size_t e2 = card_name[i].second;
			double dist = card_hll[fn1]->jaccard_index (*card_hll[fn2]);
			//std::cout << fn1 << " " << fn2 << " " << dist << " ads " << e1 << " " << e2 << "\n";

			out_str += fn1 + " " + fn2 + " " + std::to_string (dist) + "\n";
		}

		out[i_processed] = out_str;
	}

	auto es3_time = std::chrono::steady_clock::now ();
	auto time_dist = std::chrono::duration_cast<std::chrono::microseconds>(es3_time - es2_time).count();
	std::cout << "time_dist " << (time_dist / 1000) << "[ms]\n";

	for (size_t i_processed = 0; i_processed < card_name.size (); ++i_processed)
	{
		output_matrix << out[i_processed];
	}


	return 0;
}

