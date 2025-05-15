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

void sketch_file (std::shared_ptr<sketch::hll_t> hll_aux, std::string filename, int k)
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
				hll_aux->addh(canonical_kmer (kmer));
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
	std::vector<std::string> files1;
	std::vector<std::string> files2;
	std::vector<std::string> files_u;

	std::string list_file1 = "";
	std::string list_file2 = "";
	std::string list_file_u = "";
	uint threads = 8;
	const uint k = 31;
	const uint sketch_bits = 14;
	uint hll_aux_bits = 7;
	float threshold = 0.9;
	float z_score = 1.96;
	int order_n = 1;

	char c;

	while ((c = getopt(argc, argv, "xu:l:L:t:p:h:z:n:")) != -1)
	{
		switch (c) {
			case 'x':
				std::cout << "-u -l -L -t -p -h -z -n -d\n";
				return 0;
			break;
			case 'u':
				list_file_u = std::string (optarg);
			break;			
			case 'l':
				list_file1 = std::string (optarg);
			break;
			case 'L':
				list_file2 = std::string (optarg);
			break;
			case 't':
				threads = std::stoi (optarg);
			break;
			case 'p':
				hll_aux_bits = std::stoi (optarg);
			break;
			case 'h':
				threshold = std::stof (optarg);
			break;
			case 'n':
				order_n = std::stoi (optarg);
			break;
			case 'z':
				z_score = std::stof (optarg);
			break;
			default:
			break;
		}
	}

	// Inicializar variables:
	omp_set_num_threads (threads);
  load_file_list (files1, list_file1);
	load_file_list (files2, list_file2);
	load_file_list (files_u, list_file_u);

	// Paso 1: Cargar los Hyperloglog de tamaño 2^14 ################
	std::vector<std::pair<std::string, double>> card_name1 (files1.size ());
	std::vector<std::pair<std::string, double>> card_name2 (files2.size ());
	//std::map<std::string, double> name2card;
	std::map<std::string, std::shared_ptr<sketch::hll_t>> name2hll;
	std::map<std::string, std::shared_ptr<sketch::hll_t>> name2hll_aux;
//	int mh_size = n_bands*n_rows;

	for (size_t i_processed = 0; i_processed < files_u.size (); ++i_processed)
	{
		std::string filename = files_u.at (i_processed);
		name2hll[filename] = std::make_shared<sketch::hll_t> (sketch_bits);
		name2hll_aux[filename] = std::make_shared<sketch::hll_t> (hll_aux_bits);
	}

	// Cargar sketches desde archivos.
	// card_hll: hll's de tamaño p=14
	// sktwo_hll: hll's auxiliares de tamaño p=r

	// card_name: Pares (nombre, estimado) para card_hll

	#pragma omp parallel for schedule(dynamic)
	for (size_t i_processed = 0; i_processed < files_u.size (); ++i_processed)
	{
		std::string filename = files_u.at (i_processed);
		// Cargar .hll de tamaño 2^14 ya antes creados
		name2hll[filename] = std::make_shared<sketch::hll_t>(filename + ".hll");
		// Cosntruir nuevos sketches auxiliares (hll y mh)
		sketch_file (name2hll_aux[filename], filename, k);

		//auto c = name2hll[filename]->report ();
		// CARD NAME PODRIA DAR PROBLEMAS!!!!!!!!!!!!!!!!!!!!
		//name2card[filename] = c;


	}

	for (size_t i_processed = 0; i_processed < files1.size (); ++i_processed)
	{
		std::string fn1 = files1.at (i_processed);
		auto c1 = name2hll[fn1]->report ();
		card_name1.at (i_processed) = std::make_pair(fn1, c1);

		std::string fn2 = files2.at (i_processed);
		auto c2 = name2hll[fn2]->report ();
		card_name2.at (i_processed) = std::make_pair(fn2, c2);

	}



	// Empezamos a realizar las comparaciones

	// Criterio directo
	long int tp_d = 0;
	long int tn_d = 0;
	long int fp_d = 0;
	long int fn_d = 0;

	// Criterio hllp de orden n
	long int tp_hll = 0;
	long int tn_hll = 0;
	long int fp_hll = 0;
	long int fn_hll = 0;

	// Criterio hllp directo
	long int tp_hlld = 0;
	long int tn_hlld = 0;
	long int fp_hlld = 0;
	long int fn_hlld = 0;

	// CB+Criterio hllp de orden n
	long int tp_cbhll = 0;
	long int tn_cbhll = 0;
	long int fp_cbhll = 0;
	long int fn_cbhll = 0;

	// CB+Criterio hllp directo
	long int tp_cbhlld = 0;
	long int tn_cbhlld = 0;
	long int fp_cbhlld = 0;
	long int fn_cbhlld = 0;

	#pragma omp parallel for schedule(dynamic) reduction(+:tp_d, tn_d, fp_d, fn_d, tp_hll, tn_hll, fp_hll, fn_hll, tp_hlld, tn_hlld, fp_hlld, fn_hlld, tp_cbhll, tn_cbhll, fp_cbhll, fn_cbhll, tp_cbhlld, tn_cbhlld, fp_cbhlld, fn_cbhlld)
	for (size_t i_processed = 0; i_processed < card_name1.size (); ++i_processed)
  {
		std::string fn1 = card_name1[i_processed].first;
    size_t e1 = card_name1[i_processed].second;

    std::string fn2 = card_name2[i_processed].first;
    size_t e2 = card_name2[i_processed].second;
    if(e2 == 0)continue;

    // Cardinalidad de la union "real" y estimada
    double t = name2hll[fn1]->union_size(*name2hll[fn2]);

    // Jaccard real y estimado
    double jacc14 = ((double)e1+(double)e2-t) / t;

    // Selección de los criterios
    bool CB_criterion = CB(threshold, e1, e2);
    bool hllan_criterion = hll_an(threshold, e1, e2, name2hll_aux[fn1], name2hll_aux[fn2], hll_aux_bits, z_score, order_n);
    bool hlla_criterion = hll_a(threshold, e1, e2, name2hll_aux[fn1], name2hll_aux[fn2], hll_aux_bits, z_score);
    bool CBhllan_criterion = CB_hll_an(threshold, e1, e2, name2hll_aux[fn1], name2hll_aux[fn2], hll_aux_bits, z_score, order_n);
    bool CBhlla_criterion = CB_hll_a(threshold, e1, e2, name2hll_aux[fn1], name2hll_aux[fn2], hll_aux_bits, z_score);

    // Comparaciones y estudio:
    if (jacc14 >= threshold){
      // TP and FN:

      // CB
      if (CB_criterion){
        tp_d++;
      }else{
        fn_d++;
      }
      // hll_an
      if (hllan_criterion){
        tp_hll++;
      }else{
        fn_hll++;
      }
      // hll_a
      if (hlla_criterion){
        tp_hlld++;
      }else {
        fn_hlld++;
      }
      // CB+hll_an
      if (CBhllan_criterion){
        tp_cbhll++;
      }else{
        fn_cbhll++;
      }
      // CB+hll_a
      if (CBhlla_criterion){
        tp_cbhlld++;
      }else {
        fn_cbhlld++;
      }
    }else {
      // TN and FP:

      // CB
      if (!CB_criterion){
        tn_d++;
      }else{
        fp_d++;
      }
      // hll_an
      if (!hllan_criterion){
        tn_hll++;
      }else{
        fp_hll++;
      }
      // hll_a
      if (!hlla_criterion){
        tn_hlld++;
      }else {
        fp_hlld++;
      }
      // CB+hll_an
      if (!CBhllan_criterion){
        tn_cbhll++;
      }else{
        fp_cbhll++;
      }
      // CB+hll_a
      if (!CBhlla_criterion){
        tn_cbhlld++;
      }else {
        fp_cbhlld++;
      }
    }
	}


	// Archivo;Método; h; TP; TN; FP; FN; parámetros extras
	std::cout << "Balanced" << ";" <<"CB" << ";" << threshold << ";"
		<< std::to_string(tp_d) << ";" << std::to_string(tn_d) << ";" << std::to_string(fp_d) << ";" << std::to_string(fn_d) << ";"
		<< "NA" << "\n";

	std::cout << "Balanced"  << ";" << "hll_an" << ";" << threshold << ";"
		<< std::to_string(tp_hll) << ";" << std::to_string(tn_hll) << ";" << std::to_string(fp_hll) << ";" << std::to_string(fn_hll) << ";"
		<< "z:" << z_score << "_" << "p:" << hll_aux_bits << "_" << "n:" << order_n << "\n";

	std::cout << "Balanced" << ";" << "hll_a" << ";" << threshold << ";"
		<< std::to_string(tp_hlld) << ";" << std::to_string(tn_hlld) << ";" << std::to_string(fp_hlld) << ";" << std::to_string(fn_hlld) << ";"
		<< "z:" << z_score << "_" << "p:" << hll_aux_bits << "\n";

	std::cout << "Balanced" << ";" << "CB+hll_an" << ";" << threshold << ";"
		<< std::to_string(tp_cbhll) << ";" << std::to_string(tn_cbhll) << ";" << std::to_string(fp_cbhll) << ";" << std::to_string(fn_cbhll) << ";"
		<< "z:" << z_score << "_" << "p:" << hll_aux_bits << "_" << "n:" << order_n << "\n";

	std::cout << "Balanced" << ";" << "CB+hll_a" << ";" << threshold << ";"
		<< std::to_string(tp_cbhlld) << ";" << std::to_string(tn_cbhlld) << ";" << std::to_string(fp_cbhlld) << ";" << std::to_string(fn_cbhlld) << ";"
		<< "z:" << z_score << "_" << "p:" << hll_aux_bits << "\n";

  
	return 0;
}
