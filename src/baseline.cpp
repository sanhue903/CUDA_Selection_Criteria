#include "sketch/sketch.h"
#include <fstream>
#include <omp.h>
#include <iostream>
#include <seqan/seq_io.h>


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

	char c;

	while ((c = getopt(argc, argv, "xl:t:")) != -1)
	{
		switch (c) {
			case 'x':
				std::cout << "l:t:\n";
				return 0;
			case 'l':
				list_file = std::string (optarg);
			break;
			case 't':
				threads = std::stoi (optarg);
			break;
			default:
			break;
		}
	}

	// Inicializar variables:
	omp_set_num_threads (threads);
	load_file_list (files, list_file);
	std::string out[files.size()];

	// Paso 1: Cargar los Hyperloglog de tamaño 2^14 ################
	std::vector<std::pair<std::string, double>> card_name (files.size ());
	std::map<std::string, std::shared_ptr<sketch::hll_t>> name2hll;

	for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
	{
		std::string filename = files.at (i_processed);
		name2hll[filename] = std::make_shared<sketch::hll_t> (sketch_bits);
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
		name2hll[filename] = std::make_shared<sketch::hll_t>(filename + ".hll");
		// Cosntruir nuevos sketches auxiliares (hll y mh)

		auto c = name2hll[filename]->report ();
		card_name.at (i_processed) = std::make_pair (filename, c);
	}

	// Ordenamos los pares de card_name según la estimación (hll's de tamaño p=14).
	std::sort (card_name.begin (), card_name.end (),
	           [](const std::pair<std::string, double> &x,
	              const std::pair<std::string, double> &y)
	{
		return x.second < y.second;
	});


	// Empezamos a realizar las comparaciones




  #pragma omp parallel for schedule(dynamic)
  for (size_t i_processed = 0; i_processed < card_name.size () - 1; ++i_processed){
    std::string out_str;

    std::string fn1 = card_name[i_processed].first;
    size_t e1 = card_name[i_processed].second;

    size_t k = i_processed + 1;
    for (; k < card_name.size (); ++k)
    {
      std::string fn2 = card_name[k].first;
      size_t e2 = card_name[k].second;

      // Cardinalidad de la union "real" y estimada
      double t = name2hll[fn1]->union_size(*name2hll[fn2]);

      // Jaccard real y estimado
      double jacc14 = ((double)e1+(double)e2-t) / t;

      out_str += fn1 + " " + fn2 + " " + std::to_string(jacc14) + "\n";
    }

    out[i_processed] = out_str;
  }


	for (size_t i_processed = 0; i_processed < card_name.size (); ++i_processed)
	{
      std::cout << out[i_processed];
	}



  return 0;
}
