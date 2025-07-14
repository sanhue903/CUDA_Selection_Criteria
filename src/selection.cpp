#include "sketch/sketch.h"
#include <fstream>
#include <omp.h>
#include <iostream>
#include <seqan/seq_io.h>
#include <include/metrictime2.hpp>
#include <include/criteria_sketch.hpp>


// read_hll()

std::vector<uint64_t> read_smh(std::string path){
  //std::vector<uint64_t> smh_vector = smh->h_;
  //uint32_t smh_size = smh->m_;
        gzFile fp(gzopen(path.data(), "rb"));
        if(fp == nullptr) throw std::runtime_error(std::string("Could not open file at '") + path + "' for reading");

  //read
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


// ----------------------------------------------------------
// ------------------------- MAIN ---------------------------
// ----------------------------------------------------------

int main(int argc, char *argv[])
{
  // Paso 0: Cargar parámetros ##############################
  std::vector<std::string> files;

  //const uint k = 31;
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
        break;
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

  // Inicializar variables:
  omp_set_num_threads (threads);
  load_file_list (files, list_file);

  std::vector<std::pair<std::string, double>> card_name (files.size ());
  std::map<std::string, std::shared_ptr<sketch::hll_t>> name2hll14;
  std::string out[files.size()];

  // hll_a is cb+hll_a
  if (criterion == "hll_a"){
    // Read sketches from disk and order according to reported cardinality of primary hll sketch
    std::map<std::string, std::shared_ptr<sketch::hll_t>> name2aux;
    uint p = __builtin_ctz (aux_bytes);

    for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
    {
      std::string filename = files.at (i_processed);
      name2hll14[filename] = std::make_shared<sketch::hll_t> (14);
      name2aux[filename] = std::make_shared<sketch::hll_t> (p);
    }

    #pragma omp parallel for schedule(dynamic)
    for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
    {
      std::string filename = files.at (i_processed);
      name2hll14[filename] = std::make_shared<sketch::hll_t>(filename + ".hll");
      name2aux[filename] = std::make_shared<sketch::hll_t>(filename + ".hll_" + std::to_string(p));
      auto c = name2hll14[filename]->report ();
      card_name.at (i_processed) = std::make_pair (filename, c);
    }

    std::sort (card_name.begin (), card_name.end (),
               [](const std::pair<std::string, double> &x,
                  const std::pair<std::string, double> &y)
               {
               return x.second < y.second;
               });

    // Comparison
    #pragma omp parallel for schedule(dynamic)
    for (size_t i_processed = 0; i_processed < card_name.size () - 1; ++i_processed)
    {
      std::string out_str;
      std::string fn1 = card_name[i_processed].first;
      size_t e1 = card_name[i_processed].second;
      size_t k = i_processed + 1;
      for (; k < card_name.size (); ++k)
      {
        std::string fn2 = card_name[k].first;
        size_t e2 = card_name[k].second;
        if(e2 == 0)continue;
        bool CB_selects = CB(threshold, e1, e2);
        if (!CB_selects) break;
        bool hlla_selects = hll_a(threshold, e1, e2, name2aux[fn1], name2aux[fn2], p, z_score);
        if (!hlla_selects) continue;
        double t = name2hll14[fn1]->union_size(*name2hll14[fn2]);
        double jacc14 = ((double)e1+(double)e2-t) / t;
        if (jacc14 >= threshold) out_str += fn1 + " " + fn2 + " " + std::to_string(jacc14) + "\n";
      }
      out[i_processed] = out_str;
    }

  }else if (criterion == "hll_an"){
    // Read sketches from disk and order according to reported cardinality of primary hll sketch
    std::map<std::string, std::shared_ptr<sketch::hll_t>> name2aux;
    uint p = __builtin_ctz (aux_bytes);

    for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
    {
      std::string filename = files.at (i_processed);
      name2hll14[filename] = std::make_shared<sketch::hll_t> (14);
      name2aux[filename] = std::make_shared<sketch::hll_t> (p);
    }

    #pragma omp parallel for schedule(dynamic)
    for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
    {
      std::string filename = files.at (i_processed);
      name2hll14[filename] = std::make_shared<sketch::hll_t>(filename + ".hll");
      name2aux[filename] = std::make_shared<sketch::hll_t>(filename + ".hll_" + std::to_string(p));
      auto c = name2hll14[filename]->report ();
      card_name.at (i_processed) = std::make_pair (filename, c);
    }

    std::sort (card_name.begin (), card_name.end (),
               [](const std::pair<std::string, double> &x,
                  const std::pair<std::string, double> &y)
               {
               return x.second < y.second;
               });


    // Comparison 
    #pragma omp parallel for schedule(dynamic)
    for (size_t i_processed = 0; i_processed < card_name.size () - 1; ++i_processed)
    {
      std::string out_str;
      std::string fn1 = card_name[i_processed].first;
      size_t e1 = card_name[i_processed].second;
      size_t k = i_processed + 1;
      for (; k < card_name.size (); ++k)
      {
        std::string fn2 = card_name[k].first;
        size_t e2 = card_name[k].second;
        if(e2 == 0)continue;
        bool CB_selects = CB(threshold, e1, e2);
        if (!CB_selects) break;
        bool hllan_selects = hll_an(threshold, e1, e2, name2aux[fn1], name2aux[fn2], p, z_score, order_n);
        if (!hllan_selects) continue;
        double t = name2hll14[fn1]->union_size(*name2hll14[fn2]);
        double jacc14 = ((double)e1+(double)e2-t) / t;
        if (jacc14 >= threshold) out_str += fn1 + " " + fn2 + " " + std::to_string(jacc14) + "\n";
      }
      out[i_processed] = out_str;
    }
  }else if (criterion == "smh_a"){
    // Read sketches from disk and order according to reported cardinality of primary hll sketch
    std::map<std::string, std::vector<uint64_t>> name2aux;
    uint m = aux_bytes/8;

    for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
    {
      std::string filename = files.at (i_processed);
      std::vector<uint64_t> v(m);
      name2hll14[filename] = std::make_shared<sketch::hll_t> (14);
      name2aux[filename] = v;
    }

    #pragma omp parallel for schedule(dynamic)
    for (size_t i_processed = 0; i_processed < files.size (); ++i_processed)
    {
      std::string filename = files.at (i_processed);
      name2hll14[filename] = std::make_shared<sketch::hll_t>(filename + ".hll");
      name2aux[filename] = read_smh(filename + ".smh"+std::to_string(m));
      auto c = name2hll14[filename]->report ();
      card_name.at (i_processed) = std::make_pair (filename, c);
    }

    std::sort (card_name.begin (), card_name.end (),
               [](const std::pair<std::string, double> &x,
                  const std::pair<std::string, double> &y)
               {
               return x.second < y.second;
               });

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

    // Comparison 
    #pragma omp parallel for schedule(dynamic)
    for (size_t i_processed = 0; i_processed < card_name.size () - 1; ++i_processed)
    {
      std::string out_str;
      std::string fn1 = card_name[i_processed].first;
      size_t e1 = card_name[i_processed].second;
      size_t k = i_processed + 1;
      for (; k < card_name.size (); ++k)
      {
        std::string fn2 = card_name[k].first;
        size_t e2 = card_name[k].second;
        if(e2 == 0)continue;
        bool CB_selects = CB(threshold, e1, e2);
        if (!CB_selects) break;
        bool smha_selects = smh_a(name2aux[fn1], name2aux[fn2], n_rows, n_bands);
        if (!smha_selects) continue;
        double t = name2hll14[fn1]->union_size(*name2hll14[fn2]);
        double jacc14 = ((double)e1+(double)e2-t) / t;
        if (jacc14 >= threshold) out_str += fn1 + " " + fn2 + " " + std::to_string(jacc14) + "\n";
      }
      out[i_processed] = out_str;
    }
  }else{
    std::cout << "Option -c invalid. The accepted criteria are hll_a, hll_an and smh_a.\n";
  }


  for (size_t i_processed = 0; i_processed < card_name.size (); ++i_processed)
  {
      std::cout << out[i_processed];
  }


  return 0;
}
