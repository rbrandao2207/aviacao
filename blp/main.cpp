#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <valarray>
#include <vector>
#include <boost/archive/text_oarchive.hpp>

#include "BLP.hpp"
#include "GenArrays.hpp"


// current run options:
// 1) argv[1] genarrays
// 2) argv[1] estimation
// 3) argv[1] genarrays & argv[2] estimation


int main(int argc, char* argv[])
{
  auto chrono_start = std::chrono::steady_clock::now();

  /* PARAMETERS */

  // estimation periods - enter all periods in tuple
  const std::vector<std::string> dates = {"201801", "201802", "201803"};
  const std::string run_id = "01";

  // price bins
  const std::valarray<double> bins = {0, 200, 400, 600, 800, 1000, 1e5};

  // population threshold
  const unsigned pop_thres = 1000000;

  // results directory
  const std::string results_dir = "results/";
  const std::string persist_file = results_dir + "arrays/" + run_id;

  /// Estimation params:
  // initial guess ((alpha, beta)_r, gamma, lambda, mu)
  const std::vector<double> init_guess = {.1, .1, .1, .1, .1, .1, .1, .1, .1, .1,\
					  .1, .1, .5, .8, .01};
  // BLP contraction tolerance
  const double contract_tol = {.01};
  // initial tetrahedron "size" for Nelder Mead procedure
  const double init_tetra_size = {.1};
  // NM coefficients
  const double alpha = {1.}; // reflection, alpha > 0
  const double beta = {.5};  // contraction, beta in [0,1]
  const double gamma = {2.}; // expansion, gamma > 1
  
  /* END OF PARAMETERS */

  if (argc > 1 && std::strcmp(argv[1], "genarrays") == 0) {
    GenArrays inst_GA(dates, bins, pop_thres);
    inst_GA.gen_instruments();
    inst_GA.gen_arrays();    
    // serialize
    {
        std::remove(persist_file.c_str());
        std::ofstream ofs(persist_file);
        assert(ofs.is_open());
        boost::archive::text_oarchive oa(ofs);
        oa << inst_GA;
    }

  } else if ((argc > 1 && std::strcmp(argv[1], "estimation") == 0) || \
	     (argc > 2 && std::strcmp(argv[1], "genarrays") == 0 && \
	      std::strcmp(argv[2], "estimation") == 0)) {

    // instantiate
    BLP inst_BLP(init_guess, init_tetra_size);
    
    // deserialize
    {
        std::ifstream ifs(persist_file);
        assert(ifs.is_open());
        boost::archive::text_iarchive ia(ifs);
        ia >> inst_BLP;
    }
    inst_BLP.allocate();
    
    // GMM
    std::vector<unsigned> points;
    for (unsigned i = 0; i < inst_BLP.params_nbr + 1; ++i) {
      points.push_back(i);
    }
    while (True) {
      // calc simplex
      std::vector<std::thread> threads;
      for (auto& pt : points) {
        threads.push_back(std::thread(&BLP::calc_objective, std::ref(inst_BLP),\
      				    contract_tol, pt));
      }
      for (auto& thread : threads) {
        thread.join();
      }
      // NM procedure
      inst_BLP.nelder_mead(contract_tol, alpha, beta, gamma, points);
      break;
    }

  } else {
    std::cout << "Invalid args!" << std::endl;
    throw std::runtime_error("aborting");
  }


  // finish chrono
  auto chrono_end = std::chrono::steady_clock::now();
  auto time_diff = chrono_end - chrono_start;
  std::string time_fpersist = results_dir + "ellapsed_time";
  std::remove(time_fpersist.c_str());
  std::ofstream fdesc_time;
  fdesc_time.open(time_fpersist);
  fdesc_time << "Last run duration: " << \
    std::chrono::duration_cast<std::chrono::minutes> (time_diff).count() << \
    " mins" << std::endl;
  fdesc_time.close();
  std::cout << "Elapsed time: " << \
    std::chrono::duration_cast<std::chrono::minutes> (time_diff).count() << \
    " mins" << std::endl;

  return 0;
}
