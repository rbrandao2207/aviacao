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

  // estimation params:
  // initial guess ((alpha, beta)_r, gamma, lambda, mu)
  // initial tetrahedron "size" for Nelder Mead procedure
  // BLP contraction tolerance
  const std::vector<double> init_guess = {.1, .1, .1, .1, .1, .1, .1, .1, .1, \
					  .1, .1, .1, .5, .8, .01};
  const double init_tetra_size = {.1};
  const double contract_tol = .01;

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
    BLP inst_BLP(init_guess, init_tetra_size);
    // deserialize
    {
        std::ifstream ifs(persist_file);
        assert(ifs.is_open());
        boost::archive::text_iarchive ia(ifs);
        ia >> inst_BLP;
    }
    inst_BLP.allocate();
    // GMM estimation
    std::vector<std::thread> threads;
    while (true) { // TODO NM stop
      for (unsigned th = 0; th < inst_BLP.params_nbr + 1; ++th) {
        threads.push_back(std::thread(&BLP::calc_objective, std::ref(inst_BLP), contract_tol, th));
      }
      for (auto& thread : threads)
        thread.join();
      break;
    }
    inst_BLP.nelder_mead();

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
