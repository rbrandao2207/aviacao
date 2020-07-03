#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstring>
#include <fstream>
#include <limits>
#include <iostream>
#include <string>
#include <thread>
#include <valarray>
#include <vector>
#include <boost/archive/text_oarchive.hpp>

#include "BLP.hpp"
#include "GenArrays.hpp"


/* current run options:
   1) argv[1] genarrays
   2) argv[1] estimation
   3) argv[1] genarrays, argv[2] estimation */


int main(int argc, char* argv[])
{
  auto chrono_start = std::chrono::steady_clock::now();

  /* PARAMETERS */

  // estimation periods - currently quarters only
  const std::vector<std::string> dates = {"201801", "201802", "201803"};
  // run identifier
  const std::string run_id = "2018Q1";

  /// Genarrays
  // price bins
  const std::valarray<double> bins = {0, 100, 200, 300, 400, 500, 600, 700, 800,\
				      900, 1e3, 1.25e3, 1.5e3, 1.75e3, 2e3, 3e3,\
				      5e5};
  // population threshold
  const unsigned pop_thres = 5e5;

  // results directory
  const std::string results_dir = "results/";
  const std::string persist_file = results_dir + "arrays/" + run_id;
  const std::string persist_file2 = results_dir + "est_params/" + run_id;

  // maximum number of iterations
  const unsigned max_iter = std::numeric_limits<unsigned>::max();
  
  /// Estimation params:
  // initial guess ((alpha, beta)_r, gamma, lambda, mu)
  std::vector<double> init_guess = {.1, -.1, .1, -.1, .1, .1, .1, .1, .1,\
					  .1,\
					  .1, -.1, .1, -.1, .1, .1, .1, .1, .1,\
					  .1,\
					  .5, .8, .1};
  // minimum 'observed shares' for numerical feasibility
  const double min_share = {1e-20};
  // BLP contraction tolerance (BJ10 suggests 1e-12)
  const double contract_tol = {1e-4};
  // constrained optimization penalty
  const double penalty_param1 = {1e6};
  const unsigned penalty_param2 = {4}; // (must be even)

  // dx for numerical gradient
  const double inc = std::sqrt(std::numeric_limits<double>::min());

  const double step_size = 1e-2;
  const double tol = 1e-4;
  
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

  } else if ((argc > 1 && std::strcmp(argv[1], "estimation") == 0) ||\
	     (argc > 2 && std::strcmp(argv[1], "genarrays") == 0 &&\
	      std::strcmp(argv[2], "estimation") == 0)) {

    // instantiate
    BLP inst_BLP(init_guess, min_share, contract_tol, penalty_param1,\
		 penalty_param2);
    
    // deserialize
    {
        std::ifstream ifs(persist_file);
        assert(ifs.is_open());
        boost::archive::text_iarchive ia(ifs);
        ia >> inst_BLP;
	ifs.close();
    }
    inst_BLP.allocate();

    // GMM - Newton Raphson
    std::vector<unsigned> points;
    for (unsigned i = 0; i < inst_BLP.params_nbr+1; ++i) {
      points.push_back(i);
    }
    unsigned iter_nbr = 0;
    while (true) {
      inst_BLP.updatePs(inc);
      std::vector<std::thread> threads = {};
      for (auto& pt : points) {
        threads.push_back(std::thread(&BLP::calc_objective, std::ref(inst_BLP),\
				      pt));
      }
      for (auto& thread : threads) {
        thread.join();
      }
      inst_BLP.step(inc, step_size, tol);
      //if (inst_BLP.halt_check(NM_tol, iter_nbr))
      if (inst_BLP.halt_check) 
        break;
      ++iter_nbr;
      if (iter_nbr == max_iter)
	break;
    }
    // persist results
    inst_BLP.persist(persist_file2);
    std::cout << "# of iterations: " << iter_nbr << std::endl;

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