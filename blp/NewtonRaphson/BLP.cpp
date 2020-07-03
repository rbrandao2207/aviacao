#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>

#include "BLP.hpp"

namespace ublas = boost::numeric::ublas;


BLP::BLP(const std::vector<double> init_guess, const double min_share_, const\
	 double contract_tol_, const double penalty_param1_, const double\
	 penalty_param2_, unsigned const max_threads)
{
  min_share = min_share_;
  contract_tol = contract_tol_;
  penalty_param1 = penalty_param1_;
  penalty_param2 = penalty_param2_;
  params_nbr = init_guess.size();
  ublas::vector<double> auxP;
  auxP.resize(params_nbr);
  for (unsigned i = 0; i <= params_nbr+1; ++i) {
    P.push_back(auxP); // init P0, P1
    y.push_back(0.); // init y
  }
  // initialize P's, P[0] at center
  for (unsigned j = 0; j < params_nbr; ++j) {
    P[0][j] = init_guess[j];
  }

  // init parallel params
  unsigned hardware_threads = std::thread::hardware_concurrency();
  num_threads = std::min(hardware_threads != 0 ? hardware_threads : 1,\
			 max_threads);
}

void BLP::allocate()
{
  grad.resize(params_nbr);
  // Initialize N, unobs util, and allocate other ublas objs
  N = s_obs_wg.size();
  ublas::vector<double> auxV;
  for (unsigned i = 0; i <= params_nbr+1; ++i) {
    xi0.push_back(auxV);
    xi0[i].resize(N);
    std::fill(xi0[i].data().begin(), xi0[i].data().end(), 0.);
    xi1.push_back(auxV);
    xi1[i].resize(N);
    std::fill(xi1[i].data().begin(), xi1[i].data().end(), 0.);
    s_aux1.push_back(auxV);
    s_aux1[i].resize(N);
    s_aux2.push_back(auxV);
    s_aux2[i].resize(N);
    D1.push_back(auxV);
    D1[i].resize(N);
    D2.push_back(auxV);
    D2[i].resize(N);
    s_calc.push_back(auxV);
    s_calc[i].resize(N);
    std::fill(s_calc[i].data().begin(), s_calc[i].data().end(), .1);
    ln_s_obs.push_back(auxV);
    ln_s_obs[i].resize(N);
  }
}

void BLP::calc_objective(unsigned pt)
{
  // initialization
  std::fill(xi0[pt].data().begin(), xi0[pt].data().end(), 0.);
  bool xi_nan = false;
  std::vector<std::thread> threads;
  unsigned j, k, block_size;
  block_size = N / num_threads;
  // observe s_obs = s_obs_wg * pop_ave * mu (s_obs_wg is within group share)
  auto ln_s_obs_L = [&] (unsigned begin, unsigned end) {
		      for (unsigned i = begin; i < end; ++i) {
			ln_s_obs[pt][i] = std::max(std::log(s_obs_wg[i] *\
							    (pop_ave[i] / 1e6)\
							    * P[pt][22]),\
						   std::log(min_share *\
							    (pop_ave[i] / 1e6)\
							    * P[pt][22]));
		    }
		  };
  j = 0;
  for (unsigned i = 0; i < (num_threads - 1); ++i) {
    k = j + block_size;
    threads.push_back(std::thread(ln_s_obs_L, j, k));
    j = k;
  }
  threads.push_back(std::thread(ln_s_obs_L, j, N));
  for (auto& thread : threads) {
    thread.join();
  }
  
  /// Berry's Contraction 
  bool conv_check = 0;
  while (!conv_check) {

    auto xi1_L = [&] (unsigned begin, unsigned end) {
		   for (unsigned i = begin; i < end; ++i) {
		     xi1[pt][i] = xi0[pt][i] + P[pt][21] * (ln_s_obs[pt][i] -\
					     std::log(s_calc[pt][i]));
		     if (std::isnan(xi1[pt][i])) {
		       xi_nan = true;
		       break;
		     }
		   }
		 };
    threads.clear();
    j = 0;
    for (unsigned i = 0; i < (num_threads - 1); ++i) {
      k = j + block_size;
      threads.push_back(std::thread(xi1_L, j, k));
      j = k;
    }
    threads.push_back(std::thread(xi1_L, j, N));
    for (auto& thread : threads) {
      thread.join();
    }
    // check for convergence
    if (xi_nan)
      break;
    for (unsigned i = 0; i <= N; ++i) {
      if (i == N) {
        conv_check = 1;
        break;
      }
      if (std::abs(xi1[pt][i] - xi0[pt][i]) < contract_tol) {
        continue;
      } else {
        break;
      }
    }
    /// calc shares
    // calc s_ind1 & s_ind2 (exp(Xbeta....))
    // check header file for params mapping to P
    auto s_aux_L = [&] (unsigned begin, unsigned end) {
		     for (unsigned i = begin; i < end; ++i) {
		       s_aux1[pt][i] = std::exp((X(i, 0) * P[pt][0] + X(i, 1) *\
						 P[pt][1] + X(i, 2) * P[pt][2] +\
						 X(i, 3) * P[pt][3] + X(i, 4) *\
						 P[pt][4] + X(i, 5) * P[pt][5] +\
						 X(i, 6) * P[pt][6] + X(i, 7) *\
						 P[pt][7] + X(i, 8) * P[pt][8] +\
						 X(i, 9) * P[pt][9] +\
						 xi1[pt][i]) / P[pt][21]);
		       s_aux2[pt][i] = std::exp((X(i, 0) * P[pt][10] + X(i, 1) *\
						 P[pt][11] + X(i, 2) * P[pt][12]\
						 + X(i, 3) * P[pt][13] + X(i, 4)\
						 * P[pt][14] + X(i, 5) *\
						 P[pt][15] + X(i, 6) *\
						 P[pt][16] + X(i, 7) * P[pt][17]\
						 + X(i, 8) * P[pt][18] + X(i, 9)\
						 * P[pt][19] + xi1[pt][i]) /\
						P[pt][21]);
		     }
		   };
    threads.clear();
    j = 0;
    for (unsigned i = 0; i < (num_threads - 1); ++i) {
      k = j + block_size;
      threads.push_back(std::thread(s_aux_L, j, k));
      j = k;
    }
    threads.push_back(std::thread(s_aux_L, j, N));
    for (auto& thread : threads) {
      thread.join();
    }
    // calc D1 and D2
    double aux_D1 = 0;
    double aux_D2 = 0;
    unsigned initial_aux_i = 0;
    unsigned aux_mkt_id = 0;
    for (unsigned i = 0; i < N; ++i) {
      if (aux_mkt_id == mkt_id[i]) {
        aux_D1 += s_aux1[pt][i];
        aux_D2 += s_aux2[pt][i];
      } else {
        for (unsigned j = initial_aux_i; j < i; ++j) {
          D1[pt][j] = aux_D1;
          D2[pt][j] = aux_D2;
        }
        aux_D1 = s_aux1[pt][i];
        aux_D2 = s_aux2[pt][i];
        initial_aux_i = i;
        aux_mkt_id = mkt_id[i];
      }
      if (i == N - 1) {
        for (unsigned j = initial_aux_i; j <= i; ++j) {
          D1[pt][j] = aux_D1;
          D2[pt][j] = aux_D2;
        }
      }
    }
  
    // compute model shares
    auto s_calc_L = [&] (unsigned begin, unsigned end) {
		      for (unsigned i = begin; i < end; ++i) {
			s_calc[pt][i] = P[pt][20] * ((s_aux1[pt][i] / D1[pt][i])\
						     * (std::pow(D1[pt][i],\
								 P[pt][21]) /\
							(1 +\
							 std::pow(D1[pt][i],\
								  P[pt][21]))))\
			  + (1 - P[pt][20]) * ((s_aux2[pt][i] / D2[pt][i]) *\
					       (std::pow(D2[pt][i], P[pt][21]) /\
						(1 + std::pow(D2[pt][i],\
							      P[pt][21]))));
		      }
		    };
    threads.clear();
    j = 0;
    for (unsigned i = 0; i < (num_threads - 1); ++i) {
      k = j + block_size;
      threads.push_back(std::thread(s_calc_L, j, k));
      j = k;
    }
    threads.push_back(std::thread(s_calc_L, j, N));
    for (auto& thread : threads) {
      thread.join();
    }
    // Update unobs util
    xi0[pt] = xi1[pt];
  }

  /// Compute objective function
  if (xi_nan) {
    y[pt] = std::numeric_limits<double>::max();
  } else {
    // y = (1/N xi'Z)I(1/N Z'xi)
    ublas::identity_matrix<double> I (Z.size2());
    y[pt] = ublas::inner_prod(ublas::prod(1./N *\
    					ublas::prod(ublas::trans(xi0[pt]),\
    						    Z), I), 1./N *\
    			    ublas::prod(ublas::trans(Z), xi0[pt]));
    
    /* Add penalty for constraints violation (gamma or lambda not in [0,1], or \
       mu < 0) */
    double penalty = {0.};
    if (P[pt][20] < 0.) {
      penalty += std::pow(1 + P[pt][20], penalty_param2);
    } else if (P[pt][20] > 1.) {
      penalty += std::pow(P[pt][20], penalty_param2);
    } else if (P[pt][21] < 0.) {
      penalty += std::pow(1 + P[pt][21], penalty_param2);
    } else if (P[pt][21] > 1.) {
      penalty += std::pow(P[pt][21], penalty_param2);
    } else if (P[pt][22] < 0.) {
      penalty += std::pow(1 + P[pt][22], penalty_param2);
    }
    penalty *= penalty_param1;
    y[pt] += penalty;
    if (std::isnan(y[pt]))
      y[pt] = std::numeric_limits<double>::max();
  }
}

void BLP::updatePs(const double inc)
{

  for (unsigned i = 1; i <= params_nbr; ++i) {
    P[i] = P[0];
    P[i][i-1] += inc;
  }
}

void BLP::step(const double inc, const double step_size, const double tol)
{
  halt_check = true;
  for (unsigned i = 0; i < params_nbr; ++i) {
    grad[i] = (y[i+1] - y[0]) / inc;
    if (halt_check == true && grad[i] > tol)
      halt_check = false;
    P[0][i] += grad[i] * step_size;
  }
}

void BLP::persist(const std::string persist_file2)
{
  std::ofstream fdesc;
  fdesc.open(persist_file2);
  assert(fdesc.is_open());
  for (unsigned i = 0; i < params_nbr; ++i) {
    fdesc << "P" << i << ": " << P[0][i] << '\n';
  }
  fdesc << "y: " << y[0] << '\n';
  fdesc.close();
  std::cout << "Finished params persistance in file " << persist_file2 <<\
    std::endl;
}
