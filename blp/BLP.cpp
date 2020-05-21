#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>

#include "BLP.hpp"

using namespace boost::numeric::ublas;


BLP::BLP(const std::vector<double> init_guess, const double init_tetra_size)
{

  params_nbr = init_guess.size();
  vector<double> auxP;
  auxP.resize(params_nbr);
  for (unsigned i = 0; i < params_nbr + 1; ++i) {
    P.push_back(auxP);
  }
  // initialize P's, P[0] at center
  for (unsigned i = 0; i < params_nbr + 1; ++i) {
    for (unsigned j = 0; j < params_nbr; ++j) {
      if (i == 0 || j != i) {
	P[i][j] = init_guess[j];
      } else {
	P[i][j] = init_guess[j] + init_tetra_size;
      }	
    }
  }
}

void BLP::allocate()
{
  // Initialize N, unobs util, and allocate other ublas objs
  N = s_obs_wg.size();
  vector<double> auxV;
  for (unsigned i = 0; i < params_nbr + 1; ++i) {
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
    ln_s_obs.push_back(auxV);
    ln_s_obs[i].resize(N);
  }
}

void BLP::calc_objective(const double contract_tol, unsigned th)
{
  // initialization
  // observe s_obs = s_obs_wg * pop_ave * mu (s_obs_wg is within group share)
  for (unsigned i = 0; i < N; ++i) {
    ln_s_obs[th][i] = std::max(std::log(s_obs_wg[i] * pop_ave[i] * P[th][14]),\
			   std::numeric_limits<double>::lowest());
  }

    
  // Berry's Contraction 
  bool conv_check = 0;
  while (!conv_check) {
    for (unsigned i = 0; i < N; ++i) {
      if (ln_s_obs[th][i] == std::numeric_limits<double>::lowest()) {
        xi1[th][i] = std::numeric_limits<double>::lowest();
      } else {
        xi1[th][i] = xi0[th][i] + P[th][13] * (ln_s_obs[th][i] -\
					       std::log(s_calc[th][i]));
      }
    }

    // check for convergence
    for (unsigned i = 0; i <= N; ++i) {
      if (i == N) {
        conv_check = 1;
        break;
      }
      if (std::abs(xi1[th][i] - xi0[th][i]) < contract_tol) {
        continue;
    		
      } else {
        break;
      }
    }

    // Begin calc shares
    
    // calc s_ind1 & s_ind2 (exp(Xbeta....))
    // check header file for params mapping to P
    for (unsigned i = 0; i < N; ++i) {
      s_aux1[th][i] = std::exp((X(i, 0) * P[th][0] + X(i, 1) * P[th][1] +\
  			      X(i, 2) * P[th][2] + X(i, 3) * P[th][3] +\
  			      X(i, 4) * P[th][4] + X(i, 5) * P[th][5] +\
  			      xi1[th][i]) / P[th][13]);
      s_aux2[th][i] = std::exp((X(i, 0) * P[th][6] + X(i, 1) * P[th][7] +\
  			      X(i, 2) * P[th][8] + X(i, 3) * P[th][9] +\
  			      X(i, 4) * P[th][10] + X(i, 5) * P[th][11] +\
  			      xi1[th][i]) / P[th][13]);
    }
  
    //calc D1 and D2
    double aux_D1 = 0;
    double aux_D2 = 0;
    unsigned initial_aux_i = 0;
    unsigned aux_mkt_id = 0;
    for (unsigned i = 0; i < N; ++i) {
      if (aux_mkt_id == mkt_id[i]) {
        aux_D1 += s_aux1[th][i];
        aux_D2 += s_aux2[th][i];
      } else {
        for (unsigned j = initial_aux_i; j < i; ++j) {
          D1[th][j] = aux_D1;
          D2[th][j] = aux_D2;
        }
        aux_D1 = s_aux1[th][i];
        aux_D2 = s_aux2[th][i];
        initial_aux_i = i;
        aux_mkt_id = mkt_id[i];
      }
      if (i == N - 1) {
        for (unsigned j = initial_aux_i; j < i; ++j) {
          D1[th][j] = aux_D1;
          D2[th][j] = aux_D2;
        }
      }
    }
  
    // compute model shares
    for (unsigned i = 0; i < N; ++i) {
      s_calc[th][i] = P[th][12] * ((s_aux1[th][i] / D1[th][i]) *\
  				 (std::pow(D1[th][i], P[th][13]) /\
  				  (1 + std::pow(D1[th][i], P[th][13])))) +\
        (1 - P[th][12]) * ((s_aux2[th][i] / D2[th][i]) *\
  			 (std::pow(D2[th][i], P[th][13]) /\
  			  (1 + std::pow(D2[th][i], P[th][13]))));
    }

    // DEBUG
    double error = 0;
    for (unsigned i = 0; i < N; ++i) {
      error += std::pow(xi0[th][i] - xi1[th][i], 2);
    }
    std::cout << error << std::endl;
    // ENDDEBUG
    
    // Update unobs util
    xi0[th] = xi1[th];
  }
}

void BLP::nelder_mead()
{
  unsigned x = 0;
}
