#include <algorithm>
#include <cmath>
#include <future>
#include <limits>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>

#include "BLP.hpp"

using namespace boost::numeric


BLP::BLP(const std::vector<double> init_guess, const double init_tetra_size)
{
  params_nbr = init_guess.size();
  ublas::vector auxP;
  auxP.resize(params_nbr);
  for (unsigned i = 0; i < params_nbr + 1; ++i) {
    P.push_back(auxP);
  }
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
  /* Initialize N, unobs util, and allocate other ublas objs */
  N = s_obs_wg.size();
  for (unsigned i = 0; i < params_nbr + 1; ++i) {
    xi0[i].resize(N);
    std::fill(xi0[i].data().begin(), xi0[i].data().end(), 0.);
    xi1[i].resize(N);
    std::fill(xi1[i].data().begin(), xi1[i].data().end(), 0.);
    s_aux1[i].resize(N);
    s_aux2[i].resize(N);
    D1[i].resize(N);
    D2[i].resize(N);
    s_calc[i].resize(N);
    ln_s_obs[i].resize(N);
  }
}

void BLP::gmm(const double& contract_tol)
{
  //TODO move all async to main, mantain two member functions with i arg (i -> P[i])
  //calc_objective and NM
  std::vector<std::future<int>> futures;
  for (unsigned i = 0; i < params_nbr + 1; ++i) {
    futures.push_back(std::async(calc_objective, std::ref(P[i])));
  }
		      
}

void BLP::calc_objective()
{
  contraction
  1 / N //continue from here
}

void BLP::contraction(const double& contract_tol)
{
  // initialization
  // observe s_obs = s_obs_wg * pop_ave * mu, where s_obs_wg is within group share
  for (unsigned i = 0; i < N; ++i) {
    ln_s_obs[i] = std::max(std::log(s_obs_wg[i] * pop_ave[i] * mu), \
			   std::numeric_limits<double>::lowest());
  }
  this->BLP::calc_shares();

  bool conv_check = 0;
  double s_calc_d; //DEBUG
  while (!conv_check) {
    for (unsigned i = 0; i < xi1.size(); ++i) {
      if (ln_s_obs[i] == std::numeric_limits<double>::lowest()) {
        xi1[i] = std::numeric_limits<double>::lowest();
      } else {
        xi1[i] = xi0[i] + lambda * (ln_s_obs[i] - std::log(s_calc[i]));
      }
      /*DEBUG
      if (i == 4807)
          s_calc_d = std::log(s_calc[i]);
      //ENDDEBUG*/
    }

    // check for convergence
    for (unsigned i = 0; i <= xi1.size(); ++i) {
      if (i == xi1.size()) {
        conv_check = 1;
        break;
      }
      if (std::abs(xi1[i] - xi0[i]) < contract_tol) {
        continue;
    		
      } else {
        break;
      }
    }

    // update vars
    this->BLP::calc_shares();
    xi0 = xi1;
  }
}

void BLP::calc_shares(const double& beta1_0, //TODO fill others)
{
  //calc s_ind1 & s_ind2 (exp(Xbeta....))
  for (unsigned i = 0; i < N; ++i) {
    /*DEBUG
    if (i == 4807)
        std::cout << X(i,0) << '\t' << X(i,1) << '\t' << X(i,2) << '\t' << \
                X(i,3) << '\t'  << X(i,4) << '\t'  << X(i,5) << '\t'  << \
                xi1[i] << std::endl;
    //ENDDEBUG*/
    s_aux1[i] = std::exp((X(i, 0) * beta1_0 + X(i, 1) * beta1_1 + X(i, 2) * \
			  beta1_2 + X(i, 3) * beta1_3 + X(i, 4) * beta1_4 + \
			  X(i, 5) * beta1_5 + xi1[i]) / lambda);
    s_aux2[i] = std::exp((X(i, 0) * beta2_0 + X(i, 1) * beta2_1 + X(i, 2) * \
			  beta2_2 + X(i, 3) * beta2_3 + X(i, 4) * beta2_4 + \
			  X(i, 5) * beta2_5 + xi1[i]) / lambda);
  }

  //calc D1 and D2
  double aux_D1 = 0;
  double aux_D2 = 0;
  unsigned initial_aux_i = 0;
  unsigned aux_mkt_id = 0;
  for (unsigned i = 0; i < mkt_id.size(); ++i) {
    if (aux_mkt_id == mkt_id[i]) {
      aux_D1 += s_aux1[i];
      aux_D2 += s_aux2[i];
    } else {
      for (unsigned j = initial_aux_i; j < i; ++j) {
        D1[j] = aux_D1;
        D2[j] = aux_D2;
      }
      aux_D1 = s_aux1[i];
      aux_D2 = s_aux2[i];
      initial_aux_i = i;
      aux_mkt_id = mkt_id[i];
    }
    if (i == mkt_id.size() - 1) {
      for (unsigned j = initial_aux_i; j < i; ++j) {
        D1[j] = aux_D1;
        D2[j] = aux_D2;
      }
    }
  }

  // calc model shares
  for (unsigned i = 0; i < s_calc.size(); ++i) {
    s_calc[i] = gamma * ((s_aux1[i] / D1[i]) * (std::pow(D1[i], lambda) / \
                (1 + std::pow(D1[i], lambda)))) + (1 - gamma) * ((s_aux2[i] / \
                D2[i]) * (std::pow(D2[i], lambda) / (1 + std::pow(D2[i], \
                lambda))));
    /*DEBUG
    if (std::isnan(s_calc[i]))
        std::cout << "nan" << std::endl;
    //ENDDEBUG*/
  }
}
