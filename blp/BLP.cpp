#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "BLP.hpp"


BLP::BLP(const std::vector<double> init_guess, const double init_tetra_size)
{
  P0.resize(init_guess.size());
  for (unsigned i = 0; i < init_guess.size(); ++i) {
    P0[i] = init_guess[i];
  }
  P1 = P0; P1[1] += init_tetra_size;
  P2 = P0; P2[2] += init_tetra_size;
  P3 = P0; P3[3] += init_tetra_size;
  P4 = P0; P4[4] += init_tetra_size;
  P5 = P0; P5[5] += init_tetra_size;
  P6 = P0; P6[6] += init_tetra_size;
  P7 = P0; P7[7] += init_tetra_size;
  P8 = P0; P8[8] += init_tetra_size;
  P9 = P0; P9[9] += init_tetra_size;
  P10 = P0; P10[10] += init_tetra_size;
  P11 = P0; P11[11] += init_tetra_size;
  P12 = P0; P12[12] += init_tetra_size;
  P13 = P0; P13[13] += init_tetra_size;
  P14 = P0; P14[14] += init_tetra_size;
  P15 = P0; P15[15] += init_tetra_size;
}

void BLP::allocate()
{
  /* Initialize N, unobs util, and allocate other ublas objs */
  N = s_obs_wg.size();
  // xi0
  P0_xi0.resize(N); P1_xi0.resize(N); P2_xi0.resize(N); P3_xi0.resize(N);
  P4_xi0.resize(N); P5_xi0.resize(N); P6_xi0.resize(N);	P7_xi0.resize(N);
  P8_xi0.resize(N); P9_xi0.resize(N); P10_xi0.resize(N); P11_xi0.resize(N);
  P12_xi0.resize(N); P13_xi0.resize(N); P14_xi0.resize(N); P15_xi0;
  std::fill(P0_xi0.data().begin(), P0_xi0.data().end(), 0.);
  std::fill(P1_xi0.data().begin(), P1_xi0.data().end(), 0.);
  std::fill(P2_xi0.data().begin(), P2_xi0.data().end(), 0.);
  std::fill(P3_xi0.data().begin(), P3_xi0.data().end(), 0.);
  std::fill(P4_xi0.data().begin(), P4_xi0.data().end(), 0.);
  std::fill(P5_xi0.data().begin(), P5_xi0.data().end(), 0.);
  std::fill(P6_xi0.data().begin(), P6_xi0.data().end(), 0.);
  std::fill(P7_xi0.data().begin(), P7_xi0.data().end(), 0.);
  std::fill(P8_xi0.data().begin(), P8_xi0.data().end(), 0.);
  std::fill(P9_xi0.data().begin(), P9_xi0.data().end(), 0.);
  std::fill(P10_xi0.data().begin(), P10_xi0.data().end(), 0.);
  std::fill(P11_xi0.data().begin(), P11_xi0.data().end(), 0.);
  std::fill(P12_xi0.data().begin(), P12_xi0.data().end(), 0.);
  std::fill(P13_xi0.data().begin(), P13_xi0.data().end(), 0.);
  std::fill(P14_xi0.data().begin(), P14_xi0.data().end(), 0.);
  std::fill(P15_xi0.data().begin(), P15_xi0.data().end(), 0.);
  // xi1
  P0_xi1.resize(N); P1_xi1.resize(N); P2_xi1.resize(N); P3_xi1.resize(N);
  P4_xi1.resize(N); P5_xi1.resize(N); P6_xi1.resize(N);	P7_xi1.resize(N);
  P8_xi1.resize(N); P9_xi1.resize(N); P10_xi1.resize(N); P11_xi1.resize(N);
  P12_xi1.resize(N); P13_xi1.resize(N); P14_xi1.resize(N); P15_xi1;
  std::fill(P0_xi1.data().begin(), P0_xi1.data().end(), 0.);
  std::fill(P1_xi1.data().begin(), P1_xi1.data().end(), 0.);
  std::fill(P2_xi1.data().begin(), P2_xi1.data().end(), 0.);
  std::fill(P3_xi1.data().begin(), P3_xi1.data().end(), 0.);
  std::fill(P4_xi1.data().begin(), P4_xi1.data().end(), 0.);
  std::fill(P5_xi1.data().begin(), P5_xi1.data().end(), 0.);
  std::fill(P6_xi1.data().begin(), P6_xi1.data().end(), 0.);
  std::fill(P7_xi1.data().begin(), P7_xi1.data().end(), 0.);
  std::fill(P8_xi1.data().begin(), P8_xi1.data().end(), 0.);
  std::fill(P9_xi1.data().begin(), P9_xi1.data().end(), 0.);
  std::fill(P10_xi1.data().begin(), P10_xi1.data().end(), 0.);
  std::fill(P11_xi1.data().begin(), P11_xi1.data().end(), 0.);
  std::fill(P12_xi1.data().begin(), P12_xi1.data().end(), 0.);
  std::fill(P13_xi1.data().begin(), P13_xi1.data().end(), 0.);
  std::fill(P14_xi1.data().begin(), P14_xi1.data().end(), 0.);
  std::fill(P15_xi1.data().begin(), P15_xi1.data().end(), 0.);
  
  xi0.resize(N);
  xi1.resize(N);
  std::fill(xi0.data().begin(), xi0.data().end(), 0.);
  std::fill(xi1.data().begin(), xi1.data().end(), 0.);
  s_aux1.resize(N);
  s_aux2.resize(N);
  D1.resize(N);
  D2.resize(N);
  s_calc.resize(N);
  ln_s_obs.resize(N);
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

void BLP::calc_objective()
{
  1 / N //continue from here
}

void BLP::gmm(const double& contract_tol)
{
}
