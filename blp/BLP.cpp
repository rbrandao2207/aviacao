#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>

#include "BLP.hpp"

namespace ublas = boost::numeric::ublas;


BLP::BLP(const std::vector<double> init_guess, const double init_tetra_size)
{

  params_nbr = init_guess.size();
  ublas::vector<double> auxP;
  auxP.resize(params_nbr);
  for (unsigned i = 0; i < params_nbr+4; ++i) {
    P.push_back(auxP); /* init P0, ..., Pn+1 (=P[15], where params_nbr=15),
                          P^bar (P[params_nbr+1]), P* (...+2), P** (...+3) */ 
    y.push_back(0.); // init y
  }
  // initialize P's, P[0] at center
  for (unsigned i = 0; i < params_nbr+1; ++i) {
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
  ublas::vector<double> auxV;
  for (unsigned i = 0; i < params_nbr+4; ++i) {
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

void BLP::calc_objective(const double contract_tol, unsigned pt)
{
  // initialization
  // observe s_obs = s_obs_wg * pop_ave * mu (s_obs_wg is within group share)
  for (unsigned i = 0; i < N; ++i) {
    ln_s_obs[pt][i] = std::max(std::log(s_obs_wg[i] * (pop_ave[i] / 1e6) *\
					P[pt][14]),\
			       std::numeric_limits<double>::lowest());
  }

  // adjust lambda in [0,1]
  if (P[pt][13] > 1.) {
    P[pt][13] = 1.;
  } else if (P[pt][13] < 0.) {
    P[pt][13] = 0;
  }
    
  /// Berry's Contraction 
  bool conv_check = 0;
  while (!conv_check) {
    for (unsigned i = 0; i < N; ++i) {
      if (ln_s_obs[pt][i] == std::numeric_limits<double>::lowest()) {
        xi1[pt][i] = std::numeric_limits<double>::lowest();
      } else {
        xi1[pt][i] = xi0[pt][i] + P[pt][13] * (ln_s_obs[pt][i] -\
					       std::log(s_calc[pt][i]));
      }
    }

    // check for convergence
    for (unsigned i = 0; i <= N; ++i) {
      if (i == N) {
        conv_check = 1;
        break;
      }
      if (std::abs(xi1[pt][i] - xi0[pt][i]) < contract_tol) {
        continue;
      } else {
	/*DEBUG
	std::cout << std::abs(xi1[pt][i] - xi0[pt][i]) << '\t' << i << '\r'\
		  << std::flush;
	ENDDEBUG*/
        break;
      }
    }

    /// calc shares
    
    // calc s_ind1 & s_ind2 (exp(Xbeta....))
    // check header file for params mapping to P
    for (unsigned i = 0; i < N; ++i) {
      s_aux1[pt][i] = std::exp((X(i, 0) * P[pt][0] + X(i, 1) * P[pt][1] +\
  			      X(i, 2) * P[pt][2] + X(i, 3) * P[pt][3] +\
  			      X(i, 4) * P[pt][4] + X(i, 5) * P[pt][5] +\
  			      xi1[pt][i]) / P[pt][13]);
      s_aux2[pt][i] = std::exp((X(i, 0) * P[pt][6] + X(i, 1) * P[pt][7] +\
  			      X(i, 2) * P[pt][8] + X(i, 3) * P[pt][9] +\
  			      X(i, 4) * P[pt][10] + X(i, 5) * P[pt][11] +\
  			      xi1[pt][i]) / P[pt][13]);
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
        for (unsigned j = initial_aux_i; j < i; ++j) {
          D1[pt][j] = aux_D1;
          D2[pt][j] = aux_D2;
        }
      }
    }
  
    // compute model shares
    for (unsigned i = 0; i < N; ++i) {
      s_calc[pt][i] = P[pt][12] * ((s_aux1[pt][i] / D1[pt][i]) *\
  				 (std::pow(D1[pt][i], P[pt][13]) /\
  				  (1 + std::pow(D1[pt][i], P[pt][13])))) +\
        (1 - P[pt][12]) * ((s_aux2[pt][i] / D2[pt][i]) *\
  			 (std::pow(D2[pt][i], P[pt][13]) /\
  			  (1 + std::pow(D2[pt][i], P[pt][13]))));
    }

    // Update unobs util
    xi0[pt] = xi1[pt];
  }

  /// Compute objective function

  // adjust xi for numerical limits
  double lwr_bound = std::numeric_limits<double>::lowest() / 1e300;
  for (unsigned i = 0; i < N; ++i) {
    if (xi0[pt][i] < lwr_bound)
      xi0[pt][i] = lwr_bound;
  }

  // y = (1/N xi'Z)I(1/N Z'xi)
  ublas::identity_matrix<double> I (Z.size2());
  y[pt] = ublas::inner_prod(ublas::prod(1./N *\
					ublas::prod(ublas::trans(xi0[pt]),\
						    Z), I), 1./N *\
			    ublas::prod(ublas::trans(Z), xi0[pt]));
}

bool BLP::halt_check(const double NM_tol, unsigned iter_nbr)
{
  double y_sum = std::accumulate(y.begin(), y.end(), 0.);
  double y_avg = y_sum / y.size();
  double y_std = 0.;
  for (auto& y_elem : y) {
    y_std += std::pow(y_elem - y_avg, 2);
  }
  y_std = y_std / y.size();
  std::cout << y_std << '\t' << "# of iterations: " << iter_nbr << '\r' <<\
    std::flush;
  /*DEBUG
  std::cout << P[1][0] << '\t' << P[1][1] << '\t' << P[1][2] << '\t' <<\
            P[1][3] << '\r' << std::flush;
  ENDDEBUG*/
  if (y_std < NM_tol) {
    return true;
  } else {
    return false;
  }
}

void BLP::nelder_mead(const double contract_tol, const double alpha, const\
		      double beta, const double gamma, std::vector<unsigned>&\
		      points)
{
  //allocate
  unsigned h;
  unsigned l;

  // Basic functions
  auto get_h = [&] () {
		   std::vector<double>::iterator y_high;
		   y_high = std::max_element(y.begin(), y.end());
		   h = std::distance(y.begin(), y_high);
		   return h;
		 };
  auto get_l = [&] () {
		   std::vector<double>::iterator y_low;
		   y_low = std::min_element(y.begin(), y.end());
		   l = std::distance(y.begin(), y_low);
		   return l;
		 };
  auto P_bar = [&] () {
		 for (unsigned i = 0; i < params_nbr+1; ++i) {
		   if (i != h)
		     P[params_nbr+1] += P[i];
		 }
		 P[params_nbr+1] /= params_nbr;
		 return;
	       };
  auto reflection = [&] () {
		   P[params_nbr+2] = (1 + alpha) * P[params_nbr+1] - alpha *\
		     P[h];
		   this->calc_objective(contract_tol, params_nbr+2);
		   return;
		 };
  auto expansion = [&] () {
		  P[params_nbr+3] = gamma * P[params_nbr+2] + (1 - gamma) *\
		    P[params_nbr+1];
		  this->calc_objective(contract_tol, params_nbr+3);
		  return;
		};
  auto contraction = [&] () {
		  P[params_nbr+3] = beta * P[h] + (1 - beta) * P[params_nbr+1];
		  this->calc_objective(contract_tol, params_nbr+3);
		  return;
		};
  /* where P^bar = P[params_nbr+1],
     P* = P[params_nbr+2],
     P** = P[params_nbr+3] */
  
  // procedural functions
  auto step_1 = [&] () {
		  get_h();
		  P_bar();
		  reflection();
		};

  // NM algorithm (see Computer Journal 1965 page 2)
  step_1();
  get_l();
  if (y[params_nbr+2] < y[l]) { // y* < yl ?
    expansion();
    if (y[params_nbr+3] < y[l]) { // y** < yl ?
      P[h] = P[params_nbr+3]; // Ph = P**
      points = {h};
    } else {
      P[h] = P[params_nbr+2]; // Ph = P*
      points = {h};
    }
  } else {
    bool test = true;
    for (unsigned i = 0; i < params_nbr+1; ++i) {
      if (y[params_nbr+2] < y[i])
	test = false;
    }
    if (test) { // y* > yi all i ?
      if (y[params_nbr+2] < y[h]) { // y* < yh ?
	P[h] = P[params_nbr+2]; // replace Ph by P* (intermediate)
	this->calc_objective(contract_tol, h);
      }
      contraction();
      if (y[params_nbr+3] > y[h]) { // y** > yh ?
	for (unsigned i = 0; i < params_nbr+1; ++i) {
	  P[i] = (P[i] + P[l]) / 2; // Pi = (Pi + Pl) /2 all i
	  points = {};
	  for (unsigned i = 0; i < params_nbr+1; ++i) {
	    points.push_back(i);
	  }
	}
      } else {
	P[h] = P[params_nbr+3];  // Ph = P**
	points = {h};
      }
    } else {
      P[h] = P[params_nbr+2]; // Ph = P*
    }
  }
}
