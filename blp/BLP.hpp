#ifndef BLPHEADERDEF 
#define BLPHEADERDEF

#include <string>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;


class BLP
{
public:
  
  BLP(const std::vector<double> init_guess, const double min_share_, const\
      double contract_tol_, const double penalty_param1_, const unsigned\
      penalty_param2_, const double init_tetra_size1, const double\
      init_tetra_size2);

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & s_obs_wg;
    ar & mkt_id;
    ar & pop_ave;
    ar & X;
    ar & Z;
  }

  unsigned params_nbr;

  void allocate();
  void calc_objective(unsigned iter_nbr, unsigned pt);
  bool halt_check(const double NM_tol, unsigned iter_nbr);
  void nelder_mead(unsigned iter_nbr, const double alpha, const double beta,\
		   const double gamma, std::vector<unsigned>& points);
  /* alpha: reflection coeff
     beta:  contraction coeff
     gamma: expansion coeff*/
  void persist(const std::string persist_file2);

private:

  // Params
  double min_share;
  double contract_tol;
  double penalty_param1;
  unsigned penalty_param2;
  // Exogenous vars 
  unsigned N;
  ublas::vector<double> s_obs_wg;
  ublas::vector<double> mkt_id;
  ublas::vector<double> pop_ave;
  ublas::matrix<double> X;
  /* X matrix structure: 0 constant
                         1 fare (bin center R$/1000)
                         2 distance (km / 1000)
                         3 distance^2
                         4 carrier dummy
                         5 time dummy */
  ublas::matrix<double> Z;
  /* Z matrix structure: 0 constant
                         1 aviation fuel price (querosene, instrument)
                         2 distance (km / 1000)
                         3 distance^2
                         4 carrier dummy
                         5 time dummy */


  std::vector<ublas::vector<double>> P;
  /* Params are 15 total:
     0: beta1_0, 1: beta1_1, 2: beta1_2, 3: beta1_3, 4: beta1_4, 5: beta1_5;
     6: beta2_0, 7: beta2_1, 8: beta2_2, 9: beta2_3, 10: beta2_4, 11: beta2_5;
     12: gamma, 13: lambda, 14: mu 
     NM procedures takes number of params + 1 = 16 P's, + 
       P^bar (params_nbr+1), P* (params_nbr+2) and P** (params_nbr+3) */

  // Objective function values
  std::vector<double> y;
  
  // Calc vars

  // unobserved utility
  std::vector<ublas::vector<double>> xi0;
  std::vector<ublas::vector<double>> xi1;
  // other vars
  std::vector<ublas::vector<double>> s_aux1;
  std::vector<ublas::vector<double>> s_aux2;
  std::vector<ublas::vector<double>> D1;
  std::vector<ublas::vector<double>> D2;
  std::vector<ublas::vector<double>> s_calc;
  std::vector<ublas::vector<double>> ln_s_obs;
};

#endif
