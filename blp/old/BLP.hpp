#ifndef BLPHEADERDEF 
#define BLPHEADERDEF

#include <string>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric


class BLP
{

public:
  BLP(const std::vector<double> init_guess, const double init_tetra_size);
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & s_obs_wg;
    ar & mkt_id;
    ar & pop_ave;
    ar & X;
    ar & Z;
  }
  void allocate();
  void calc_shares();
  void contraction(const& double contract_tol);
  void calc_objective();
  void gmm(const& double contract_tol);

private:

  /* Exogenous vars */

  unsigned params_nbr;
  unsigned N;
  ublas::vector<double> s_obs_wg;
  ublas::vector<double> mkt_id;
  ublas::vector<double> pop_ave;
  /* X matrix structure: 0 constant
                         1 fare (bin center R$/1000)
                         2 distance (km / 1000)
                         3 distance^2
                         4 carrier dummy
                         5 time dummy
  */
  ublas::matrix<double> X;
  /* Z matrix structure: 0 constant
                         1 aviation fuel price (querosene, instrument)
                         2 distance (km / 1000)
                         3 distance^2
                         4 carrier dummy
                         5 time dummy
  */
  ublas::matrix<double> Z;

  /* PARAMS (15 total): */
  
  // beta1_0, beta1_1, beta1_2, beta1_3, beta1_4, beta1_5;
  // beta2_0, beta2_1, beta2_2, beta2_3, beta2_4o, beta2_5;
  // gamma, lambda, mu;
  // Nelder mead points contain params, N+1=16 total, std container for those:
  std::vector<ublas::vector<double>> P;

  /* Calc vars */
  
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
  // objective functions values
  std::vector<double> y;
};

#endif
