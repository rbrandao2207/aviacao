#ifndef BLPHEADERDEF 
#define BLPHEADERDEF

#include <string>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;


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

  unsigned params_nbr;

  void allocate();
  void calc_objective(const double contract_tol, unsigned th);
  void nelder_mead();

private:
  // Exogenous vars 

  unsigned N;

  vector<double> s_obs_wg;
  vector<double> mkt_id;
  vector<double> pop_ave;
  /* X matrix structure: 0 constant
                         1 fare (bin center R$/1000)
                         2 distance (km / 1000)
                         3 distance^2
                         4 carrier dummy
                         5 time dummy */
  matrix<double> X;
  /* Z matrix structure: 0 constant
                         1 aviation fuel price (querosene, instrument)
                         2 distance (km / 1000)
                         3 distance^2
                         4 carrier dummy
                         5 time dummy */
  matrix<double> Z;

  // PARAMS (15 total): 
  
  // 0: beta1_0, 1: beta1_1, 2: beta1_2, 3: beta1_3, 4: beta1_4, 5: beta1_5;
  // 6: beta2_0, 7: beta2_1, 8: beta2_2, 9: beta2_3, 10: beta2_4, 11: beta2_5;
  // 12: gamma, 13: lambda, 14: mu;
  // Nelder mead points contain params, N+1=16 total, std container for those:
  std::vector<vector<double>> P;

  // Calc vars
  
  // unobserved utility
  std::vector<vector<double>> xi0;
  std::vector<vector<double>> xi1;
  // other vars
  std::vector<vector<double>> s_aux1;
  std::vector<vector<double>> s_aux2;
  std::vector<vector<double>> D1;
  std::vector<vector<double>> D2;
  std::vector<vector<double>> s_calc;
  std::vector<vector<double>> ln_s_obs;
  // objective functions values
  std::vector<double> y;
};

#endif
