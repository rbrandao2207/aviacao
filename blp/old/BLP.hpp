#ifndef BLPHEADERDEF 
#define BLPHEADERDEF

#include <string>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas


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
  
  unsigned N;
  vector<double> s_obs_wg;
  vector<double> mkt_id;
  vector<double> pop_ave;
  /* X matrix structure: 0 constant
                         1 fare (bin center R$/1000)
                         2 distance (km / 1000)
                         3 distance^2
                         4 carrier dummy
                         5 time dummy
  */
  matrix<double> X;
  /* Z matrix structure: 0 constant
                         1 aviation fuel price (querosene, instrument)
                         2 distance (km / 1000)
                         3 distance^2
                         4 carrier dummy
                         5 time dummy
  */
  matrix<double> Z;

  /* PARAMS (15 total): */
  
  // beta1_0, beta1_1, beta1_2, beta1_3, beta1_4, beta1_5;
  // beta2_0, beta2_1, beta2_2, beta2_3, beta2_4, beta2_5;
  // gamma, lambda, mu;
  // Nelder mead points contain params, N + 1 = 16 total:
  //PREVIOUS vector<double> P0, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12, \
    P13, P14, P15;
  std::vector<vector<double>> P;

  /* Calc vars */
  
  // unobserved utility
  vector<double> P0_xi0, P1_xi0, P2_xi0, P3_xi0, P4_xi0, P5_xi0, P6_xi0, \
    P7_xi0, P8_xi0, P9_xi0, P10_xi0, P11_xi0, P12_xi0, P13_xi0, P14_xi0, \
    P15_xi0;
  vector<double> P0_xi1, P1_xi1, P2_xi1, P3_xi1, P4_xi1, P5_xi1, P6_xi1, \
    P7_xi1, P8_xi1, P9_xi1, P10_xi1, P11_xi1, P12_xi1, P13_xi1, P14_xi1, \
    P15_xi1;

  // other vars
  vector<double> P0_s_aux1, P1_s_aux1, P2_s_aux1, P3_s_aux1, P4_s_aux1, \
    P5_s_aux1, P6_s_aux1, P7_s_aux1, P8_s_aux1, P9_s_aux1, P10_s_aux1, \
    P11_s_aux1, P12_s_aux1, P13_s_aux1, P14_s_aux1, P15_s_aux1;
  vector<double> P0_s_aux2, P1_s_aux2, P2_s_aux2, P3_s_aux2, P4_s_aux2, \
    P5_s_aux2, P6_s_aux2, P7_s_aux2, P8_s_aux2, P9_s_aux2, P10_s_aux2, \
    P11_s_aux2, P12_s_aux2, P13_s_aux2, P14_s_aux2, P15_s_aux2;
  vector<double> P0_D1, P1_D1, P2_D1, P3_D1, P4_D1, P5_D1, P6_D1, P7_D1, \
    P8_D1, P9_D1, P10_D1, P11_D1, P12_D1, P13_D1, P14_D1, P15_D1;
  vector<double> P0_D2, P1_D2, P2_D2, P3_D2, P4_D2, P5_D2, P6_D2, P7_D2, \
    P8_D2, P9_D2, P10_D2, P11_D2, P12_D2, P13_D2, P14_D2, P15_D2;
  vector<double> P0_s_calc, P1_s_calc, P2_s_calc, P3_s_calc, P4_s_calc, \
    P5_s_calc, P6_s_calc, P7_s_calc, P8_s_calc, P9_s_calc, P10_s_calc, \
    P11_s_calc, P12_s_calc, P13_s_calc, P14_s_calc, P15_s_calc;
  vector<double> P0_ln_s_obs, P1_ln_s_obs, P2_ln_s_obs, P3_ln_s_obs, \
    P4_ln_s_obs, P5_ln_s_obs, P6_ln_s_obs, P7_ln_s_obs, P8_ln_s_obs, \
    P9_ln_s_obs, P10_ln_s_obs, P11_ln_s_obs, P12_ln_s_obs, P13_ln_s_obs, \
    P14_ln_s_obs, P15_ln_s_obs;
};

#endif
