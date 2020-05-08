#ifndef BLPHEADERDEF 
#define BLPHEADERDEF

#include <string>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>


class BLP
{
public:
    BLP(const std::vector<double> init_guess);
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & s_obs_wg;
        ar & mkt_id;
        ar & pop_ave;
        ar & X;
    }
    void allocate();
    void calc_shares();
    void contraction(const double contract_tol);
private:

    // params
    double beta1_0, beta1_1, beta1_2, beta1_3, beta1_4, beta1_5;
    double beta2_0, beta2_1, beta2_2, beta2_3, beta2_4, beta2_5;
    double gamma, lambda, mu;

    // exogenous vars
    boost::numeric::ublas::vector<double> s_obs_wg;
    boost::numeric::ublas::vector<double> mkt_id;
    boost::numeric::ublas::vector<double> pop_ave;
    /* X matrix structure: 0 constant
                           1 fare (bin center R$/1000)
                           2 distance (km / 1000)
                           3 distance^2
                           4 carrier dummy
                           5 time dummy
    */
    boost::numeric::ublas::matrix<double> X;

    // unobserved utility
    boost::numeric::ublas::vector<double> qsi0;
    boost::numeric::ublas::vector<double> qsi1;

    // calculated vars (params dependent)
    boost::numeric::ublas::vector<double> s_aux1;
    boost::numeric::ublas::vector<double> s_aux2;
    boost::numeric::ublas::vector<double> D1;
    boost::numeric::ublas::vector<double> D2;
    boost::numeric::ublas::vector<double> s_calc;
    boost::numeric::ublas::vector<double> ln_s_obs;
};

#endif
