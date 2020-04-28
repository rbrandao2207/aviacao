#ifndef ESTIMATIONHEADERDEF
#define ESTIMATIONHEADERDEF

#include <iostream>
#include <pqxx/pqxx>
#include <string>
#include <vector>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>


class Estimation
{
public:
    Estimation(const std::string start_date_, const std::string end_date_, const unsigned periodicity_, const unsigned int num_threads);
    ~Estimation()
    {
    }
private:
    std::string start_date;
    std::string end_date;
    std::string periodicity;
    boost::numeric::ublas::vector<double> alpha;
    boost::numeric::ublas::vector<double> beta;
    boost::numeric::ublas::vector<double> gamma;
    double lambda;
    boost::numeric::ublas::matrix<double> obs_shares;
    boost::numeric::ublas::vector<double> partial_calc_shares; // for paralelism
    boost::numeric::ublas::matrix<double> calc_shares;
    boost::numeric::ublas::matrix<double> error_term;
};

#endif
