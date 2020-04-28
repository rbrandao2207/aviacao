#ifndef LOGITESTIMATIONHEADERDEF
#define LOGITESTIMATIONHEADERDEF

#include <pqxx/pqxx>
#include <string>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <vector>


class LogitEstimation
{
public:
    LogitEstimation(pqxx::nontransaction& N, std::string const& tblname_, std::vector<int> const& arr_sizes, unsigned long const& max_threads=16);
    void NewtonRaphson(double const& tol, double& step, double& step2, double const& maxiter, std::string const& method="Hninvpar");
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & beta;
        ar & grad;
        ar & LL;
        ar & aux_CS;
    }
    void ConsSurplus();
    void PrintRes(std::string const& res_file);
private:
    unsigned long hardware_threads;
    unsigned long num_threads;
    std::string tblname;
    int numrows;
    int numcoeffs;
    int numprods;
    int iterations;
    double aux_LL;
    double LL;
    double aux_CS;
    double CS;
    boost::numeric::ublas::matrix<double> X;
    boost::numeric::ublas::vector<double> y;
    boost::numeric::ublas::vector<double> aux_P;
    boost::numeric::ublas::vector<double> P;      // choice probabilities
    boost::numeric::ublas::vector<double> aux_beta;
    boost::numeric::ublas::vector<double> beta;
    boost::numeric::ublas::vector<double> aux_grad;
    boost::numeric::ublas::vector<double> grad;
    boost::numeric::ublas::matrix<double> aux_A;
    boost::numeric::ublas::matrix<double> A;
    void grad_calc();
    void gradpar_calc();
    void Hninv_calc();
    void Hninvpar_calc();
    void DFP_calc();
};

#endif
