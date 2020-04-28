#ifndef AUXFUNCTIONSHEADERDEF
#define AUXFUNCTIONSHEADERDEF

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cmath>
#include <pqxx/pqxx>
#include <string>
#include <vector>


std::vector<int> getArraySizes(pqxx::nontransaction& N, std::string tblname);

std::string getSyscallOutput(const char* cmd);

double Exp(double x);

struct compare_gt
{
    double tol;
    compare_gt(double const& tol_): tol(tol_) { }

    bool operator()(double const& i)
    {
        return (std::abs(i) > tol);
    }
};

struct P_calc
{
    boost::numeric::ublas::vector<double>& P;
    boost::numeric::ublas::vector<double>& aux_P;
    boost::numeric::ublas::vector<double>& y;
    double aux_sum;
    double& LL_partial;
    P_calc(boost::numeric::ublas::vector<double>& P_, boost::numeric::ublas::vector<double>& aux_P_, boost::numeric::ublas::vector<double>& y_, double& LL_partial_) : P(P_), aux_P(aux_P_), y(y_), LL_partial(LL_partial_) {}
    void operator() (long unsigned int const& start_row, long unsigned int const& end_row, int const& numprods)
    {
        LL_partial = 0;
        for (long unsigned int n1 = start_row; n1 < end_row; n1 += numprods) {
            aux_sum = 0.0;
            for (int i = 0; i < numprods; ++i) {
                aux_sum += aux_P(n1 + i);
            }
            for (int i = 0; i < numprods; ++i) {
                P(n1 + i) /= aux_sum;
                LL_partial += y(n1 + i) * std::log(P(n1 + i));
            }
        }

    }
};

struct Hn_calc
{
    const int& numrows;
    const int& numcoeffs;
    boost::numeric::ublas::matrix<double> const& X;
    boost::numeric::ublas::vector<double> const& P;
    boost::numeric::ublas::matrix<double>& aux_A;
    Hn_calc(int const& numrows_, int const& numcoeffs_, boost::numeric::ublas::matrix<double> const& X_, boost::numeric::ublas::vector<double> const& P_, boost::numeric::ublas::matrix<double>& aux_A_) : numrows(numrows_), numcoeffs(numcoeffs_), X(X_), P(P_), aux_A(aux_A_) {}
    void operator() (int const& start_row, int const& end_row)
    {
        for (int j = start_row; j < end_row; ++j) {
            for (int k = 0; k < numcoeffs; ++k) {
                if (j == k) {
                    for (int i = 0; i < numrows; ++i) {
                        double xi_dotprod = 0;
                        for (int l = 0; l < numcoeffs; ++l)
                            xi_dotprod += X(i, l) * X(i, l);
                        aux_A(j, k) += -1 * P(i) * (1 - P(i)) * xi_dotprod;
                    }
                } else {
                    for (int i = 0; i < numrows; ++i) {
                        double xi_dotprod = 0;
                        for (int l = 0; l < numcoeffs; ++l)
                            xi_dotprod += X(i, l) * X(i, l);
                        aux_A(j, k) += P(i) * P(k) * xi_dotprod;
                    }
                }
            }
        }
    }
};

#endif
