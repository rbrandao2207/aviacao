#include <algorithm>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <pqxx/pqxx>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "auxFunctions.hpp"
#include "LogitEstimation.hpp"

namespace ublas = boost::numeric::ublas;


LogitEstimation::LogitEstimation(pqxx::nontransaction& N, std::string const& tblname_, std::vector<int> const& arr_sizes, unsigned long const& max_threads)
{
    // initialize y and X
    tblname = tblname_;
    std::string sql_query = "SELECT dni, preco";
    for (int i = 0; i != arr_sizes[1] - 1; ++i) {
        sql_query += ", asc_" + std::to_string(i + 1);
    }
    sql_query += " FROM " + tblname + " ORDER BY id;";
    pqxx::result R(N.exec(sql_query));
    int row_n = 0;
    X.resize(arr_sizes[0], arr_sizes[1]);
    y.resize(arr_sizes[0]);
    for (auto c = R.begin(); c != R.end(); ++c ) {
        y(row_n) = c[0].as<double>();
        for (pqxx::row::size_type i = 1; i != c.size(); ++i) {
            if (i == 1) {
                X(row_n, i - 1) = c[i].as<double>() / 1000;
            } else {
                X(row_n, i - 1) = c[i].as<double>();
            }
        }
        ++row_n;
    }

    // allocate P, aux_P
    P.resize(arr_sizes[0]);
    aux_P.resize(arr_sizes[0]);

    // initialize beta, allocat aux_beta
    beta.resize(arr_sizes[1]);
    beta(0) = .1;
    for (auto i = 1; i != arr_sizes[1]; ++i)
        beta(i) = .1;
    aux_beta.resize(arr_sizes[1]);

    // allocate grad
    grad.resize(arr_sizes[1]);

    // allocate A and aux_A
    A.resize(arr_sizes[1], arr_sizes[1]);
    aux_A.resize(arr_sizes[1], arr_sizes[1]);

    // init numrows, numcoeffs and numprods
    numrows = arr_sizes[0];
    numcoeffs = arr_sizes[1];
    numprods = arr_sizes[2];

    // init parallel params
    hardware_threads = std::thread::hardware_concurrency();
    num_threads = std::min(hardware_threads != 0 ? hardware_threads : 2, max_threads);
}

void LogitEstimation::grad_calc()
{
    aux_P = ublas::prod(X, beta);
    std::transform(aux_P.begin(), aux_P.end(), aux_P.begin(), Exp);
    P = aux_P;
    LL = 0;
    aux_CS = 0;
    double aux_sum;
    // calc Pni for every n and every i (n -> obs, i -> alt)
    for (long unsigned int n1 = 0; n1 < P.size(); n1 += numprods) {
        aux_sum = 0.0;
        for (int i = 0; i < numprods; ++i) {
            aux_sum += aux_P(n1 + i);
        }
        for (int i = 0; i < numprods; ++i) {
            P(n1 + i) /= aux_sum;
            LL += y(n1 + i) * std::log(P(n1 + i));
            aux_CS += std::log(aux_sum);
        }
    }
    grad = ublas::trans(ublas::prod(ublas::trans(y - P), X));
}

void LogitEstimation::gradpar_calc()
{
    aux_P = ublas::prod(X, beta);
    std::transform(aux_P.begin(), aux_P.end(), aux_P.begin(), Exp);
    P = aux_P;
    std::vector<std::thread> threads(num_threads - 1);
    std::vector<double> LL_partials;
    std::vector<P_calc> p_calcs;
    unsigned long const block_size = numrows / numprods / num_threads;
    int j = 0, k;
    for (unsigned long i = 0; i < (num_threads - 1); ++i)
    {
        k = j + (block_size * numprods);
        LL_partials.push_back(0);
        p_calcs.push_back(P_calc(P, aux_P, y, std::ref(LL_partials[i])));
        threads[i] = std::thread(p_calcs[i], j, k, numprods);
        j = k;
    }
    double LL_partial;
    P_calc p_calc(P, aux_P, y, std::ref(LL_partial));
    p_calc(j, numrows, numprods);
    LL = std::accumulate(LL_partials.begin(), LL_partials.end(), LL_partial);
    for (auto& entry : threads)
        entry.join();
    grad = ublas::trans(ublas::prod(ublas::trans(y - P), X));
}

void LogitEstimation::Hninv_calc()
{
    std::fill(aux_A.data().begin(), aux_A.data().end(), 0.0);
    Hn_calc hn_calc(numrows, numcoeffs, X, P, aux_A);
    hn_calc(0, numcoeffs);
    A = ublas::identity_matrix<double> (aux_A.size1());
    ublas::permutation_matrix<size_t> pm(aux_A.size1());
    ublas::lu_factorize(aux_A, pm);
    ublas::lu_substitute(aux_A, pm, A);
    A *= -1;
}

void LogitEstimation::Hninvpar_calc()
{
    unsigned long const block_size = numcoeffs / num_threads;
    std::vector<std::thread> threads(num_threads - 1);
    std::fill(aux_A.data().begin(), aux_A.data().end(), 0.0);
    Hn_calc hn_calc(numrows, numcoeffs, X, P, aux_A);
    int j = 0, k;
    for (unsigned long i = 0; i < (num_threads - 1); ++i)
    {
        k = j + block_size;
        threads[i] = std::thread(hn_calc, j, k);
        j = k;
    }
    hn_calc(j, numcoeffs);
    for (auto& entry : threads)
        entry.join();

    A = ublas::identity_matrix<double> (aux_A.size1());
    ublas::permutation_matrix<size_t> pm(aux_A.size1());
    ublas::lu_factorize(aux_A, pm);
    ublas::lu_substitute(aux_A, pm, A);
    A *= -1;
}

void LogitEstimation::DFP_calc()
{
    // implement
}

void LogitEstimation::NewtonRaphson(double const& tol, double& step, double& step2, double const& maxiter, std::string const& method)
{
    iterations = 0;
    int iterations2 = 0;
    grad_calc();
    A = ublas::identity_matrix<double> (numcoeffs);
    while (true) {
        if (std::none_of(grad.begin(), grad.end(), compare_gt(tol))) {
            std::cout << "Convergence attained\n";
            break;
        }
        if (method == "steep") {
            ;
        } else if (method == "Hninv") {
            Hninv_calc();
        } else if (method == "Hninvpar") {
            Hninvpar_calc();
        } else if (method == "DFP") {
            aux_A = A;
            DFP_calc();
        } else {
            throw std::runtime_error("NR curvature method not available");
        }
        aux_beta = beta;
        aux_grad = grad;
        aux_LL = LL;
        if (method != "steep") {
            beta = aux_beta + step * ublas::prod(A, grad);
        } else {
            beta = aux_beta + step * grad;
        }
        grad_calc();
        /*// adjust step size
        double step2 = step;
        if (LL > aux_LL) {
            while (LL > aux_LL) {
                step2 *= 2;
                aux_LL = LL;
                beta += step2 * ublas::prod(A, grad);
                grad_calc();
            }
        } else if (LL < aux_LL) {
            while (LL < aux_LL) {
                step2 /= 2;
                aux_LL = LL;
                beta += step2 * ublas::prod(A, grad);
                grad_calc();
            }
        }
        */
        ++iterations;
        std::cout << "Iteration: " << iterations << std::endl;
        std::cout << "grad: " << grad << std::endl;
        std::cout << "beta: " << beta << std::endl;
        std::cout << "LL: " << LL << std::endl;

        /*// check maxiter and switch method
        ++iterations2;
        if (iterations2 == maxiter) {
            for (auto i = 1; i != numcoeffs; ++i)
                beta(i) = .1;
            step = step2;
            while (true) {
                if (std::none_of(grad.begin(), grad.end(), compare_gt(tol))) {
                    std::cout << "Convergence attained\n";
                    break;
                }
                Hninv_calc();
                aux_beta = beta;
                aux_grad = grad;
                aux_LL = LL;
                beta = aux_beta + step * ublas::prod(A, grad);
                grad_calc();
                ++iterations;
                std::cout << "Iteration: " << iterations << std::endl;
                std::cout << "grad: " << grad << std::endl;
                std::cout << "beta: " << beta << std::endl;
                std::cout << "LL: " << LL << std::endl;
            }
        }
        */

        /*// check for deadlock and tremble
        ++iterations2;
        if (iterations2 >= maxiter) {
            typedef std::mt19937 MyRNG;
            MyRNG rng;
            uint32_t seed_val = 0;
            rng.seed(seed_val);
            std::uniform_int_distribution<uint32_t> uint_dist(0, 10);
            for (auto i = 0; i != numcoeffs; ++i) {
                double rand_disturb = uint_dist(rng) / 1e3;
                if (i == 0) {
                    beta(i) += rand_disturb;
                } else {
                    beta(i) = rand_disturb;
                }
            }
            iterations2 = 0;
        }
        */
        /*// check for deadlock and tremble (first version)
        aux_beta = beta - aux_beta;
        if (std::none_of(aux_beta.begin(), aux_beta.end(), compare_gt(tol*1000))) {
            typedef std::mt19937 MyRNG;
            MyRNG rng;
            uint32_t seed_val = 0;
            rng.seed(seed_val);
            std::uniform_int_distribution<uint32_t> uint_dist(0, 10);
            for (auto i = 0; i != numcoeffs; ++i) {
                double rand_disturb = uint_dist(rng) / 1e3;
                std::cout << rand_disturb;
                if (i == 0) 
                beta(i) += rand_disturb;
            }
        }
        */
    }
}

void LogitEstimation::ConsSurplus()
{
    CS = -1 * aux_CS / beta(0);
}

void LogitEstimation::PrintRes(std::string const& res_file)
{
    std::ofstream file_dsc;
    file_dsc.open(res_file);
    assert(file_dsc.is_open());
    file_dsc << "iterations: " << iterations << '\n';
    file_dsc << "beta: " << beta << '\n';
    file_dsc << "grad: " << grad << '\n';
    file_dsc << "LL: " << LL << '\n';
    file_dsc << "CS: " << CS << '\n';
}
