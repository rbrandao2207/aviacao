#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "BLP.hpp"


BLP::BLP(const std::vector<double> init_guess)
{ 
    beta1_0 = init_guess[0];
    beta1_1 = init_guess[1];
    beta1_2 = init_guess[2];
    beta1_3 = init_guess[3];
    beta1_4 = init_guess[4];
    beta1_5 = init_guess[5];
    beta2_0 = init_guess[6];
    beta2_1 = init_guess[7];
    beta2_2 = init_guess[8];
    beta2_3 = init_guess[9];
    beta2_4 = init_guess[10];
    beta2_5 = init_guess[11];
    gamma = init_guess[12];
    lambda = init_guess[13];
    mu = init_guess[14];
}

void BLP::allocate()
{
    // initialize unobs util and allocate ublas objs
    qsi0.resize(s_obs_wg.size());
    qsi1.resize(s_obs_wg.size());
    std::fill(qsi0.data().begin(), qsi0.data().end(), 0.);
    std::fill(qsi1.data().begin(), qsi1.data().end(), 0.);
    s_aux1.resize(s_obs_wg.size());
    s_aux2.resize(s_obs_wg.size());
    D1.resize(s_obs_wg.size());
    D2.resize(s_obs_wg.size());
    s_calc.resize(s_obs_wg.size());
    ln_s_obs.resize(s_obs_wg.size());
}

void BLP::calc_shares()
{
    //calc s_ind1 & s_ind2 (exp(Xbeta....))
    for (unsigned i = 0; i < s_aux1.size(); ++i) {
        //DEBUG
        if (i == 4807)
            std::cout << X(i,0) << '\t' << X(i,1) << '\t' << X(i,2) << '\t'  << X(i,3) << '\t'  << X(i,4) << '\t'  << X(i,5) << '\t'  << qsi1[i] << std::endl;
        //ENDDEBUG
        s_aux1[i] = std::exp((X(i, 0) * beta1_0 + X(i, 1) * beta1_1 + X(i, 2) * beta1_2 + X(i, 3) * beta1_3 + X(i, 4) * beta1_4 + X(i, 5) * beta1_5 + qsi1[i]) / lambda);
        s_aux2[i] = std::exp((X(i, 0) * beta2_0 + X(i, 1) * beta2_1 + X(i, 2) * beta2_2 + X(i, 3) * beta2_3 + X(i, 4) * beta2_4 + X(i, 5) * beta2_5 + qsi1[i]) / lambda);
    }

    //calc D1 and D2
    double aux_D1 = 0;
    double aux_D2 = 0;
    unsigned initial_aux_i = 0;
    unsigned aux_mkt_id = 0;
    for (unsigned i = 0; i < mkt_id.size(); ++i) {
        if (aux_mkt_id == mkt_id[i]) {
            aux_D1 += s_aux1[i];
            aux_D2 += s_aux2[i];
        } else {
            for (unsigned j = initial_aux_i; j < i; ++j) {
              D1[j] = aux_D1;
              D2[j] = aux_D2;
            }
            aux_D1 = s_aux1[i];
            aux_D2 = s_aux2[i];
            initial_aux_i = i;
            aux_mkt_id = mkt_id[i];
        }
        if (i == mkt_id.size() - 1) {
            for (unsigned j = initial_aux_i; j < i; ++j) {
              D1[j] = aux_D1;
              D2[j] = aux_D2;
            }
        }
    }

    // calc model shares
    for (unsigned i = 0; i < s_calc.size(); ++i) {
        s_calc[i] = gamma * ((s_aux1[i] / D1[i]) * (std::pow(D1[i], lambda) / (1 + std::pow(D1[i], lambda)))) + (1 - gamma) * ((s_aux2[i] / D2[i]) * (std::pow(D2[i], lambda) / (1 + std::pow(D2[i], lambda))));
        //DEBUG
        if (std::isnan(s_calc[i]))
            std::cout << "nan" << std::endl;
        //ENDDEBUG
    }
}

void BLP::contraction(const double contract_tol)
{
    // initialization
    // observe s_obs = s_obs_wg * mu, where s_obs_wg is within group share
    for (unsigned i = 0; i < s_obs_wg.size(); ++i) {
        ln_s_obs[i] = std::max(std::log(s_obs_wg[i] * mu), std::numeric_limits<double>::lowest());
    }
    this->BLP::calc_shares();

    bool conv_check = 0;
    double s_calc_d; //DEBUG
    while (!conv_check) {
        for (unsigned i = 0; i < qsi1.size(); ++i) {
            if (ln_s_obs[i] == std::numeric_limits<double>::lowest()) {
                qsi1[i] = std::numeric_limits<double>::lowest();
            } else {
                qsi1[i] = qsi0[i] + lambda * (ln_s_obs[i] - std::log(s_calc[i]));
            }
            //DEBUG
            if (i == 4807)
                s_calc_d = std::log(s_calc[i]);
            //ENDDEBUG
        }

        // check for convergence
        for (unsigned i = 0; i <= qsi1.size(); ++i) {
            if (i == qsi1.size()) {
                conv_check = 1;
                break;
            }
            if (std::abs(qsi1[i] - qsi0[i]) < contract_tol) {
                continue;
            } else {
                break;
            }
        }

        // update vars
        this->BLP::calc_shares();
        qsi0 = qsi1;
    }
}
