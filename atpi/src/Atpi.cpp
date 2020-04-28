#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <pqxx/pqxx>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "Atpi.hpp"

namespace ublas = boost::numeric::ublas;


Atpi::Atpi(const std::string& prods_file, const std::string& tbls_file, const std::string& nperiods_file, const unsigned int& inst_nbr, const unsigned int& num_threads, const std::string& data_dir_)
{
    // fill products
    std::ifstream fdesc_prods;
    std::string product;
    fdesc_prods.open(prods_file);
    assert(fdesc_prods.is_open());
    while (std::getline(fdesc_prods, product))
        products.push_back(product);
    fdesc_prods.close();

    // get nperiods
    std::ifstream fdesc_nperiods;
    fdesc_nperiods.open(nperiods_file);
    assert(fdesc_nperiods.is_open());
    std::string aux_nperiods;
    std::getline(fdesc_nperiods, aux_nperiods);
    fdesc_nperiods.close();
    unsigned int nperiods = std::stoi(aux_nperiods);

    // find out indexes for instance tables
    unsigned int quotient = nperiods / num_threads;
    unsigned int remainder = nperiods % num_threads;
    unsigned int count_quo = 0;
    unsigned int count_rem = 0;
    unsigned int count_inst = 0;
    unsigned int count_line = 0;

    // fill tables
    std::ifstream fdesc_tbls;
    std::ifstream fdesc2_tbls;
    std::string table;
    std::string table2;
    fdesc_tbls.open(tbls_file);
    assert(fdesc_tbls.is_open());
    if (quotient == 0) {
        while (std::getline(fdesc_tbls, table)) {
            tables.push_back(table);
        }
    } else {
        while (count_line < nperiods) {
            for (count_quo = 0; count_quo < quotient; ++count_quo) {
                if (count_inst == inst_nbr) { 
                    std::getline(fdesc_tbls, table);
                    tables.push_back(table);
                    ++count_line;
                } else {
                    std::getline(fdesc_tbls, table);
                    ++count_line;
                }
            }
            if (count_inst < remainder) {
                ++count_rem;
                if (count_inst == inst_nbr) {
                    std::getline(fdesc_tbls, table);
                    tables.push_back(table);
                    ++count_line;
                } else {
                    std::getline(fdesc_tbls, table);
                    ++count_line;
                }
            }
            ++count_inst;
        }
    }
    fdesc_tbls.close();
    if (inst_nbr < num_threads - 1) {
        fdesc2_tbls.open(tbls_file);
        assert(fdesc2_tbls.is_open());
        while (std::getline(fdesc2_tbls, table2)) {
            if (table2 == tables.back()) {
                std::getline(fdesc2_tbls, table2);
                tables.push_back(table2);
                fdesc2_tbls.close();
                break;
            }
        }
    }
    fdesc2_tbls.close();

    // fill vars
    data_dir = data_dir_;
    persist_file = data_dir_ + "/" + tables[0] + "-" + tables.back();
}

int Atpi::FillMats(unsigned qty_threshold)
{
    // allocate memory
    uvals.resize(products.size(), tables.size());
    qties.resize(products.size(), tables.size());

    // connect to database
    try {
        pqxx::connection C("dbname = aviacao user = postgres password = passwd hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            std::cout << "Opened database successfully: " << C.dbname() << std::endl;
        } else {
            std::cout << "Can't open database" << std::endl;
            return 1;
        }
        pqxx::nontransaction N(C);
        // perform computations
        for (unsigned i = 0; i < uvals.size1(); ++i) {
            std::string product = products[i];
            std::vector<std::string> split_prod;
            boost::split(split_prod, product, [](char c){return c == ',';});
            for (unsigned j = 0; j < uvals.size2(); ++j) {
                std::string table = tables[j];
                std::string query = "SELECT tarifa, assentos FROM " + table + " WHERE origem = '" + split_prod[0] + "' AND destino = '" + split_prod[1] + "' AND empresa = '" + split_prod[2] + "';";
                pqxx::result R(N.exec(query));
                if (R.size() == 0) {
                    uvals(i, j) = 0.0;
                    qties(i, j) = 0;
                } else {
                    std::vector<double> receitas;
                    std::vector<int> assentos;
                    for (auto c = R.begin(); c != R.end(); ++c) {
                        receitas.push_back(c[0].as<double>() * c[1].as<int>());
                        assentos.push_back(c[1].as<int>());
                    }
                    int qty;
                    qty = std::accumulate(assentos.begin(), assentos.end(), 0);
                    qties(i, j) = qty;
                    double uval;
                    uval = std::accumulate(receitas.begin(), receitas.end(), 0.0) / qty;
                    uvals(i, j) = uval;
                }
                // time check
                if ((i == (uvals.size1() / 5)) && (j == uvals.size2() - 1))
                    std::cout << "finished one fifth of uvals and qties for instance ending on table " << table << std::endl;
            }
        }
        std::cout << "finished 100 percent of uvals and qties for instance ending on table " << tables.back() << std::endl;
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    // Calc indicator matrix c12:
    // traverse qties; if > 0 and if col-1 > 0 store +1
    // else store -1
    c12.resize(products.size(), tables.size() - 1);
    for (unsigned col = 1; col < qties.size2(); ++col) {
        for (unsigned row = 0; row < qties.size1(); ++row) {
            if ((qties(row, col) > qty_threshold) && (qties(row, col - 1) > qty_threshold)) {
                c12(row, col - 1) = 1;
            } else {
                c12(row, col - 1) = -1;
            }
        }
    }

    // Calc weights
    // allocate memory
    weights.resize(products.size(), tables.size());

    // calc
    weights = element_prod(uvals, qties);
    for (unsigned col = 1; col < weights.size2(); ++col) {
        long double column_sum = 0;
        for (unsigned row1 = 0; row1 < weights.size1(); ++row1) {
            if (c12(row1, col - 1) > 0)
                column_sum += weights(row1, col);
        }
        for (unsigned row2 = 0; row2 < weights.size1(); ++row2) {
            if (c12(row2, col - 1) > 0) {
                weights(row2, col) /= column_sum;
                assert(("Weights check failed, increase qty_threshold", weights(row2, col) * weights(row2, col) > std::numeric_limits<long double>::min()));
            } else {
                weights(row2, col) = 0;
            }
        }
    }
    // calc weights for first col
    unsigned col = 0;
    for (unsigned row = 0; row < weights.size1(); ++row) {
        if (c12(row, col) > 0) {
            weights(row, col) = weights(row, col + 1);
        } else {
            weights(row, col) = 0;
        }
    }
    std::cout << "finished weights for tables ending on " << tables.back() << std::endl;

    // Calc uidxs
    uidxs.resize(products.size(), tables.size() - 1);
    for (unsigned col = 0; col < uidxs.size2(); ++col) {
        for (unsigned row = 0; row < uidxs.size1(); ++row) {
            if (c12(row, col) > 0) {
                uidxs(row, col) = uvals(row, col + 1) / uvals(row, col);
            }
        }
    }
    std::cout << "finished uidxs for tables ending on " << tables.back() << std::endl;

    return 0;
}

int Atpi::CalcIdxs()
{
    unsigned vsize = uidxs.size2();
    std::vector<long double> laspeyres(vsize);
    std::vector<long double> paasche(vsize);
    std::vector<long double> fisher(vsize);
    std::vector<long double> jevons(vsize);
    std::vector<long double> tornqvist(vsize);
    std::vector<long double> walsh(vsize);

    for (unsigned col = 0; col < uidxs.size2(); ++col) {
        long double aux_laspeyres = 0;
        long double aux_paasche = 0;
        long double aux_jevons = 1;
        long double aux_tornqvist = 1;
        long double aux1_walsh = 0;
        long double aux2_walsh = 0;
        for (unsigned row = 0; row < uidxs.size1(); ++row) {
            if (weights(row, col + 1) > 0) {
                aux_laspeyres += weights(row, col) * uidxs(row, col);
                aux_paasche += weights(row, col + 1) / uidxs(row, col);
                aux_jevons *= pow(uidxs(row, col), weights(row, col));
                aux_tornqvist *= pow(uidxs(row, col), (weights(row, col) + weights(row, col + 1)) / 2);
                long double aux_walsh_wgts = pow(weights(row, col) * weights(row, col + 1), 0.5);
                aux1_walsh += aux_walsh_wgts * pow(uidxs(row, col), 0.5);
                aux2_walsh += aux_walsh_wgts * pow(1 / uidxs(row, col), 0.5);
            }
        }
        laspeyres[col] = aux_laspeyres; 
        paasche[col] = 1 / aux_paasche;
        long double aux_fisher = laspeyres[col] * paasche[col];
        fisher[col] = pow(aux_fisher, 0.5);
        jevons[col] = aux_jevons;
        tornqvist[col] = aux_tornqvist;
        walsh[col] = aux1_walsh / aux2_walsh;
    }

    std::string persistf = data_dir + "/idxs_" + tables[0].substr(4) + "-" + tables.back().substr(4);
    std::ofstream fdesc; 
    fdesc.open(persistf);
    assert(fdesc.is_open());
    for (unsigned it = 0; it < uidxs.size2(); ++it)
        fdesc << laspeyres[it] << '\t' << paasche[it] << '\t' << fisher[it] << '\t' << jevons[it] << '\t' << tornqvist[it] << '\t' << walsh[it] << '\t' << '\n';
    fdesc.close();
    std::cout << "finished idxs for tables ending on " << tables.back() << std::endl;

    return 0;
}
