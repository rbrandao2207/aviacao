#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <pqxx/pqxx>
#include <string>
#include <tuple>
#include <vector>

#include "Estimation.hpp"
#include "AuxFunctions.hpp"


// generates postgres tables mkts_yymm and prods_yymm 
int genTbls(unsigned qty_threshold, std::vector<double> bins, std::vector<std::string> thread_tbls)
{
    try {
        pqxx::connection C("dbname = aviacao user = postgres password = passwd hostaddr = 127.0.0.1 port = 5432");
        pqxx::connection C2("dbname = aviacao user = postgres password = passwd hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open() && C2.is_open()) {
            std::cout << "Opened database successfully: " << C.dbname() << std::endl;
        } else {
            std::cout << "Can't open database" << std::endl;
            return 1;
        }
        for (auto& tbl : thread_tbls) {

            // create tbls
            pqxx::work W(C);
            std::string tblname1 = "mkts_" + tbl;
            std::string sql_createtbl1 = "DROP TABLE IF EXISTS " + tblname1 + "; "
                "CREATE TABLE " + tblname1 + " ("
                "id         int, "
                "name       text, "
                "min        real, "
                "max        real, "
                "mean       real, "
                "std_dev    real, "
                "tot_qty    int);";
            W.exec(sql_createtbl1);
            std::string tblname2 = "prodQties_" + tbl;
            std::string sql_createtbl2 = "DROP TABLE IF EXISTS " + tblname2 + "; "
                "CREATE TABLE " + tblname2 + " ("
                "mkt_id     int, "
                "empresa    text";
            for (unsigned bin_nbr = 0; bin_nbr < bins.size(); ++bin_nbr)
                sql_createtbl2 += ", bin" + std::to_string(bin_nbr) + " int";
            sql_createtbl2 += ");";
            W.exec(sql_createtbl2);
            W.commit();

            // get data from csvs
            pqxx::nontransaction N(C);
            std::string sql_getmkts = "SELECT DISTINCT origem, destino FROM csv_" + tbl + ";";
            pqxx::result R(N.exec(sql_getmkts));

            // init aux var and connection for data insertion
            unsigned mkt_id = 0;
            pqxx::work W2(C2);

            // loop mkts
            for (auto mkt = R.begin(); mkt != R.end(); ++mkt) {

                std::string sql_gettarifs = "SELECT empresa, tarifa, assentos FROM csv_" + tbl + " WHERE origem = '" + mkt[0].as<std::string>() + "' AND destino = '" + mkt[1].as<std::string>() + "';";
                pqxx::result R2(N.exec(sql_gettarifs));
                unsigned tot_qty = 0;
                double tot_revenue = 0;
                std::vector<std::tuple<std::string, double, unsigned>> x;
                for (auto mkt_obs = R2.begin(); mkt_obs != R2.end(); ++mkt_obs) {
                    tot_qty += mkt_obs[2].as<unsigned>();
                    tot_revenue += mkt_obs[1].as<double>() * mkt_obs[2].as<double>();
                    x.push_back(std::make_tuple(mkt_obs[0].as<std::string>(), mkt_obs[1].as<double>(), mkt_obs[2].as<int>()));
                }

                // process and fill data for mkts
                // (qty_threshold param allows for a priori exclusion of thin markets)
                if (tot_qty >= qty_threshold) {

                    // calc mean
                    double mean = tot_revenue / tot_qty;

                    // sort and calc other stats
                    double min;
                    double max;
                    unsigned cum_qty = 0;
                    double std_dev;
                    long double cum_std_dev = 0;
                    //
                    std::sort(x.begin(), x.end(), [](std::tuple<std::string, double, unsigned> const &t1, std::tuple<std::string, double, unsigned> const &t2) {
                            return std::get<1>(t1) < std::get<1>(t2); });
                    min = std::get<1>(x[0]);
                    max = std::get<1>(x.back());
                    for (long unsigned i = 0; i < x.size(); ++i) {
                        cum_std_dev += pow(std::get<1>(x[i]) - mean, 2) * std::get<2>(x[i]) / (double) tot_qty;
                        cum_qty += std::get<2>(x[i]);
                    }
                    std_dev = pow(cum_std_dev, 0.5);

                    // Insert market data in database
                    // TODO last insert field below wrong duplicate (std_dev); check if tot_qty or cum_qty applies and substitute
                    std::string insert_mkt = "INSERT INTO " + tblname1 + " VALUES (" + std::to_string(mkt_id) + ", '" + mkt[0].as<std::string>() + "_" + mkt[1].as<std::string>() + "', " + std::to_string(min) + ", " + std::to_string(max) + ", " + std::to_string(mean) + ", " + std::to_string(std_dev) + ", " + std::to_string(std_dev) + ");";
                    W2.exec(insert_mkt);
                    ++mkt_id;


                    // Insert prodQties **** implementar tirando subvector relativo à cada empresa
                    // fazer sort de x em get<0>, usar var aux_empresa, loopar x dando pushback em vetor de empresas se get<0> != aux_empresa
                    // no loop acima ir somando bins
                }
            }
            W2.commit();
        }
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}
