#include <numeric>
#include <pqxx/pqxx>
#include <stdexcept>
#include <string>
#include <vector>

#include "GenMktTbls.hpp"

using namespace std::string_literals;

using std::string;
using std::vector;


GenMktTbls::GenMktTbls(pqxx::nontransaction& N, pqxx::work& W, const string& cidade_orig, const string& cidade_dest, const string& period, bool createtbl)
{
    // filled in airport icao codes 
    if (cidade_orig == "BR") {
        origens = {"SBBR"};
    } else if (cidade_orig == "CT") {
        origens = {"SBCT"};
    } else if (cidade_orig == "FZ") {
        origens = {"SBFZ"};
    } else if (cidade_orig == "PA") {
        origens = {"SBPA"};
    } else if (cidade_orig == "RJ") {
        origens = {"SBRJ", "SBGL"};
    } else if (cidade_orig == "SP") {
        origens = {"SBSP", "SBGR", "SBKP"};
    } else if (cidade_orig == "SV") {
        origens = {"SBSV"};
    } else {
        throw std::runtime_error("Cidade de origem inválida (" + cidade_orig + ")");
    }
    if (cidade_dest == "BR") {
        destinos = {"SBBR"};
    } else if (cidade_dest == "CT") {
        destinos = {"SBCT"};
    } else if (cidade_dest == "FZ") {
        destinos = {"SBFZ"};
    } else if (cidade_dest == "PA") {
        destinos = {"SBPA"};
    } else if (cidade_dest == "RJ") {
        destinos = {"SBRJ", "SBGL"};
    } else if (cidade_dest == "SP") {
        destinos = {"SBSP", "SBGR", "SBKP"};
    } else if (cidade_dest == "SV") {
        destinos = {"SBSV"};
    } else {
        throw std::runtime_error("Cidade de destino inválida (" + cidade_orig + ")");
    }

    // fill produtos
    string aux_ano = period.substr(0, 4);
    string aux_mes = period.substr(5, 7);
    for (auto orig : origens) {
        for (auto dest : destinos) {
            string sql_empresas = "SELECT DISTINCT empresa FROM csv_anac WHERE origem = '" + orig + "' AND destino = '" + dest + "' AND ano = " + aux_ano + " AND mes = " + aux_mes + " AND empresa IN ('GLO', 'AZU', 'TAM', 'ONE');";
            pqxx::result R(N.exec(sql_empresas));
            for (auto c = R.begin(); c != R.end(); ++c) {
                produtos.push_back(orig + "_" + dest + "_" + c[0].as<string>());
            }
        }
    }

    // create tbl
    tblname = cidade_orig + "_" + cidade_dest + "_" + period;
    if (createtbl == true) {
        string sql_createtbl = "DROP TABLE IF EXISTS " + tblname + ";"
            "CREATE TABLE " + tblname + " ("
            "id     bigserial primary key,"
            "obs    int,"
            "alt    text,"
            "dni    int,"
            "preco  real,";
        for (long unsigned int i = 0; i != produtos.size() - 1; ++i) {
            sql_createtbl += "asc_" + std::to_string(i + 1) + "  int,";
        }
        sql_createtbl += "assentos2_4  int,";
        sql_createtbl += "assentos5p   int);";
        W.exec(sql_createtbl);
    }
}

void GenMktTbls::datafill(pqxx::nontransaction& N, pqxx::work& W, const string& period)
{
    string aux_ano = period.substr(0, 4);
    string aux_mes = period.substr(5, 7);
    string aux_origens {};
    for (auto orig : origens) {
        aux_origens += "'" + orig + "',";
    }
    aux_origens.pop_back();
    string aux_destinos {};
    for (auto dest : destinos) {
        aux_destinos += "'" + dest + "',";
    }
    aux_destinos.pop_back();

    // query relevant input
    string sql_obs = "SELECT empresa, origem, destino, tarifa, assentos FROM csv_anac WHERE ano = " + aux_ano + "AND mes = " + aux_mes + " AND origem IN (" + aux_origens + ") AND destino IN (" + aux_destinos + ") AND empresa IN ('GLO', 'AZU', 'TAM', 'ONE');";
    pqxx::result R(N.exec(sql_obs));
    long unsigned int obs_id = 0;
    for (auto c = R.begin(); c != R.end(); ++c, ++obs_id) {
        string sql_insert;
        string dni;
        for (auto produto : produtos) {
            // fill var dni
            if (c[1].as<string>() + "_" + c[2].as<string>() + "_" + c[0].as<string>() == produto) {
                dni = "1";
            } else {
                dni = "0";
            }
            // fill var preco
            string preco = std::to_string(c[3].as<double>());
            // fill vars asc_x
            long unsigned int aux_produto = 0;
            for (auto& produto2 : produtos) {
                if (produto == produto2)
                    break;
                ++aux_produto;
            }
            string produtos_dummies = {};
            for (long unsigned int i = 1; i != produtos.size(); ++i) {
                if (aux_produto == i) {
                    produtos_dummies += "1"s + ","s;
                } else {
                    produtos_dummies += "0"s + ","s;
                }
            }
            produtos_dummies.pop_back();
            // fill vars assentos
            string assentos2_4 = "0";
            string assentos5p = "0";
            if (c[4].as<int>() > 1 && c[4].as<int>() < 5) {
                assentos2_4 = "1";
            } else if (c[4].as<int>() >= 5) {
                assentos5p = "1";
            }
            sql_insert = "INSERT INTO " + tblname + " (obs, alt, dni, preco";
            for (long unsigned int i = 0; i != produtos.size() -1; ++i) {
                sql_insert += ", asc_" + std::to_string(i + 1);
            }
            sql_insert += ", assentos2_4, assentos5p) VALUES (" + std::to_string(obs_id) + ",'" + produto + "'," + dni + "," + preco + "," + produtos_dummies + "," + assentos2_4 + "," + assentos5p + ");";
            W.exec(sql_insert);
        }
    }
}

void GenMktTbls::priceavgs(pqxx::nontransaction& N, pqxx::work& W)
{
    for (auto& produto : produtos) {
        for (int aux_dummy = 0; aux_dummy <= 2; ++aux_dummy) {
            string assentos2_4, assentos5p;
            if (aux_dummy == 0) {
                assentos2_4 = "0", assentos5p = "0";
            } else if (aux_dummy == 1) {
                assentos2_4 = "1", assentos5p = "0";
            } else if (aux_dummy == 2) {
                assentos2_4 = "0", assentos5p = "1";
            }
            // check and delete if num obs equals zero for some assentos dummy
            string sql_count = "SELECT COUNT(*) FROM " + tblname + " WHERE alt = '" + produto + "' AND dni = 1 AND assentos2_4 = " + assentos2_4 + " AND assentos5p = " + assentos5p + ";";
            pqxx::result R0(N.exec(sql_count));
            auto c = R0.begin();
            if (c[0].as<int>() == 0) {
                string sql_delete = "DELETE FROM " + tblname + " WHERE alt = '" + produto + "';";
                W.exec(sql_delete);
                break;
            }

            string sql_query = "SELECT preco FROM " + tblname + " WHERE alt = '" + produto + "' AND dni = 1 AND assentos2_4 = " + assentos2_4 + " AND assentos5p = " + assentos5p + ";";
            vector<double> prices;
            pqxx::result R(N.exec(sql_query));
            for (auto c = R.begin(); c != R.end(); ++c) {
                prices.push_back(c[0].as<double>());
            }
            double avgprice;
            if (prices.size() > 0) {
                avgprice = std::accumulate(prices.begin(), prices.end(), 0.0) / prices.size();
            } else {
                throw std::runtime_error("Found 0 prices in table " + tblname + " case aux_dummy = " + std::to_string(aux_dummy) + " and produto = " + produto + " (check GenMktTbls::priceavgs)");
            }
            string sql_update = "UPDATE " + tblname + " SET preco = " + std::to_string(avgprice) + " WHERE alt = '" + produto + "' AND dni = 0 AND assentos2_4 = " + assentos2_4 + " AND assentos5p = " + assentos5p + ";";
            W.exec(sql_update);
        }
    }
}
