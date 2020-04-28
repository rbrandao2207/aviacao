#include <iostream>
#include <pqxx/pqxx>
#include <string>
#include <vector>

#include "dbfunctions.hpp"
#include "GenMktTbls.hpp"

using std::string;
using std::vector;


int main(int argc, char* argv[])
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

        // csv import
        if (argc > 1 && std::strcmp(argv[1], "csvimport") == 0) {
            pqxx::work W0(C);
            //import anac
            vector<string> anos = {"2014", "2015", "2016", "2017", "2018", "2019"};
            vector<string> meses = {"01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"};
            string datadir = "/var/lib/postgresql/data/pgdata/databases/anac/";
            for (string ano : anos) {
                for (string mes : meses) {
                    if (ano == "2019" && !(mes == "01" || mes == "02"))
                        continue;
                    if (ano == "2014" && mes == "01") {
                        csvimport_anac(W0, datadir, ano, mes, true);  // create tbl if first
                    } else {
                        csvimport_anac(W0, datadir, ano, mes);
                    }
                }
            }
            //import abear
            //datadir = "/var/lib/postgresql/data/pgdata/databases/abear/";
            //csvimport_abear(W0, datadir);
            W0.commit();
        }

        // build var tables, e.g. SP_RJ_yymm
        //vector<string> cidades_orig = {"RJ", "SP", "BR", "CT", "PA", "FZ", "SV"};
        vector<string> cidades_orig = {"CT", "PA", "FZ", "SV"};
        vector<string> cidades_dest = {"RJ", "SP", "BR", "CT", "PA", "FZ", "SV"};
        vector<string> anos = {"2014", "2015", "2016", "2017", "2018", "2019"};
        vector<string> meses = {"01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"};
        for (auto cidade_orig : cidades_orig) {
            for (auto cidade_dest : cidades_dest) {
                if (cidade_orig == cidade_dest)
                    continue;
                if ((cidade_orig == "CT" && cidade_dest == "SP") || (cidade_orig == "CT" && cidade_dest == "RJ") || (cidade_orig == "CT" && cidade_dest == "SV") || (cidade_orig == "CT" && cidade_dest == "FZ") || (cidade_orig == "PA" && cidade_dest == "SV") || (cidade_orig == "PA" && cidade_dest == "FZ"))
                    continue;
                for (auto ano : anos) {
                    for (auto mes: meses) {
                        if ((cidade_orig == "RJ" && cidade_dest == "SP") && ((ano == "2014") || ((ano == "2015") && (!(mes == "09" || mes == "10" || mes == "11" || mes == "12")))))
                            continue;
                        if (ano == "2019" && !(mes == "01" || mes == "02"))
                            continue;
                        pqxx::nontransaction N1(C);
                        pqxx::work W1(C2);
                        string period = ano + "_" + mes;
                        GenMktTbls mkt_instance(N1, W1, cidade_orig, cidade_dest, period);
                        //turn off comment below and on further down to skip tbl creation
                        //GenMktTbls mkt_instance(N1, W1, cidade_orig, cidade_dest, period, false);
                        W1.commit();
                        /**/
                        pqxx::work W2(C2);
                        mkt_instance.datafill(N1, W2, period);
                        W2.commit();
                        std::cout << "Dados " << cidade_orig << cidade_dest << ano << mes << " carregados c/ sucesso" << std::endl;
                        /**/
                        pqxx::work W3(C2);
                        mkt_instance.priceavgs(N1, W3);
                        W3.commit();
                        std::cout << "Preços médios " << cidade_orig << cidade_dest << ano << mes << " imputados c/ sucesso" << std::endl;
                    }
                }
            }
        }

        C.disconnect();
        C2.disconnect();
        std::cout << "Database operation successfull" << std::endl;
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}
