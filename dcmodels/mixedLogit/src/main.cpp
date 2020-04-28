#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <fstream>
#include <iostream>
#include <pqxx/pqxx>
#include <string>
#include <unistd.h>
#include <vector>

#include "auxFunctions.hpp"
#include "mixedLogitEstimation.hpp"


int main(int argc, char* argv[])
{
    double const tol = 0.01;
    double step = .0005; // default = .0005
    double step2 = .05;
    double maxiter = 5000;
    std::string datadir = "results/";

    try {
        pqxx::connection C("dbname = aviacao user = postgres password = passwd hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            std::cout << "Opened database successfully: " << C.dbname() << std::endl;
        } else {
            std::cout << "Can't open database" << std::endl;
            return 1;
        }
        //std::vector<std::string> cidades_orig = {"CT", "PA", "FZ", "SV", "RJ", "BR", "SP"};
        //std::vector<std::string> cidades_dest = {"CT", "PA", "FZ", "SV", "RJ", "BR", "SP"};
        std::vector<std::string> cidades_orig = {"SP", "BR", "RJ"};
        std::vector<std::string> cidades_dest = {"RJ", "BR", "SP"};
        std::vector<std::string> anos = {"2014", "2015", "2016", "2017", "2018", "2019"};
        std::vector<std::string> meses = {"01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"};
        //std::vector<std::string> anos = {"2018"};
        //std::vector<std::string> meses = {"05"};
        for (auto cidade_orig : cidades_orig) {
            for (auto cidade_dest : cidades_dest) {
                if (cidade_orig == cidade_dest)
                    continue;
                if (cidade_orig == "SP" && !(cidade_dest == "RJ"))
                    continue;
                for (auto& ano : anos) {
                    for (auto& mes : meses) {
                        if (ano == "2019" && !(mes == "01" || mes == "02"))
                            continue;
                        pqxx::nontransaction N(C);
                        std::string tblname = cidade_orig + "_" + cidade_dest + "_" + ano + "_" + mes;
                        std::string le_file = datadir + "LE/LE_" + tblname;
                        std::vector<int> arr_sizes = getArraySizes(N, tblname);
                        LogitEstimation inst_LE(N, tblname, arr_sizes);

                        // Estimation
                        if (argc > 1 && std::strcmp(argv[1], "estimation") == 0) {
                            inst_LE.NewtonRaphson(tol, step, step2, maxiter, "steep");
                            std::ofstream ofs(le_file);
                            {
                                boost::archive::text_oarchive oa(ofs);
                                oa << inst_LE;
                            }
                        } else {
                            std::ifstream ifs(le_file);
                            assert(ifs.is_open());
                            boost::archive::text_iarchive ia(ifs);
                            ia >> inst_LE;
                        }

                        inst_LE.ConsSurplus();
                        std::string res_file = datadir + tblname;
                        inst_LE.PrintRes(res_file);
                    }
                }
            }
        }

        C.disconnect();
        std::cout << "Database operation successfull" << std::endl;
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}
