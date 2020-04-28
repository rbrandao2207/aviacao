#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <pqxx/pqxx>
#include <string>
#include <vector>

#include "Atpi.hpp"
#include "AuxFunctions.hpp"


// persists two files containing products and tables data
int getProds(const std::string& start_date, const std::string& end_date, const std::string& tbls_file, const std::string& prods_file, const std::string& nperiods_file)
{
    try {
        pqxx::connection C("dbname = aviacao user = postgres password = passwd hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            std::cout << "Opened database successfully: " << C.dbname() << std::endl;
        } else {
            std::cout << "Can't open database" << std::endl;
            return 1;
        }
        // hardcoded periods bellow, to be updated w/ new ANAC releases
        std::vector<std::string> yrs = {"2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019"};
        std::vector<std::string> mths = {"01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"};

        pqxx::nontransaction N(C);
        int bool_start_date= -1;
        unsigned int count_nperiods = 0;
        for (std::string yr : yrs) {
            for (std::string mth : mths) {
                if (yr + mth == start_date)
                    bool_start_date = 0;
                if (bool_start_date == -1)
                    continue;
                if (yr + mth == end_date)
                    bool_start_date = -1;
                ++count_nperiods;

                // persist tables
                std::ofstream fdesc_tbls;
                fdesc_tbls.open(tbls_file, std::ios_base::app);
                assert(fdesc_tbls.is_open());
                fdesc_tbls << "csv_" + yr + mth << '\n';
                fdesc_tbls.close();

                // get products and persist
                std::string query = "SELECT DISTINCT origem, destino, empresa FROM csv_" + yr + mth + ";";
                pqxx::result R(N.exec(query));
                std::ofstream fdesc_prods;
                fdesc_prods.open(prods_file, std::ios_base::app);
                assert(fdesc_prods.is_open());
                for (auto c = R.begin(); c != R.end(); ++c) {
                    fdesc_prods << c[0].as<std::string>() << ',' << c[1].as<std::string>() << ',' << c[2].as<std::string>() << '\n';
                }
                fdesc_prods.close();
            }
        }

        // persist nperiods
        std::ofstream fdesc_npers;
        fdesc_npers.open(nperiods_file);
        assert(fdesc_npers.is_open());
        fdesc_npers << count_nperiods << '\n';
        fdesc_npers.close();

        // sort and eliminate duplicates: prods_file
        std::string syscall_sort = "sort -u " + prods_file + " > " + prods_file + "2; mv " + prods_file + "2 " + prods_file;
        system(syscall_sort.c_str());
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}

int asyncFillMats(std::string prods_file, std::string tbls_file, std::string nperiods_file, unsigned int inst_nbr, unsigned int num_threads, std::string data_dir, unsigned int qty_threshold)
{
    Atpi instance(prods_file, tbls_file, nperiods_file, inst_nbr, num_threads, data_dir);
    instance.FillMats(qty_threshold);
    // serialize
    {
        std::string aux_file;
        aux_file = instance.persist_file;
        std::remove(aux_file.c_str());
        std::ofstream ofs(aux_file);
        assert(ofs.is_open());
        boost::archive::text_oarchive oa(ofs);
        oa << instance;
    }

    return 0;
}

int asyncCalcIdxs(std::string prods_file, std::string tbls_file, std::string nperiods_file, unsigned int inst_nbr, unsigned int num_threads, std::string data_dir)
{
    Atpi instance(prods_file, tbls_file, nperiods_file, inst_nbr, num_threads, data_dir);
    // deserialize
    {
        std::string aux_file = instance.persist_file;
        std::ifstream ifs(aux_file);
        assert(ifs.is_open());
        boost::archive::text_iarchive ia(ifs);
        ia >> instance;
    }
    instance.CalcIdxs();

    return 0;
}
