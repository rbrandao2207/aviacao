#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <pqxx/pqxx>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "Atpi.hpp"
#include "AuxFunctions.hpp"


// current run options:
// 1) argv[1] getprods & argv[2] fillmats
// 2) argv[1] fillmats (needs getprods files)
// 3) argv[1] idxs (needs previous serialization)

int main(int argc, char* argv[])
{
    // params and initialization
    // available dates are hardcoded in function getProds (current first/last = Jan2002/May2019)
    std::string start_date = "200201";
    std::string end_date = "201905";
    unsigned int qty_threshold = 0;
    unsigned int max_threads = 64;
    unsigned int hardware_threads = std::thread::hardware_concurrency();
    unsigned int num_real_threads = std::min(hardware_threads != 0 ? hardware_threads : 2, max_threads);
    unsigned int num_threads = 26; //** change between num_real_threads/fixed number, as desired
    auto chrono_start = std::chrono::steady_clock::now();
    std::string data_dir = "/app/results/" + start_date + "-" + end_date + "_threads-" + std::to_string(num_threads);
    std::filesystem::path aux_path(data_dir);
    if (std::filesystem::exists(aux_path) == 0) {
        std::string syscall_mkdir = "mkdir " + data_dir;
        system(syscall_mkdir.c_str());
    }
    std::string tbls_file = data_dir + "/tables";
    std::string prods_file = data_dir + "/products";
    std::string nperiods_file = data_dir + "/nperiods";

    // get products and tables persisted
    if (argc > 1 && std::strcmp(argv[1], "getprods") == 0) {
        std::remove(tbls_file.c_str());
        std::remove(prods_file.c_str());
        std::remove(nperiods_file.c_str());
        getProds(start_date, end_date, tbls_file, prods_file, nperiods_file);
    } else {
        if (!((access(tbls_file.c_str(), F_OK) != -1) && (access(prods_file.c_str(), F_OK) != -1)) && (access(prods_file.c_str(), F_OK) != -1)) // check if files are really there
            throw std::runtime_error("some file not found; run w/ getprods arg instead");

    }

    //** asynchronous computation of member function chosen
    std::vector<std::future<int>> futures;

    // compute and serialize matrices uvals and qties
    if ((argc > 2 && std::strcmp(argv[2], "fillmats") == 0) || (argc > 1 && std::strcmp(argv[1], "fillmats") == 0)) {
        for (unsigned int inst_nbr = 0; inst_nbr < num_threads; ++inst_nbr) {
            futures.push_back(std::async(asyncFillMats, prods_file, tbls_file, nperiods_file, inst_nbr, num_threads, data_dir, qty_threshold));
        }

        for (auto& entry : futures)
            entry.wait();

    // compute and serialize matrices weights
    } else if (argc > 1 && std::strcmp(argv[1], "idxs") == 0) {
        for (unsigned int inst_nbr = 0; inst_nbr < num_threads; ++inst_nbr) {
            futures.push_back(std::async(asyncCalcIdxs, prods_file, tbls_file, nperiods_file, inst_nbr, num_threads, data_dir));
        }

        for (auto& entry : futures)
            entry.wait();

    // else args are not accepted
    } else {
        std::cout << "Invalid main args..." << std::endl;
    }

    // finish chrono
    auto chrono_end = std::chrono::steady_clock::now();
    auto time_diff = chrono_end - chrono_start;
    std::string time_fpersist = data_dir + "/ellapsed_time";
    if (argc > 1 && std::strcmp(argv[1], "idxs") == 0)
        time_fpersist += "_idxs";
    std::remove(time_fpersist.c_str());
    std::ofstream fdesc_time;
    fdesc_time.open(time_fpersist);
    fdesc_time << std::chrono::duration_cast<std::chrono::minutes> (time_diff).count();
    fdesc_time.close();
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::minutes> (time_diff).count() << " mins" << std::endl;

    return 0;
}
