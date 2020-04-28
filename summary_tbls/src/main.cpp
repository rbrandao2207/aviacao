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

#include "Estimation.hpp"
#include "AuxFunctions.hpp"


// current run options:
// 1) argv[1] gentbls

int main(int argc, char* argv[])
{
    auto chrono_start = std::chrono::steady_clock::now();

    // hardcoded periods bellow, to be updated w/ new ANAC releases
    // current first = Jan2002; current last = May2019
    std::vector<std::string> yrs = {"2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019"};
    std::vector<std::string> mths = {"01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"};

    //** global params (additional params for gentbls below)
    std::string start_date = "200201";
    std::string end_date = "201911";
    unsigned periodicity = 3; // # of months to aggregate
    unsigned max_threads = 64;
    unsigned hardware_threads = std::thread::hardware_concurrency();
    unsigned num_real_threads = std::min(hardware_threads != 0 ? hardware_threads : 2, max_threads);
    unsigned num_threads = 1; //** change between num_real_threads/fixed number, as desired
    std::string results_dir = "/app/results/" + start_date + "-" + end_date + "-" + std::to_string(periodicity);

    // initialization
    std::vector<std::string> tables;
    bool bool_tbl = false;
    for (auto& yr : yrs) {
        for (auto& mth : mths) {
            if (yr + mth == start_date)
                bool_tbl = true;
            if (!bool_tbl)
                continue;
            if (yr + mth == end_date)
                bool_tbl = false;
            tables.push_back(yr + mth);
        }
    }

    // argv[1] gentbls
    if (argc > 1 && std::strcmp(argv[1], "gentbls") == 0) {
        //** params
        unsigned qty_threshold = 0;
        std::vector<double> bins = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};

        // initialization
        std::vector<std::future<int>> futures;
        unsigned quotient = tables.size() / num_threads;
        unsigned remainder = tables.size() % num_threads;
        unsigned count_quo = 0;
        unsigned count_thread= 0;
        unsigned count_arg = 0;

        if (quotient == 0) {
            std::vector<std::string> thread_tbls;
            thread_tbls = tables;
            futures.push_back(std::async(genTbls, qty_threshold, bins, thread_tbls));
        } else {
            while (count_arg < tables.size()) {
                std::vector<std::string> thread_tbls;
                for (count_quo = 0; count_quo < quotient; ++count_quo) {
                    thread_tbls.push_back(tables[count_arg]);
                    ++count_arg;
                }
                // (remainder is distributed evenly to foremost threads)
                if (count_thread < remainder) {
                    thread_tbls.push_back(tables[count_arg]);
                    ++count_arg;
                }
                // do async job
                futures.push_back(std::async(genTbls, qty_threshold, bins, thread_tbls));
                ++count_thread;
            }
        }

        for (auto& entry : futures)
            entry.wait();

    }

    // finish chrono
    auto chrono_end = std::chrono::steady_clock::now();
    auto time_diff = chrono_end - chrono_start;
    std::string time_fpersist = results_dir + "/ellapsed_time";
    std::remove(time_fpersist.c_str());
    std::ofstream fdesc_time;
    fdesc_time.open(time_fpersist);
    fdesc_time << std::chrono::duration_cast<std::chrono::minutes> (time_diff).count();
    fdesc_time.close();
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::minutes> (time_diff).count() << " mins" << std::endl;

    return 0;
}
