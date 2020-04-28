#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "auxFunctions.hpp"


std::vector<int> getArraySizes(pqxx::nontransaction& N, std::string tblname)
{
    std::vector<int> ret;
    // get number of rows
    std::string sql_count = "SELECT COUNT(*) FROM " + tblname + ";";
    pqxx::result R(N.exec(sql_count));
    auto c = R.begin();
    ret.push_back(c[0].as<int>());

    // get number of coefficients
    std::string sql_count2 = "SELECT * FROM " + tblname + " LIMIT 1;";
    pqxx::result R2(N.exec(sql_count2));
    ret.push_back(R2.begin().size() - 6);

    // get number of products
    std::string sql_count3 = "SELECT id FROM " + tblname + " WHERE obs = 0;";
    pqxx::result R3(N.exec(sql_count3));
    ret.push_back(R3.size());

    return ret;
}

std::string getSyscallOutput(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

double Exp(double x)
{
    return std::exp(x);
}
