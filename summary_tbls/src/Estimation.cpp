#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
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

#include "Estimation.hpp"

namespace ublas = boost::numeric::ublas;


Estimation::Estimation(const std::string start_date_, const std::string end_date_, const unsigned periodicity_, const unsigned int num_threads)
{
}
