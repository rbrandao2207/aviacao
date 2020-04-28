#ifndef ATPIHEADERDEF
#define ATPIHEADERDEF

#include <iostream>
#include <pqxx/pqxx>
#include <string>
#include <vector>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>


class Atpi
{
public:
    Atpi(const std::string& prods_file, const std::string& tbls_file, const std::string& nperiods_file, const unsigned int& inst_nbr, const unsigned int& num_threads, const std::string& data_dir_);
    ~Atpi()
    {
        uvals.clear();
        qties.clear();
        c12.clear();
        weights.clear();
        uidxs.clear();
        std::cout << "Instance ending on table " << tables.back() << " gone!" << std::endl;
    }
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & weights;
        ar & uidxs;
    }
    int FillMats(unsigned qty_threshold);
    int CalcIdxs();
    std::string persist_file;
private:
    std::vector<std::string> products;
    std::vector<std::string> tables;
    std::string data_dir;
    boost::numeric::ublas::matrix<double> uvals;
    boost::numeric::ublas::matrix<unsigned> qties;
    boost::numeric::ublas::matrix<short> c12;
    boost::numeric::ublas::matrix<long double> weights;
    boost::numeric::ublas::matrix<double> uidxs;
};

#endif
