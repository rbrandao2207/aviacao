#ifndef AUXFUNCTIONSHEADERDEF
#define AUXFUNCTIONSHEADERDEF

#include <string>

int getProds(const std::string& start_date, const std::string& end_date, const std::string& tbls_file, const std::string& prods_file, const std::string& nperiods_file);

int asyncFillMats(std::string prods_file, std::string tbls_file, std::string nperiods_file, unsigned int inst_nbr, unsigned int num_threads, std::string data_dir, unsigned int qty_threshold);

int asyncCalcIdxs(std::string prods_file, std::string tbls_file, std::string nperiods_file, unsigned int inst_nbr, unsigned int num_threads, std::string data_dir);

#endif
