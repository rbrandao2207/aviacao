#ifndef DBFUNCTIONSHEADERDEF
#define DBFUNCTIONSHEADERDEF

#include <pqxx/pqxx>
#include <string>

void csvimport_anac(pqxx::work& W, std::string datadir, std::string year, std::string month, bool createtbl=false);

void csvimport_abear(pqxx::work& W, std::string datadir);

#endif
