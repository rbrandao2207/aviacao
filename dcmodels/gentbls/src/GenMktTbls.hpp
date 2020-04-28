#ifndef GENMKTTBLSHEADERDEF
#define GENMKTTBLSHEADERDEF

#include <pqxx/pqxx>
#include <string>
#include <vector>


class GenMktTbls
{
public:
    GenMktTbls(pqxx::nontransaction& N, pqxx::work& W, const std::string& cidade_orig, const std::string& cidade_dest, const std::string& period, bool createtbl=true);
    void datafill(pqxx::nontransaction& N, pqxx::work& W, const std::string& period);
    void priceavgs(pqxx::nontransaction& N, pqxx::work& W);

private:
    std::vector<std::string> origens;
    std::vector<std::string> destinos;
    std::vector<std::string> produtos;
    std::string tblname;
}; 
#endif
