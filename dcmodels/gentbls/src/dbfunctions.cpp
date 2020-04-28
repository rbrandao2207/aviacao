#include <iostream>
#include <pqxx/pqxx>
#include <string>

#include "dbfunctions.hpp"

using std::string;


void csvimport_anac(pqxx::work& W, string datadir, string year, string month, bool createtbl)
{
    if (createtbl) {
        string sql_createtbl = "DROP TABLE IF EXISTS csv_anac;"
            "CREATE TABLE csv_anac ("
            "ano         integer,"
            "mes         integer,"
            "empresa     text,"
            "origem      text,"
            "destino     text,"
            "tarifa      real,"
            "assentos    integer"
            ");";

        W.exec(sql_createtbl);
    }

    string sql_copy = "COPY csv_anac FROM '" + datadir + year + month + ".csv' delimiter ';' csv header;";
    W.exec(sql_copy);
    std::cout << "Arquivo csv carregado com sucesso: anac " << year << month << std::endl;
}

void csvimport_abear(pqxx::work& W, string datadir)
{
    string sql_createtbl = "DROP TABLE IF EXISTS csv_abear;"
        "CREATE TABLE csv_abear("
        "ano           integer,"
        "mes           integer,"
        "origem        text,"
        "destino       text,"
        "retorno       text,"
        "tarifa_rt     real,"
        "passageiros   integer,"
        "antecedencia  integer,"
        "estadia       integer"
        ");";

    W.exec(sql_createtbl);

    string sql_copy = "COPY csv_abear FROM '" + datadir + "abear.csv' delimiter ';' csv header;";
    W.exec(sql_copy);
    std::cout << "Arquivo csv carregado com sucesso: abear.csv" << std::endl;
}
