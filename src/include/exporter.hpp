#ifndef INCLUDED_EXPORT
#define INCLUDED_EXPORT

#include "solver.hpp"

/************************************************
 * Classe pour exporter les donnees
 */
class Exporter{
    public:
        // solveur dont on eporte les resultats
        const Solver *solver;
        /***************
         * Constructeur
         */
        Exporter(const Solver *new_solver);

        /**************
         * Exporte une dataframe de donnees spatiales et temporelles
         */
        void write_dataframe(std::string file_name, std::string mode);

        /**************
         * Exporte n_simu simulations au format binaire (format sds - source, densite, signal)
         */
        void write_binary(std::string file_name, std::string mode, std::string n_simu);

        /***************
         * Destructeur
         */
        virtual ~Exporter(){};
};

#endif
