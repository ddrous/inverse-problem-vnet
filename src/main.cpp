#include <iostream>
#include <fstream>

#include "include/config.hpp"
#include "include/mesh.hpp"
#include "include/solver.hpp"
#include "include/exporter.hpp"

using namespace std;
using namespace mu;


int main(int argc, char * argv[]){
    try {
        cout << "\n============= Equation du transfer radiatif en 2D =============" << endl;

        if(argc < 2 || argc > 2)
            throw string("ERREUR: Fournissez un (et un seul) fichier de configuration");


        /* Lecture du fichier config */
        Config cfg = Config(argv[1]);
        cout << "\nLecture des parametres en cours ... " << endl;
        cfg.read();
        cfg.check_numpy();
        cout << "Lecture de " << cfg.n_param << " paramètres OK !" << endl;


        /* Creation d'un maillage uniforme */
        Mesh m = Mesh(cfg);
        cout << "\nCreation des mailles en cours ... " << endl;
        cout << "  -- Suivant l'horizontale: N = " << m.N << endl;
        cout << "  -- Suivant la verticale:  M = " << m.M << endl;
        cout << "  -- Nombre de mailles total (y compris les mailles fantomes): " << m.n_cells << endl;
        m.create_cells();
        // m.display();
        cout << "Creation OK !" << endl;


        /* Resolution du probleme */
        Solver s = Solver(&m, cfg);
        cout << "\nResolution (" << s.step_count << " iterations) en cours ..." << endl;
        s.solve();
        // s.display();
        cout << "Resolution OK !" << endl;


        /* Exportation */
        Exporter ex = Exporter(&s);
        cout << "\nExportation des resultats en cours ..." << endl;
        if (cfg.values["export_mode"] == "dataframe")
            ex.write_dataframe(cfg.values["export_file"], cfg.values["write_mode"]);
        else if (cfg.values["export_mode"] == "binary")
            ex.write_binary(cfg.values["export_file"], cfg.values["write_mode"], cfg.values["simu_count"]);
        else
            throw string ("ERREUR: Mode d'export '"+cfg.values["export_mode"]+"' inconnu. Valeurs possibles: 'datafrane' ou 'binary'");
        
        cout << "  -- Signaux exportés dans le fichier '" << cfg.values["export_file"] << "'"  << endl;
        cout << "Export OK !" << endl;

        cout << "\n======================================================================"  << endl << endl;
    }
    catch(const string &e){
        cerr << endl << e << endl << endl;
        exit(1);
    }
    catch (Parser::exception_type &e){
        cerr << endl << "ERREUR: Dans l'analyseur de fonctions: " << e.GetMsg() << endl << endl;
        exit(1);
    }
    catch(const exception &e){
        cerr << endl << "ERREUR: \n" << e.what() << endl << endl;
        exit(1);
    }

    return 0;
}
