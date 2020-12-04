#ifndef INCLUDED_MESH
#define INCLUDED_MESH

#include <vector>

#include "config.hpp"

/*****************************************
 * Classe pour créer le maillage
 */
class Mesh{
    public:
        // Parametres du maillage
        double x_min;            // borne gauche du domaine
        double x_max;            // borne droite du domaine
        double y_min;            // borne du haut du domaine
        double y_max;            // borne du bas du domaine
        int N;                  // nombre de mailles/volumes en verticales
        int M;                  // nombre de mailles/volumes en horizontale
        int n_cells;            // nombre de mailles/volumes totales

        double dx;              // Delta x
        double dy;              // Delta y

        std::vector<double> x;         // abscisses des centres des mailles
        std::vector<double> y;         // ordonnees des centres des mailles

        int **coord;       // indices i,j identifiants chaque maille
        int **neighb;        // numero des 4 mailles voisines

        // Constructeur
        Mesh(const Config &cfg);

        // Creation des differents volumes uniformes
        void create_cells();

        // Affichage du malliage
        void display();
        
        // Destructeur
        virtual ~ Mesh();
};

/***************************************
 * @brief  Fonction qui calcule le numero d'une cellule
 * @param  i: abscisse
 * @param  j: ordonnee
 * @param  n_rows: nombre total d'abscisses
 * @param  n_cols: nombre total d'ordonnees
 * @retval identifiant de la cellule dans le maillage global (contenant les cellules fantomes)
 */
int cell_id(int i, int j, int n_rows, int n_cols);

#endif
