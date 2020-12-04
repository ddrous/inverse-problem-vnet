#include <string>
#include <iostream>
#include <iomanip> 

#include "include/mesh.hpp"

using namespace std;

Mesh::Mesh(const Config &cfg){
    x_min = atof(cfg.values.at("x_min").c_str());
    x_max = atof(cfg.values.at("x_max").c_str());
    y_min = atof(cfg.values.at("y_min").c_str());
    y_max = atof(cfg.values.at("y_max").c_str());
    N = atoi(cfg.values.at("N").c_str());
    M = atoi(cfg.values.at("M").c_str());

    if (x_min >= x_max)
        throw std::string("ERREUR: Vérifiez que x_min < x_max");
    if (y_min >= y_max)
        throw std::string("ERREUR: Vérifiez que y_min < y_max");
    if (N <= 0)
        throw std::string("ERREUR: Vérifiez que N > 0");
    if (M <= 0)
        throw std::string("ERREUR: Vérifiez que M > 0");

    x = std::vector<double>(N+2);
    y = std::vector<double>(M+2);

    dx = (x_max - x_min) / N;
    dy = (y_max - y_min) / M;

    n_cells = (N+2) * (M+2);
    coord = new int*[n_cells];
    neighb = new int*[n_cells];
    for (int k = 0; k < n_cells; k++){
        coord[k] = new int[2];
        neighb[k] = new int[4];
    }
}


int cell_id(int i, int j, int n_rows, int n_cols){
    return i + j*n_rows;
}


void Mesh::create_cells(){
    for (int i = 0; i < N+2; i++)
        x[i] = x_min + (i-1)*dx + dx/2.;     // centres des mailles suivant l'horizontale i

    for (int j = 0; j < M+2; j++)
        y[j] = y_min + (j-1)*dy + dy/2.;     // centres des mailles suivant la verticale j

    /* Remplissage des coordonnes et des voisins (-1 indique un voisin absent) */
    for (int i = 0; i < N+2; i++)
        for (int j = 0; j < M+2; j++){

            int k = cell_id(i, j, N+2, M+2);        // numero (identifiant) de la maille
            coord[k][0] = i;
            coord[k][1] = j;

            if (j == M+1){            // Les mailles fantomes du haut
                neighb[k][0] = -1;               // voisin du haut
                neighb[k][1] = k - (N+2);        // voisin du bas
                neighb[k][2] = -1;               // voisin de gauche
                neighb[k][3] = -1;               // voisin de droite
            }

            else if (j == 0){            // Les mailles fantomes du bas
                neighb[k][0] = k + (N+2);
                neighb[k][1] = -1;
                neighb[k][2] = -1;
                neighb[k][3] = -1;
            }

            else if (i == 0){            // Les mailles fantomes de gauche
                neighb[k][0] = -1;
                neighb[k][1] = -1;
                neighb[k][2] = -1;
                neighb[k][3] = k + 1;
            }

            else if (i == N+1){            // Les mailles fantomes de droite
                neighb[k][0] = -1;
                neighb[k][1] = -1;
                neighb[k][2] = k - 1;
                neighb[k][3] = -1;
            }

            else {      // Les mailles intermediaires
                neighb[k][0] = k + (N+2);        // voisin du haut
                neighb[k][1] = k - (N+2);      // voisin du bas
                neighb[k][2] = k - 1;            // voisin de gauche
                neighb[k][3] = k + 1;            // voisin de droite
            }
        }
}


void Mesh::display(){
    cout << "-----------  Maillage  -----------\n" ;     // up-down-left-right, -1 pour les voisins manquants
    for (int j = M+1; j >= 0; j--){
        for (int i = 0; i <= N+1; i++){
            int k = cell_id(i, j, N+2, M+2);
            cout << setw(8) << "(" << setw(1) << coord[k][0] << "," << setw(1) << coord[k][1] << ") : ";
            cout << setw(2) << k << " : (" << setw(2) << neighb[k][0] << "," << setw(2) << neighb[k][1] << "," << setw(2) << neighb[k][2] << "," << setw(2) << neighb[k][3] << ")";
        }
        if (j != 0) cout << "\n\n";
        else cout << "\n";
    }
}


Mesh::~Mesh(){
    for (int k = 0; k < n_cells; k++){
        delete[] coord[k];
        delete[] neighb[k];
    }
    delete[] coord;
    delete[] neighb;
}
