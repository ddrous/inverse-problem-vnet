
#ifndef INCLUDED_SOLVER
#define INCLUDED_SOLVER

#include <vector>
#include <map>

#include "mesh.hpp"
#include "config.hpp"

#include "muParser.h"       // Pour transformer des expressions (chaines de caracteres) en fonctions
#include "cnpy.h"           // Pour directement lire les tenseurs Numpy

/************************************************
 * Classe pour resoudre le probleme
 */
class Solver{
    public:
        // Maillage des volumes finis pour le probleme
        const Mesh * mesh;

        // Parametres physiques
        double c;                   // Vitesse de la lumiere
        double a;                   // Constante de Stefen-Boltzmann

        // Propriete du domaine
        double C_v;                 // Capacite calorifique a volume constant

        // Parametres du probleme
        double CFL;                 // Condition CFL: dt=CFL*dx/c
        double precision;           // Precision sur les calculs
        double t_0;                 // Temps de depart de la simulation
        double t_f;                 // Temps de la simulation

        // Autres proprietes du domaine
        std::string rho_expr;           // Expression de la densite du domaine
        std::string sigma_a_expr;       // Expression de l'opacite d'absorption
        std::string sigma_c_expr;       // Expression de l'opacite de scatering

        // La densite sous forme vectorielle pour l'export
        std::vector<double> rho_vec;

        // Attributs des crenaux places sur la densite
        int n_niche;
        double ** attr;

        // Attributs indiquant si on savegarde les signaux sur tout le domaine en tout temps
        bool save_anim;

        // Solution numerique
        std::vector<double> E;              // Solution numerique du problme: E(t, x)
        std::vector<std::vector<double>> F;              // Solution numerique du problme: F(t, x)
        std::vector<double> T;              // Solution numerique du problme: T(t, x)

        // Au temps initial
        std::string E_0_expr;               // Au temps initial: E(t_0, x, y)
        std::string F_0_x_expr;               // Au temps initial: F(t_0, x, y)
        std::string F_0_y_expr;               // Au temps initial: F(t_0, x, y)
        std::string T_0_expr;               // Au temps initial: T(t_0, x, y)

        // Conditions imposees sur le bord du haut
        std::string E_u_expr;               // E(t, x, y_min)
        std::string F_u_x_expr;               // F(t, x, y_min)
        std::string F_u_y_expr;               // F(t, x, y_min)
        std::string T_u_expr;               // T(t, x, y_min)

        // Conditions imposees sur le bord du bas
        std::string E_d_expr;               // E(t, x, y_max)
        std::string F_d_x_expr;               // F(t, x, y_max)
        std::string F_d_y_expr;               // F(t, x, y_max)
        std::string T_d_expr;               // T(t, x, y_max)

        // Conditions imposees sur le bord gauche
        std::string E_l_expr;               // E(t, x_min, y)
        std::string F_l_x_expr;               // F(t, x_min, y)
        std::string F_l_y_expr;               // F(t, x_min, y)
        std::string T_l_expr;               // T(t, x_min, y)

        // Conditions imposees sur le bord droit
        std::string E_r_expr;               // E(t, x_max, y)
        std::string F_r_x_expr;               // F(t, x_max, y)
        std::string F_r_y_expr;               // F(t, x_max, y)
        std::string T_r_expr;               // T(t, x_max, y)

        // Solution exacte
        std::string E_exact_expr;           // Solution exacte: E_exact(t, x, y)
        std::string F_exact_x_expr;           // Solution exacte: F_exact(t, x, y)
        std::string F_exact_y_expr;           // Solution exacte: F_exact(t, x, y)
        std::string T_exact_expr;           // Solution exacte: T_exact(t, x, y)


        // Resultats a exporter en haut
        double** E_up;         // Sur le bord du haut
        double** F_up;         // Sur le bord du haut
        double** T_up;         // Sur le bord du haut

        // Resultats a exporter en bas
        double** E_down;        // Sur le bord du bas
        double** F_down;        // Sur le bord du bas
        double** T_down;        // Sur le bord du bas

        // Resultats a exporter a gauche
        double** E_left;         // Sur le bord gauche
        double** F_left;         // Sur le bord gauche
        double** T_left;         // Sur le bord gauche

        // Resultats a exporter a gauche
        double** E_right;        // Sur le bord droit
        double** F_right;        // Sur le bord droit
        double** T_right;        // Sur le bord droit

        // Autres parametres du probleme
        double dt;                      // Pas de temps
        int step_count;                 // Nombre d'iterations en temps
        std::vector<double> time_steps; // Les temps a chaque pas

        /***************
         * Constructeur
         * @param params: Tous les parametres du probleme
         */
        Solver(const Mesh *new_mesh, const Config &cfg);

        /**
         * Importe rho depuis un tenseur numpy
         * param @filename nom du fichier numpy
         * retourne un vecteur contenant le signal rho
         */ 
        std::vector<double> rho_numpy(std::string filename);

        /**
         * Calcule rho sous forme de fonction crenaux
         * param @n_niche nombre de crenaux
         * param @n_smooth nombre de lissage a effectuer sur le signal
         * retourne un vecteur contenant le signal
         */ 
        std::vector<double> niche(int nb_niche, int nb_smooth);

        /***************
         * Fonction pour calculer rho a partir de son expression rho_expr
         */
        double rho(int i, int j);

        /***************
         * Fonction pour calculer sigma_a
         */
        double sigma_a(double rho, double T);

        /***************
         * Fonction pour calculer sigma_c
         */        
        double sigma_c(double rho, double T);

        /***************
         * Calcule E(t_0, x), energie a la position x au temps initial
         */ 
        double E_0(double x, double y);

        /***************
         * Calcule F(t_0, x)
         */ 
        std::vector<double> F_0(double x, double y);

        /***************
         * Calcule T(t_0, x)        // Eviter T_0 = 0 !
         */ 
        double T_0(double x, double y);

        /**
         * Calcule la source E sous forme de fonction ponctuelle
         * param @t temps de calcul de ce signal
         * param @i la position de calcul du signal
         * param @start debut de l'onde sur le haut
         * param @end fin de l'onde sur le haut
         * param @edge cote sur lequel la source s'applique
         * param @edge_length longueur du cote
         * retourne un vecteur contenant le signal
         */ 
        double ponctual_source(int edge_id, double start, double end, double t, int i);

        /***************
         * Calcule E(t, x, y_max), energie a la position y_max en tout temps
         */ 
        double E_u(double t, int i);

        /***************
         * Calcule F(t, x, y_max), flux a la position y_max en tout temps
         */ 
        std::vector<double> F_u(double t, int i);

        /***************
         * Calcule T(t, x, y_max), temperature a la position y_max en tout temps
         */ 
        double T_u(double t, int i);

        /***************
         * Calcule E(t, x, y_min), energie a la position y_min en tout temps
         */ 
        double E_d(double t, int i);

        /***************
         * Calcule F(t, x, y_min), flux a la position y_min en tout temps
         */ 
        std::vector<double> F_d(double t, int i);

        /***************
         * Calcule T(t, x, y_min), temperature a la position x_min en tout temps
         */ 
        double T_d(double t, int i);

        /***************
         * Calcule E(t, x_min, y), energie a la position x_min en tout temps
         */ 
        double E_l(double t, int j);

        /***************
         * Calcule F(t, x_min, y), flux a la position x_min en tout temps
         */ 
        std::vector<double> F_l(double t, int j);

        /***************
         * Calcule T(t, x_min, y), temperature a la position x_min en tout temps
         */ 
        double T_l(double t, int j);

        /***************
         * Calcule E(t, x_max, y), energuie a la position x_max en tout temps
         */ 
        double E_r(double t, int j);

        /***************
         * Calcule F(t, x_max, y), flux a la position x_max en tout temps
         */ 
        std::vector<double> F_r(double t, int j);

        /***************
         * Calcule T(t, x_max, y), temperature a la position x_max en tout temps
         */ 
        double T_r(double t, int j);

        /***************
         * Calcule E(t, x, y), solution exacte
         */ 
        double E_exact(double t, double x, double y);

        /***************
         * Calcule F(t, x, y), solution exacte
         */ 
        std::vector<double> F_exact(double t, double x, double y);

        /***************
         * Calcule T(t, x, y), solution exacte
         */ 
        double T_exact(double t, double x, double y);

        /**
         * Exporte les solution spatiales pour chaque pas de temps
         */
        void save_animation(int time_step);

        /***************
         * Etape 1 de la resolution su probleme
         */
        void phase_1();

        /***************
         * Etape 1 qui preserve l'equilibre radiatif
         */
        void phase_1_eq();

        /***************
         * Etape 2 de la resolution su probleme
         */
        void phase_2();

        /***************
         * Resout le probleme et exporte les resultats au fur et a mesure
         */
        void solve();

        /***************
         * Affiche les resultats sur la console
         */
        void display();

        /***************
         * Destructeur vide
         */
        virtual ~ Solver();
};

/**
 * extrait le debut de la source et sa fin a l'interieur de l'expression d'un signal ponctuel (ex: pontuel(0.1,0.3))
 */ 
std::vector<double> parse_ponctual(std::string expr);

#endif
