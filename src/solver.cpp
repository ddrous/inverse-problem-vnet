#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <omp.h>

#include "include/solver.hpp"
#include "include/vector_ops.hpp"

using namespace mu;
using namespace std;

/* Types utilisés pour le stockage des signaux 1d ou 2d */
typedef std::vector<double> vector_1row;                // de taille n_cells ou time_steps
typedef std::vector<double> vector_2d;                  // de taille 2 (x et y)
typedef std::vector<vector_2d> vector_2rows;            // de taille n_cells * 2


/**
 * Allocate space for 2D matrix
 */ 
double ** allocate(int nrows, int ncols){
    double ** matrix = new double*[nrows];
    for (int i = 0; i < nrows; i++)
        matrix[i] = new double[ncols];
    return matrix;
}


/**
 * Free space from 2D matrix
 */ 
void free(double **matrix, int nrows){
    for (int i = 0; i < nrows; i++)
        delete[] matrix[i];
    delete[] matrix;
}


Solver::Solver(const Mesh *new_mesh, const Config &cfg){
    mesh = new_mesh;

    c = atof(cfg.values.at("c").c_str());
    a = atof(cfg.values.at("a").c_str());

    C_v = atof(cfg.values.at("C_v").c_str());

    CFL = atof(cfg.values.at("CFL").c_str());
    precision = atof(cfg.values.at("precision").c_str());
    t_0 = atof(cfg.values.at("t_0").c_str());
    t_f = atof(cfg.values.at("t_f").c_str());

    rho_expr = cfg.values.at("rho");
    sigma_a_expr = cfg.values.at("sigma_a");
    sigma_c_expr = cfg.values.at("sigma_c");

    rho_vec = vector_1row(mesh->n_cells);

    E = vector_1row(mesh->n_cells);
    F = vector_2rows(mesh->n_cells, vector_2d(2));
    T = vector_1row(mesh->n_cells);

    E_0_expr = cfg.values.at("E_0");
    F_0_x_expr = cfg.values.at("F_0_x");
    F_0_y_expr = cfg.values.at("F_0_y");
    T_0_expr = cfg.values.at("T_0");

    E_u_expr = cfg.values.at("E_u");
    F_u_x_expr = cfg.values.at("F_u_x");
    F_u_y_expr = cfg.values.at("F_u_y");
    T_u_expr = cfg.values.at("T_u");

    E_d_expr = cfg.values.at("E_d");
    F_d_x_expr = cfg.values.at("F_d_x");
    F_d_y_expr = cfg.values.at("F_d_y");
    T_d_expr = cfg.values.at("T_d");

    E_l_expr = cfg.values.at("E_l");
    F_l_x_expr = cfg.values.at("F_l_x");
    F_l_y_expr = cfg.values.at("F_l_y");
    T_l_expr = cfg.values.at("T_l");

    E_r_expr = cfg.values.at("E_r");
    F_r_x_expr = cfg.values.at("F_r_x");
    F_r_y_expr = cfg.values.at("F_r_y");
    T_r_expr = cfg.values.at("T_r");

    E_exact_expr = cfg.values.at("E_exact");
    F_exact_x_expr = cfg.values.at("F_exact_x");
    F_exact_y_expr = cfg.values.at("F_exact_y");
    T_exact_expr = cfg.values.at("T_exact");

    /* Verifications preliminaires */
    if (c <= 0)
        throw string("ERREUR: Vérifiez que c > 0");
    if (a <= 0)
        throw string("ERREUR: Vérifiez que a > 0");
    if (C_v <= 0)
        throw string("ERREUR: Vérifiez que C_v > 0");
    if (CFL <= 0)
        throw string("ERREUR: Vérifiez que CFL > 0");
    if (precision <= 0)
        throw string("ERREUR: Vérifiez que la precision est > 0");
    if (t_0 < 0)
        throw string("ERREUR: Vérifiez que t_0 >= 0");
    if (t_f <= 0)
        throw string("ERREUR: Vérifiez que t_f > 0");

    /* Les variables de temps */
    dt = CFL * mesh->dx/c;
    // dt = CFL * mesh->dx * mesh->dy /c;                  // condition CFL comme ceci plutot?

    double tmp = t_f / dt;
    if (tmp == int(tmp))
        step_count = int(tmp);
    else
        step_count = int(tmp) + 1;

    time_steps = vector_1row(step_count);

    /* Nombre de créneaux par défauts */
    n_niche = 0;        

    /* Des animations seulement quand il s'agit du mode dataframe*/
    save_anim = (cfg.values.at("export_mode") == "dataframe");

    /* A exporter */
    E_up = allocate(step_count, mesh->N);
    F_up = allocate(step_count, mesh->N);
    T_up = allocate(step_count, mesh->N);

    E_down = allocate(step_count, mesh->N);
    F_down = allocate(step_count, mesh->N);
    T_down = allocate(step_count, mesh->N);

    E_left = allocate(step_count, mesh->M);
    F_left = allocate(step_count, mesh->M);
    T_left = allocate(step_count, mesh->M);

    E_right = allocate(step_count, mesh->M);
    F_right = allocate(step_count, mesh->M);
    T_right = allocate(step_count, mesh->M);
} 


vector_1row Solver::rho_numpy(std::string filename){
    cnpy::NpyArray npydata = cnpy::npy_load(filename);

    double* data = npydata.data<double>();

    // cout << "shape: " << npydata.shape[0] << ", " << npydata.shape[1] << endl; 

    int N = npydata.shape[0];   	    // = mesh->M d'apres le fichier config 
    int M = npydata.shape[1];       	// = mesh->N d'apres le fichier config
    vector_1row my_vec = vector_1row(mesh->n_cells);

    for (int i = 0; i < N; i++){
        for (int j = 0; j < M; j++){
            int k = j + i*M;                            // formule de lecture numpy (matrices)
            int k_prime = cell_id(j+1, N-i, M+2, N+2);      // conventionn repere orthonormé pour mon maillage
            my_vec[k_prime] = data[k];

            // cout << "data ["<< i << ", " << j << "]: " << "  k=" << k << "   "<< data[k] << "\n";
            // cout << "my vec ["<< j+1 << ", " << 90-i << "]: " << "  k_prime=" << k_prime << "   "<< my_vec[k] << "\n";  
        }
    }
    return my_vec;
}


/**
 * Laplacian smoothing (k-means) d'un vecteur, avec k=3
 */ 
vector_1row smooth(vector_1row& input, int **neighb){
    int size = input.size();
    vector_1row output(size);

    for (int k = 0; k < size; k++){
        double sum = input[k];
        int count_nbr = 1;              // nbr de voisins importants
        for (int l = 0; l < 4; l++){
            int id_nbr = neighb[k][l];
            if (id_nbr != -1){
                sum += input[id_nbr];
                count_nbr += 1;
            }
        }
        output[k] = sum / count_nbr;
    }

    return output;
}


vector_1row Solver::niche(int nb_niche, int nb_smooth){
    /* Vecteur qui va contenir le signal en crenaux */
    vector_1row signal(mesh->n_cells);

    /* Declaration Les attributs du (des) crenaux */
    n_niche = nb_niche;                // nombre de crenaux
    attr = allocate(n_niche, 4);

    /* Parsing de l'expression d'un crenau */
    size_t comma1 = rho_expr.find(",");
    size_t comma2 = rho_expr.find(",", comma1+1);
    size_t comma3 = rho_expr.find(",", comma2+1);

    string pos_x_expr = rho_expr.substr(7, comma1-7);
    string pos_y_expr = rho_expr.substr(comma1+1, comma2-comma1-1);
    string rho_min_expr = rho_expr.substr(comma2+1, comma3-comma2-1);
    string rho_max_expr = rho_expr.substr(comma3+1, rho_expr.length()-comma3-2);

    double pos_x = atof(pos_x_expr.c_str());
    double pos_y = atof(pos_y_expr.c_str());
    double rho_min = atof(rho_min_expr.c_str());
    double rho_max = atof(rho_max_expr.c_str());

    if (pos_x < mesh->x_min || pos_x > mesh->x_max)
        throw string("ERREUR: Abcisse du crénau invalide");
    if (pos_y < mesh->y_min || pos_y > mesh->y_max)
        throw string("ERREUR: Ordonnée du crénau invalide");

    for (int l = 0; l < nb_niche; l++){
        /* Les memes attributs a chaque fois */
        attr[l][0] = pos_x;    // abcisse
        attr[l][1] = pos_y;    // ordonnee
        attr[l][2] = 0.1*((mesh->N + mesh->M)/2) * ((mesh->dx + mesh->dy)/2);   // diametre
        attr[l][3] = rho_max;                                                   // hauteur
    }

    /* Placement des crenaux */
    for (int i = 0; i < mesh->N+2; i++){
        for (int j = 0; j < mesh->M+2; j++){
            int k = cell_id(i, j, mesh->N+2, mesh->M+2);
            signal[k] = rho_min;
            for (int l = 0; l < n_niche; l++){
                if (sqrt(pow(mesh->x[i] - attr[l][0], 2) + pow(mesh->y[j] - attr[l][1], 2)) <= 1*attr[l][2]/2.){                       // crenau circulaire
                // if ((abs(mesh->x[i] - attr[l][0]) <= 2*attr[l][2]/2.) && ( abs(mesh->y[j] - attr[l][1]) <= 5*attr[l][2]/2.) ){            // crenau rectangulaire
                    signal[k] = attr[l][3];
                    break;
                }
            }
        }
    }

    /* Lissage du signal */
    for (int i = 0; i < nb_smooth; i++)
        signal = smooth(signal, mesh->neighb);

    /* Mettre a jour la hauteur de crenau */
    for (int l = 0; l < nb_niche; l++){
        attr[l][2] = *min_element(signal.begin(), signal.end());
        attr[l][3] = *max_element(signal.begin(), signal.end());
    }


    return signal;
}


double Solver::rho(int i, int j){
    static int first_call = true;
    static bool niche_bool = (rho_expr.find("crenau") != string::npos);
    static bool numpy_bool = (rho_expr.find("numpy") != string::npos);

    int k = cell_id(i, j, mesh->N+2, mesh->M+2);

    if (first_call){
        if (niche_bool){
            rho_vec = niche(1, 0.1*(mesh->N + mesh->M)/2.);
        }

        else if (numpy_bool){
            rho_vec = rho_numpy(rho_expr.substr(6));
        }
        
        else{
            static Parser p;
                p.SetExpr(rho_expr);
                for (int i = 0; i < mesh->N+2; i++){
                    double x = mesh->x[i];
                    p.DefineVar("x", &x);
                    for (int j = 0; j < mesh->M+2; j++){
                        double y = mesh->y[j];
                        p.DefineVar("y", &y);
                        int k_prime = cell_id(i, j, mesh->N+2, mesh->M+2);
                        rho_vec[k_prime] = p.Eval();
                    }
                }
        }
        first_call = false;
    }
    return rho_vec[k];
}


/** CETTE FONCTION EST IGNORËE AFIN DE PARALLËLISER: Le probleme etant que les
 * parseur muparser ne sont pas facilement copiable (qualité indispendable pour
 * etre utiliser dans firstprivate). Une solution aurait été de remplacer 
 * completement les parseur par nos propres objects. On pourra alors définir 
 * ce qu'est une copie. 
 * Pou le moment, cette focntion et la suivante ne sont ne sont pas utilisée
 * , (voir) phase_1 et phase_2
*/
double Solver::sigma_a(double rho, double T){
    /* VERSION SINGLE CORE */
    static int first_call = 1;
    static Parser p;
    p.DefineVar("rho", &rho); 
    p.DefineVar("T", &T); 
    if (first_call == 1){ p.SetExpr(sigma_a_expr); first_call = 0; }
    return p.Eval();

    /* VERSION PARALLELE (éviter les variables statiques) */
    // Parser p;
    // p.DefineVar("rho", &rho); 
    // p.DefineVar("T", &T); 
    // p.SetExpr(sigma_a_expr); 
    // return p.Eval();

    /* ANOTHER PARALLELIZABLE VERSION */
	// double ret;
    // try
	// {
    //     var1 = rho;
    //     var2 = T;
    //     ret = sigma_a_func.Eval();
	// }
	// catch (mu::Parser::exception_type &e)
	// {
	// 	std::cerr << e.GetMsg() << std::endl;
	// }
    // return ret;
}


double Solver::sigma_c(double rho, double T){
    /* VERSION SINGLE CORE */
    static int first_call = 1;
    static Parser p;
    p.DefineVar("rho", &rho);
    p.DefineVar("T", &T); 
    if (first_call == 1){ p.SetExpr(sigma_c_expr); first_call = 0; }
    return p.Eval();

    /* VERSION PARALLELE */
    // Parser p;
    // p.DefineVar("rho", &rho); 
    // p.DefineVar("T", &T); 
    // p.SetExpr(sigma_c_expr); 
    // return p.Eval();

    /* ANOTHER PARALLELIZABLE VERSION */
	// double ret;
    // try
	// {
    //     var3 = rho;
    //     var4 = T;
    //     ret = sigma_c_func.Eval();
	// }

	// catch (mu::Parser::exception_type &e)
	// {
	// 	std::cerr << e.GetMsg() << std::endl;
	// }
    // return ret;
}


double Solver::E_0(double x, double y){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("x", &x);
    p.DefineVar("y", &y);
    p.DefineVar("t_0", &t_0);
    if (first_call == 1){ p.SetExpr(E_0_expr); first_call = 0; }

    return p.Eval();
}


vector_2d Solver::F_0(double x, double y){
    static int first_call = 1;

    static Parser p1, p2;
    p1.DefineVar("x", &x); p2.DefineVar("x", &x);
    p1.DefineVar("y", &y); p2.DefineVar("y", &y);
    p1.DefineVar("t_0", &t_0); p2.DefineVar("t_0", &t_0);
    if (first_call == 1){ p1.SetExpr(F_0_x_expr); p2.SetExpr(F_0_y_expr); first_call = 0; }

    return {p1.Eval(), p2.Eval()};
}


double Solver::T_0(double x, double y){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("x", &x); 
    p.DefineVar("y", &y);
    p.DefineVar("t_0", &t_0);
    if (first_call == 1){ p.SetExpr(T_0_expr); first_call = 0; }

    return p.Eval();
}


vector_1row parse_ponctual(string expr){
    vector_1row edges(2);

    size_t comma = expr.find(",");
    string start = expr.substr(9, comma-9);
    string end = expr.substr(comma+1, expr.length()-comma-2);
    edges[0] = atof(start.c_str());
    edges[1] = atof(end.c_str());

    return edges;
}


double Solver::ponctual_source(int edge_id, double start, double end, double t, int i){

    int edge_length;        // longeur de cette extremite
    double current_pos;     // position courante (abcisse ou ordonnee en fonction du cas)
    int k;                  // identifiant de la cellule correspondant a i
    switch (edge_id){
        case 0:                     // en haut
            k = cell_id(i, mesh->M+1, mesh->N+2, mesh->M+2);
            edge_length = mesh->N;
            current_pos = mesh->x[i];
            break;
        case 1:                     // en bas
            k = cell_id(i, 0, mesh->N+2, mesh->M+2);
            edge_length = mesh->N;
            current_pos = mesh->x[i];
            break;
        case 2:                     // a gauche
            k = cell_id(0, i, mesh->N+2, mesh->M+2);
            edge_length = mesh->M;
            current_pos = mesh->y[i];
            break;
        case 3:                     // a droite
            k = cell_id(mesh->N+1, i, mesh->N+2, mesh->M+2);
            edge_length = mesh->M;
            current_pos = mesh->y[i];
            break;
        default:
            break;
    }

    if (start < mesh->x_min || start < mesh->y_min)
        throw string("ERREUR: Position inférieure de la source ponctuelle invalide");
    if (end > mesh->x_max || end > mesh->y_max)
        throw string("ERREUR: Position supérieure de la source ponctuelle invalide");

    if (start <= current_pos && current_pos <= end)
        return a*pow(T[k], 4) + 5*sin(2*500*M_PI*t);
    else
        return a*pow(T[k], 4);
}


double Solver::E_u(double t, int i){
    static int first_call = 1;
    static int neumann = E_u_expr.compare("neumann");
    static bool ponctual_bool = (E_u_expr.find("ponctuel") != string::npos);

    if (neumann == 0){
        int k = cell_id(i, mesh->M, mesh->N+2, mesh->M+2);
        return E[k];
    }
    else if (ponctual_bool == true){
        static vector_1row start_end = parse_ponctual(E_u_expr);
        return ponctual_source(0, start_end[0], start_end[1], t, i);
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        double x = mesh->x[i];
        p.DefineVar("x", &x);
        if (first_call == 1){ p.SetExpr(E_u_expr); first_call = 0; }
        return p.Eval();
    }
}


vector_2d Solver::F_u(double t, int i){
    static int first_call = 1;
    static int neumann = F_u_x_expr.compare("neumann");

    if (neumann == 0){
        int k = cell_id(i, mesh->M, mesh->N+2, mesh->M+2);
        return F[k];
    }
    else{ 
        static Parser p1, p2;
        p1.DefineVar("t", &t); p2.DefineVar("t", &t);
        double x = mesh->x[i];
        p1.DefineVar("x", &x); p2.DefineVar("x", &x);
        if (first_call == 1){ p1.SetExpr(F_u_x_expr); p2.SetExpr(F_u_y_expr); first_call = 0; }
        return {p1.Eval(), p2.Eval()};
    }
}


double Solver::T_u(double t, int i){
    static int first_call = 1;
    static int neumann = T_u_expr.compare("neumann");

    if (neumann == 0){
        int k = cell_id(i, mesh->M, mesh->N+2, mesh->M+2);
        return T[k];
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        double x = mesh->x[i];
        p.DefineVar("x", &x);
        if (first_call == 1){ p.SetExpr(T_u_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::E_d(double t, int i){
    static int first_call = 1;
    static int neumann = E_d_expr.compare("neumann");
    static bool ponctual_bool = (E_d_expr.find("ponctuel") != string::npos);

    if (neumann == 0){
        int k = cell_id(i, 1, mesh->N+2, mesh->M+2);
        return E[k];
    }
    else if (ponctual_bool == true){
        static vector_1row start_end = parse_ponctual(E_d_expr);
        return ponctual_source(1, start_end[0], start_end[1], t, i);
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        double x = mesh->x[i];
        p.DefineVar("x", &x);
        if (first_call == 1){ p.SetExpr(E_d_expr); first_call = 0; }
        return p.Eval();
    }
}


vector_2d Solver::F_d(double t, int i){
    static int first_call = 1;
    static int neumann = F_d_x_expr.compare("neumann");

    if (neumann == 0){
        int k = cell_id(i, 1, mesh->N+2, mesh->M+2);
        return F[k];
    }
    else{ 
        static Parser p1, p2;
        p1.DefineVar("t", &t); p2.DefineVar("t", &t);
        double x = mesh->x[i];
        p1.DefineVar("x", &x); p2.DefineVar("x", &x);
        if (first_call == 1){ p1.SetExpr(F_d_x_expr); p2.SetExpr(F_d_y_expr); first_call = 0; }
        return {p1.Eval(), p2.Eval()};
    }
}


double Solver::T_d(double t, int i){
    static int first_call = 1;
    static int neumann = T_d_expr.compare("neumann");

    if (neumann == 0){
        int k = cell_id(i, 1, mesh->N+2, mesh->M+2);
        return T[k];
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        double x = mesh->x[i];
        p.DefineVar("x", &x);
        if (first_call == 1){ p.SetExpr(T_d_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::E_l(double t, int j){
    static int first_call = 1;
    static int neumann = E_l_expr.compare("neumann");
    static bool ponctual_bool = (E_l_expr.find("ponctuel") != string::npos);

    if (neumann == 0){
        int k = cell_id(1, j, mesh->N+2, mesh->M+2);
        return E[k];
    }
    else if (ponctual_bool == true){
        static vector_1row start_end = parse_ponctual(E_l_expr);
        return ponctual_source(2, start_end[0], start_end[1], t, j);
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        double y = mesh->y[j];
        p.DefineVar("y", &y);
        if (first_call == 1){ p.SetExpr(E_l_expr); first_call = 0; }
        return p.Eval();
    }
}


vector_2d Solver::F_l(double t, int j){
    static int first_call = 1;
    static int neumann = F_l_x_expr.compare("neumann");

    if (neumann == 0){
        int k = cell_id(1, j, mesh->N+2, mesh->M+2);
        return F[k];
    }
    else{ 
        static Parser p1, p2;
        p1.DefineVar("t", &t); p2.DefineVar("t", &t);
        double y = mesh->y[j];
        p1.DefineVar("y", &y); p2.DefineVar("y", &y);
        if (first_call == 1){ p1.SetExpr(F_l_x_expr); p2.SetExpr(F_l_y_expr); first_call = 0; }
        return {p1.Eval(), p2.Eval()};
    }
}


double Solver::T_l(double t, int j){
    static int first_call = 1;
    static int neumann = T_l_expr.compare("neumann");

    if (neumann == 0){
        int k = cell_id(1, j, mesh->N+2, mesh->M+2);
        return T[k];
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        double y = mesh->y[j];
        p.DefineVar("y", &y);
        if (first_call == 1){ p.SetExpr(T_l_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::E_r(double t, int j){
    static int first_call = 1;
    static int neumann = E_r_expr.compare("neumann");
    static bool ponctual_bool = (E_r_expr.find("ponctuel") != string::npos);

    if (neumann == 0){
        int k = cell_id(mesh->N, j, mesh->N+2, mesh->M+2);
        return E[k];
    }
    else if (ponctual_bool == true){
        static vector_1row start_end = parse_ponctual(E_r_expr);
        return ponctual_source(3, start_end[0], start_end[1], t, j);
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        double y = mesh->y[j];
        p.DefineVar("y", &y);
        if (first_call == 1){ p.SetExpr(E_r_expr); first_call = 0; }
        return p.Eval();
    }
}


vector_2d Solver::F_r(double t, int j){
    static int first_call = 1;
    static int neumann = F_r_x_expr.compare("neumann");

    if (neumann == 0){
        int k = cell_id(mesh->N, j, mesh->N+2, mesh->M+2);
        return F[k];
    }
    else{ 
        static Parser p1, p2;
        p1.DefineVar("t", &t); p2.DefineVar("t", &t);
        double y = mesh->y[j];
        p1.DefineVar("y", &y); p2.DefineVar("y", &y);
        if (first_call == 1){ p1.SetExpr(F_r_x_expr); p2.SetExpr(F_r_y_expr); first_call = 0; }
        return {p1.Eval(), p2.Eval()};
    }
}


double Solver::T_r(double t, int j){
    static int first_call = 1;
    static int neumann = T_r_expr.compare("neumann");

    if (neumann == 0){
        int k = cell_id(mesh->N, j, mesh->N+2, mesh->M+2);
        return T[k];
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        double y = mesh->y[j];
        p.DefineVar("y", &y);
        if (first_call == 1){ p.SetExpr(T_r_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::E_exact(double t, double x, double y){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("t", &t);
    p.DefineVar("x", &x);
    p.DefineVar("y", &y);
    p.DefineVar("t_0", &t_0);
    if (first_call == 1){ p.SetExpr(E_exact_expr); first_call = 0; }

    return p.Eval();
}


vector_2d Solver::F_exact(double t, double x, double y){
    static int first_call = 1;

    static Parser p1, p2;
    p1.DefineVar("t", &t); p2.DefineVar("t", &t);
    p1.DefineVar("x", &x); p2.DefineVar("x", &x);
    p1.DefineVar("y", &y); p2.DefineVar("y", &y);
    p1.DefineVar("t_0", &t_0); p2.DefineVar("t_0", &t_0);
    if (first_call == 1){ p1.SetExpr(F_exact_x_expr); p2.SetExpr(F_exact_y_expr); first_call = 0; }

    return {p1.Eval(), p2.Eval()};
}


double Solver::T_exact(double t, double x, double y){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("t", &t);
    p.DefineVar("x", &x);
    p.DefineVar("y", &y);
    p.DefineVar("t_0", &t_0);
    if (first_call == 1){ p.SetExpr(T_exact_expr); first_call = 0; }

    return p.Eval();
}


/**
 * Calcule sigma_{k, l}, moyenne de sigma_c entre les cellules k et l
 */ 
double compute_sigma(double sigma_k, double sigma_l){
    return 0.5 * (sigma_k + sigma_l);
}


/**
 * Calcule les M_{k, l}
 */ 
double compute_M(double dx, double sigma_kl){
    return 2 / (2 + dx * sigma_kl);
}


/**
 * Calcule les E_{k, l}n_{k, l}, flux de E entre les cellules k et l
 */ 
vector_2d flux_E(double l_kl, double M_kl, double E_k, double E_l, vector_2d F_k, vector_2d F_l, vector_2d n_kl){
    double tmp1 = 0.5 * (E_k + E_l);
    double tmp2 = 0.5 * (dot(F_l, n_kl) - dot(F_k, n_kl));
    return prod(l_kl*M_kl * (tmp1 - tmp2), n_kl);
}


/**
 * Calcule les F_{k, l} . n_{k, l}, flux de F entre les cellules k et l
 */ 
double flux_F(double l_kl, double M_kl, double E_k, double E_l, vector_2d F_k, vector_2d F_l, vector_2d n_kl){
    double tmp1 = 0.5 * (dot(F_k, n_kl) + dot(F_l, n_kl));
    double tmp2 = 0.5 * (E_l - E_k);
    return l_kl*M_kl * (tmp1 - tmp2);
}


void Solver::save_animation(int time_step){
    string file_name = "data/anim/animation." + to_string(time_step) + ".csv";
    ofstream file;
    file.open(file_name, ios_base::trunc);

    if(!file)
        throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_name + "'");

    file << "E,F_x,F_y,T,Tr\n";

    for (int j = 1; j < mesh->M+1; j++){
        for (int i = 1; i < mesh->N+1; i++){
            int k = cell_id(i, j, mesh->N+2, mesh->M+2);
            file << E[k] << "," << F[k][0] << "," << F[k][1] << "," << T[k] << "," << pow(E[k]/a, 0.25) << "\n";
        }
    }

    file.close();
}


void Solver::phase_1(){
    /* Des variables necessaires pour cette etape */
    double Theta;       // Theta = a*T^4 
    double E_n, T_n, Theta_n; 
    double E_next, Theta_next;


    omp_set_num_threads(1);
    // #pragma omp parallel for collapse(2) default(shared) firstprivate(var1,var2,sigma_a_func)
    #pragma omp parallel for collapse(2) default(shared)
    for (int i = 1; i < mesh->N+1; i++){
        for (int j = 1; j < mesh->M+1; j++){
    /*Remplacement de sigma_a et sigma_c en vue de la parallélisation*/
    mu::Parser sigma_a_func;
    double var1=0;
    double var2=0;
    sigma_a_func.DefineVar("rho", &var1); 
    sigma_a_func.DefineVar("T", &var2); 
    sigma_a_func.SetExpr(sigma_a_expr);

            int k = cell_id(i, j, mesh->N+2, mesh->M+2);
            // Initialisation etape 1 
            E_n = E[k];
            T_n = T[k];
            Theta_n = a * pow(T_n, 4);
            Theta = Theta_n;

            E_next = E[k];
            Theta_next = Theta;

            do{
                E[k] = E_next;
                Theta = Theta_next;
                
                T[k] = pow(Theta/a, 0.25);
                double mu_q = 1/ (pow(T_n, 3) + T_n*pow(T[k], 2) + T[k]*pow(T_n, 2) + pow(T[k], 3));  // + 1e-16);
                if (isnan(mu_q))
                    cerr << "ATTENTION! mu = nan" << " en k = " << k << endl;

                double rho_tmp = rho(i, j);
                // double sigma_a_tmp = sigma_a(rho_tmp, T[k]);
                var1 = rho_tmp; 
                var2 = T[k];
                double sigma_a_tmp = sigma_a_func.Eval();

                double tmp_1 = (1/dt) + c*sigma_a_tmp;
                double alpha = 1/dt/tmp_1;
                double beta = c*sigma_a_tmp/tmp_1;

                double tmp_2 = (rho_tmp*C_v*mu_q/dt) + c*sigma_a_tmp;
                double gamma = rho_tmp*C_v*mu_q/dt/tmp_2;
                double delta = c*sigma_a_tmp/tmp_2;

                E_next = (alpha*E_n + gamma*beta*Theta_n) / (1 - beta*delta);
                Theta_next = (gamma*Theta_n + alpha*delta*E_n) / (1 - beta*delta);

            } while (abs(E_next-E[k]) > precision && abs(Theta_next-Theta) > precision);
        }
    }
};


void Solver::phase_1_eq(){
    /* Des variables necessaires pour cette etape */
    double Theta;       // Theta = a*T^4
    double E_n, T_n, Theta_n;
    double E_next, Theta_next;


    for (int i = 1; i < mesh->N+1; i++){
        for (int j = 1; j < mesh->M+1; j++){
            int k = cell_id(i, j, mesh->N+2, mesh->M+2);

            // Initialisation
            E_n = E[k];
            T_n = T[k];
            Theta_n = a * pow(T_n, 4);
            Theta = Theta_n;

            E_next = E[k];
            Theta_next = Theta;
                // cout << "E = "<< E_next << " Theta = " << Theta_next<< endl;

            do{
                E[k] = E_next;
                Theta = Theta_next;

                T[k] = pow(Theta/a, 0.25);
                double mu_q = 1/ (pow(T_n, 3) + T_n*pow(T[k], 2) + T[k]*pow(T_n, 2) + pow(T[k], 3));
                if (isnan(mu_q))
                    // cerr << "ATTENTION! mu = nan" << " en k = " << k << endl;
                    ;

                double rho_tmp = rho(i, j);
                double alpha = c * sigma_a(rho_tmp, T[k]) * dt;
                double beta = rho_tmp * C_v * mu_q;

                double X_n = E_n - Theta_n;
                double Y_n = E_n + beta*Theta_n;

                double X = X_n / (1 + alpha*(1 + (1./beta)));
                double Y = Y_n;

                E_next = (beta*X + Y) / (1 + beta);
                Theta_next = (-X + Y) / (1 + beta);

            } while (abs(E_next-E[k]) > precision && abs(Theta_next-Theta) > precision);
        }
    }
};


void Solver::phase_2(){
    /* Vecteurs necessaires pour cette etape */
    vector_1row E_etoile(mesh->n_cells);
    vector_1row E_suiv(mesh->n_cells);
    vector_2rows F_etoile(mesh->n_cells, vector_2d(2));
    vector_2rows F_suiv(mesh->n_cells, vector_2d(2));

    /* Initialisation de l'etape */
    E_etoile = E;
    F_etoile = F;

    omp_set_num_threads(1);
    // #pragma omp parallel for collapse(2) firstprivate(var3,var4,sigma_c_func)
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < mesh->N+1; i++){
        for (int j = 1; j < mesh->M+1; j++){
    /*Remplacement de la fonction sigma_c en vue de la parallélisation*/
    mu::Parser sigma_c_func;
    double var3=0;
    double var4=0;
    sigma_c_func.DefineVar("rho", &var3); 
    sigma_c_func.DefineVar("T", &var4);
    sigma_c_func.SetExpr(sigma_c_expr);
            int k = cell_id(i, j, mesh->N+2, mesh->M+2);

            vector_2d sum_flux_E {0, 0};
            double sum_flux_F = 0;
            double sum_M_sigma = 0;
            vector_2d sum_l_M_n {0, 0};

            double rho_k = rho(i, j);
            vector_2d n_kl;                 // vecteur normal entre k et l
            vector_2d flux_E_kl;
            double flux_F_kl;

            for (int neighbor = 0; neighbor < 4; neighbor++){
                int l = mesh->neighb[k][neighbor];
                int i_prime = mesh->coord[l][0];
                int j_prime = mesh->coord[l][1];

                double rho_l = rho(i_prime, j_prime);
                // double sigma_kl = compute_sigma(sigma_c(rho_k, T[k]), sigma_c(rho_l, T[l]));     // Version parallélisée

                // Version parallélisée
                var3 = rho_k;
                var4 = T[k];
                double sigma1 = sigma_c_func.Eval(); // var3*var4
                var3 = rho_l;
                var4 = T[l];
                double sigma2 = sigma_c_func.Eval();
                double sigma_kl = compute_sigma(sigma1, sigma2);

                double M_kl = compute_M(mesh->dx, sigma_kl);

                double l_kl;
                if (neighbor == 0){                        // Voisin du haut
                    n_kl = {0, 1};
                    l_kl = mesh->dy;
                }
                else if (neighbor == 1){                        // Voisin du bas
                    n_kl = {0, -1};
                    l_kl = mesh->dy;
                }
                else if (neighbor == 2){                        // Voisin de gauche
                    n_kl = {-1, 0};
                    l_kl = mesh->dx;
                }
                else if (neighbor == 3){                        // Voisin de droite
                    n_kl = {1, 0};
                    l_kl = mesh->dx;
                }

                flux_E_kl = flux_E(l_kl, M_kl, E[k], E[l], F[k], F[l], n_kl);
                flux_F_kl = flux_F(l_kl, M_kl, E[k], E[l], F[k], F[l], n_kl);

                sum_flux_E = add(sum_flux_E, flux_E_kl);
                sum_flux_F += flux_F_kl;
                sum_M_sigma += (M_kl*sigma_kl);
                sum_l_M_n = add(sum_l_M_n, prod(l_kl*M_kl, n_kl));
            }

            double mes_omega = mesh->dx * mesh->dy;
            double tmp = (1./dt) + c*sum_M_sigma;
            // double tmp = (1./dt) + c*sum_M_sigma/4;         //************* ALTERNATIVE?
            double alpha = -c*dt / mes_omega;
            double beta = 1./dt / tmp;
            vector_2d gamma = prod(c/mes_omega / tmp, sum_l_M_n);
            double delta = -c/mes_omega / tmp;

            E_suiv[k] = E_etoile[k] + alpha*sum_flux_F;

            F_suiv[k] = add(add(prod(beta, F_etoile[k]), prod(E[k], gamma)), prod(delta, sum_flux_E));
        }
    }

    E = E_suiv;
    F = F_suiv;
};


void Solver::solve(){
    /* Initialisation de la doucle de resolution */
    for (int k = 0; k < mesh->n_cells; k++){
        int i = mesh->coord[k][0];
        int j = mesh->coord[k][1];
        double x_k = mesh->x[i];
        double y_k = mesh->y[j];

        E[k] = E_0(x_k, y_k);
        F[k] = F_0(x_k, y_k);
        T[k] = T_0(x_k, y_k);
    }

    /* Temps courant (translaté de t_0) et indice de l'iteration */
    double t = 0;
    int n = 0;

    /**
     * Boucle de resolution
     */
    while (t <= t_f){
        /* Enregistrement des signaux pour ce pas de temps */
        if (save_anim == true)
            save_animation(n);

        /* Affichage du progres */
        // cout << "  -- iteration " << n+1 << " sur " << step_count << " en cours ..." << endl;
        double progress = (n+1) * 100.0 / step_count;
        if (int(progress) % 5 == 0 && (progress - int(progress) < 1e-1) && int(progress) != 0)
            cout << "  -- " << setw(3) << int(progress) << " %" << endl;

        /* Signaux exportés */
        for (int i = 1; i < mesh->N+1; i++){
            int k = cell_id(i, mesh->M, mesh->N+2, mesh->M+2);
            E_up[n][i-1] = E[k];
            F_up[n][i-1] = l2_norm(F[k]);
            T_up[n][i-1] = T[k];

            k = cell_id(i, 1, mesh->N+2, mesh->M+2);
            E_down[n][i-1] = E[k];
            F_down[n][i-1] = l2_norm(F[k]);
            T_down[n][i-1] = T[k];
        }
        for (int j = 1; j < mesh->M+1; j++){
            int k = cell_id(1, j, mesh->N+2, mesh->M+2);
            E_left[n][j-1] = E[k];
            F_left[n][j-1] = l2_norm(F[k]);
            T_left[n][j-1] = T[k];

            k = cell_id(mesh->N, j, mesh->N+2, mesh->M+2);
            E_right[n][j-1] = E[k];
            F_right[n][j-1] = l2_norm(F[k]);
            T_right[n][j-1] = T[k];
        }

        /* *************** etape 1 ******************/
        phase_1();
        // phase_1_eq();

        /* Remplissage des mailles fantomes */
        for (int i = 1; i < mesh->N+1; i++){
            int k = cell_id(i, mesh->M+1, mesh->N+2, mesh->M+2);
            E[k] = E_u(t, i);
            F[k] = F_u(t, i);
            T[k] = T_u(t, i);

            k = cell_id(i, 0, mesh->N+2, mesh->M+2);
            E[k] = E_d(t, i);
            F[k] = F_d(t, i);
            T[k] = T_d(t, i);
        }
        for (int j = 1; j < mesh->M+1; j++){
            int k = cell_id(0, j, mesh->N+2, mesh->M+2);
            E[k] = E_l(t, j);
            F[k] = F_l(t, j);
            T[k] = T_l(t, j);

            k = cell_id(mesh->N+1, j, mesh->N+2, mesh->M+2);
            E[k] = E_r(t, j);
            F[k] = F_r(t, j);
            T[k] = T_r(t, j);
        }

        /* *************** etape 2 ******************/
        phase_2();

        time_steps[n] = t;
        t += dt;
        n += 1;
    }
};


void Solver::display(){
    cout << "-----------  E  -----------\n" ;
    for (int j = mesh->M; j > 0; j--){
        for (int i = 1; i < mesh->N+1; i++){
            int k = cell_id(i, j, mesh->N+2, mesh->M+2);
            cout << E[k] << "  ";
        }
        cout << "\n";
    }

    // cout << "-----------  F  -----------\n" ;
    // for (int j = mesh->M; j > 0; j--){
    //     for (int i = 1; i < mesh->N+1; i++){
    //         int k = cell_id(i, j, mesh->N+2, mesh->M+2);
    //         cout << l2_norm(F[k]) << "  ";
    //     }
    //     cout << "\n";
    // }

    // cout << "-----------  T  -----------\n" ;
    // for (int j = mesh->M; j > 0; j--){
    //     for (int i = 1; i < mesh->N+1; i++){
    //         int k = cell_id(i, j, mesh->N+2, mesh->M+2);
    //         cout << T[k] << "  ";
    //     }
    //     cout << "\n";
    // }

    cout << "\n";
};

Solver::~Solver(){
    free(E_up, step_count);
    free(F_up, step_count);
    free(T_up, step_count);

    free(E_down, step_count);
    free(F_down, step_count);
    free(T_down, step_count);

    free(E_left, step_count);
    free(F_left, step_count);
    free(T_left, step_count);

    free(E_right, step_count);
    free(F_right, step_count);
    free(T_right, step_count);

    if (rho_expr.find("crenau") != string::npos)
        free(attr, n_niche);
};
