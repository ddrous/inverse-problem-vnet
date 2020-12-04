#include <iostream>
#include <fstream>
#include <map>
#include <cmath>
#include <inttypes.h>

#include "include/exporter.hpp"

using namespace std;


Exporter::Exporter(const Solver *new_solver){
    if (new_solver == NULL)
        throw string("ERREUR: Pas de solveur a exporter");

    solver = new_solver;
}


/**
 * Ecris un signal spatial a n_rows lignes et n_cols colones dans le fichier file
 */
void write_spatial(ofstream &file, vector<double> const &signal, int n_rows, int n_cols){
    file << "\"[";
    for (int j = 0; j < n_cols; j++){
        file << "[";
        for (int i = 0; i < n_rows; i++){
            int k = cell_id(i+1, j+1, n_rows+2, n_cols+2);
            file << signal[k];
            if (i != n_rows-1) file << ", ";
        }
        file << "]";
        if (j != n_cols-1) file << ", ";
    }
    file << "]\"";
}


/**
 * Ecris un signal temporel a n_rows lignes et n_cols colones dans le fichier file
 */
void write_temporal(ofstream &file, double **signal, int n_rows, int n_cols){
    file << "\"[";
    for (int n = 0; n < n_rows; n++){
        file << "[";
        for (int i = 0; i < n_cols; i++){
            file << signal[n][i];
            if (i != n_cols-1) file << ", ";
        }
        file << "]";
        if (n != n_rows-1) file << ", ";
    }
    file << "]\"";
}


void Exporter::write_dataframe(std::string file_name, std::string mode){
    ofstream file;
    if (mode.compare("append") == 0)
        file.open(file_name, ios_base::app);
    else if (mode.compare("truncate") == 0)
        file.open(file_name, ios_base::trunc);
    else
        throw string ("ERREUR: Mode d'ouverture de fichier '" + mode + "' non reconnu");

    if(!file)
        throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_name + "'");

    //**************************** En mode trunc, il faut rajouter l'en tete
    if (mode.compare("append") != 0)
        file << "x_min,x_max,y_min,y_max,N,M,c,a,C_v,CFL,precision,t_0,t_f,step_count,rho_expr,sigma_a_expr,sigma_c_expr,E_0_expr,F_0_x_expr,F_0_y_expr,T_0_expr,E_u_expr,F_u_x_expr,F_u_y_expr,T_u_expr,E_d_expr,F_d_x_expr,F_d_y_expr,T_d_expr,E_l_expr,F_l_x_expr,F_l_y_expr,T_l_expr,E_r_expr,F_r_x_expr,F_r_y_expr,T_r_expr,rho_attr,rho,E_f,F_f_x,F_f_y,T_f,E_u,F_u,T_u,E_d,F_d,T_d,E_l,F_l,T_l,E_r,F_r,T_r" << "\n";

    const Solver *s = solver;           // Pour raccourcir

    //**************************** Parametres du probleme
    file << s->mesh->x_min << "," << s->mesh->x_max << "," << s->mesh->y_min << "," << s->mesh->y_max << "," << s->mesh->N << "," << s->mesh->M << "," << s->c << "," << s->a << "," << s->C_v << ","<< s->CFL << "," << s->precision << "," << s->t_0 << "," << s->t_f << "," << s->step_count << ",\"" << s->rho_expr << "\",\"" << s->sigma_a_expr << "\",\"" << s->sigma_c_expr << "\",\"" << s->E_0_expr << "\",\"" << s->F_0_x_expr << "\",\"" << s->F_0_y_expr << "\",\"" << s->T_0_expr << "\",\"" << s->E_u_expr << "\",\"" << s->F_u_x_expr << "\",\"" << s->F_u_y_expr << "\",\"" << s->T_u_expr << "\",\"" << s->E_d_expr << "\",\"" << s->F_d_x_expr << "\",\"" << s->F_d_y_expr << "\",\"" << s->T_d_expr << "\",\"" << s->E_l_expr << "\",\"" << s->F_l_x_expr << "\",\"" << s->F_l_y_expr << "\",\"" << s->T_l_expr << "\",\"" << s->E_r_expr << "\",\"" << s->F_r_x_expr << "\",\"" << s->F_r_y_expr << "\",\"" << s->T_r_expr << "\",";

    //************************** Attributs du crenau de densite
    file << "\"[";
    for (int l = 0; l < s->n_niche; l++){
        if (l != 0) file << ", ";
        file << "(" << s->attr[l][0] << ", " << s->attr[l][1] << ", " << s->attr[l][2] << ", " << s->attr[l][3] << ")";
    }
    file << "]\",";

    /* Donnes spatialles */
    write_spatial(file, s->rho_vec, s->mesh->N, s->mesh->M);    file << ",";

    write_spatial(file, s->E, s->mesh->N, s->mesh->M);    file << ",";

    vector<double> F_x(s->mesh->n_cells);
    vector<double> F_y(s->mesh->n_cells);
    for (int k = 0; k < s->mesh->n_cells; k++){
        F_x[k] = s->F[k][0];
        F_y[k] = s->F[k][1];
    }
    write_spatial(file, F_x, s->mesh->N, s->mesh->M);    file << ",";
    write_spatial(file, F_y, s->mesh->N, s->mesh->M);    file << ",";

    write_spatial(file, s->T, s->mesh->N, s->mesh->M);    file << ",";

    /* Donnees temporelles */
    write_temporal(file, s->E_up, s->step_count, s->mesh->N);    file << ",";
    write_temporal(file, s->F_up, s->step_count, s->mesh->N);    file << ",";
    write_temporal(file, s->T_up, s->step_count, s->mesh->N);    file << ",";

    write_temporal(file, s->E_down, s->step_count, s->mesh->N);    file << ",";
    write_temporal(file, s->F_down, s->step_count, s->mesh->N);    file << ",";
    write_temporal(file, s->T_down, s->step_count, s->mesh->N);    file << ",";

    write_temporal(file, s->E_left, s->step_count, s->mesh->M);    file << ",";
    write_temporal(file, s->F_left, s->step_count, s->mesh->M);    file << ",";
    write_temporal(file, s->T_left, s->step_count, s->mesh->M);    file << ",";

    write_temporal(file, s->E_right, s->step_count, s->mesh->M);    file << ",";
    write_temporal(file, s->F_right, s->step_count, s->mesh->M);    file << ",";
    write_temporal(file, s->T_right, s->step_count, s->mesh->M);    file << "\n";

    file.close();
}


/**
 * permet de localiser l'identifiant d'une position de source
 */
short int id_source_pos(string expr){
        double start = parse_ponctual(expr)[0];
        if (start == 0.1) return 0;
        else if (start == 0.3) return 1;
        else if (start == 0.5) return 2;
        else if (start == 0.7) return 3;
        else return -1;
}


/**
 * Ecris un signal spatial a n_rows lignes et n_cols colones dans le fichier file en binaire
 */
void write_spatial_signal_bin(ofstream &file, vector<double> const &signal, int n_rows, int n_cols){
    for (int j = 0; j < n_cols; j++){
        for (int i = 0; i < n_rows; i++){
            int k = cell_id(i+1, j+1, n_rows+2, n_cols+2);
            float val = float(signal[k]);
            file.write((char *)&val, sizeof(float));
        }
    }
}


/**
 * Ecris un signal temporel a n_rows lignes et n_cols colones dans le fichier file en binaire
 */
void write_temporal_signal_bin(ofstream &file, double **signal, int n_rows, int n_cols){
    for (int n = 0; n < n_rows; n++){
        for (int i = 0; i < n_cols; i++){
            float val = float(signal[n][i]);
            file.write((char *)&val, sizeof(float));
        }
    }
}


void Exporter::write_binary(std::string file_name, std::string mode, std::string n_simu){
    ofstream file;
    if (mode == "append")
        file.open(file_name, ios_base::app | ios_base::binary);
    else if (mode == "truncate")
        file.open(file_name, ios_base::trunc | ios_base::binary);
    else
        throw string ("ERREUR: Mode d'ouverture de fichier '" + mode + "' non reconnu. Modes possibles: 'append' ou 'truncate'");

    if(!file)
        throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_name + "'");

    const Solver *s = solver;           // Pour raccourcir

    /* Entete du fichier binaire */
    if (mode == "truncate") {
        /* Ecriture de la version du format sds */ 
        uint8_t magic[] = {0x73, 0x64, 0x73, 0x30, 0x32};       // sds version 02
        file.write((char *)magic, sizeof(magic));

        uint16_t n_simulations = atoi(n_simu.c_str());      // nombre de simulations dans ce fichier (4 en temps normal)
        file.write((char *)&n_simulations, sizeof(uint16_t));

        /* Ecriture des tailles de matrices */
        uint16_t N = s->mesh->N;
        uint16_t M = s->mesh->M;
        uint16_t step_count = s->step_count;
        file.write((char *)&N, sizeof(uint16_t));
        file.write((char *)&M, sizeof(uint16_t));
        file.write((char *)&step_count, sizeof(uint16_t));

        int8_t newline = 0xA;
        file.write((char *)&newline, sizeof(int8_t));
    }

    /* Localisation de la source */ 
    int8_t source_edge;      // le bord sur lequel se trouve la source
    int8_t source_pos;       // l'id de la position de la source sur son bord
    if (s->E_u_expr.find("ponctuel") != string::npos) {
        source_edge = 0;
        source_pos = id_source_pos(s->E_u_expr);
    }
    else if (s->E_d_expr.find("ponctuel") != string::npos) {
        source_edge = 1;
        source_pos = id_source_pos(s->E_d_expr);
    }
    else if (s->E_l_expr.find("ponctuel") != string::npos) {
        source_edge = 2;
        source_pos = id_source_pos(s->E_l_expr);
    }
    else if (s->E_r_expr.find("ponctuel") != string::npos) {
        source_edge = 3;
        source_pos = id_source_pos(s->E_r_expr);
    }
    else {
        source_edge = -1;
        source_pos = -1;
    }

    /* Ecriture de la source dans le fichier */
    file.write((char *)&source_edge, sizeof(int8_t));
    file.write((char *)&source_pos, sizeof(int8_t));

    /* Ecriture de la densite dans le fichier */
    write_spatial_signal_bin(file, s->rho_vec, s->mesh->M, s->mesh->N);
    
    /* Ecriture des matrices */
    write_temporal_signal_bin(file, s->E_up, s->step_count, s->mesh->N);
    write_temporal_signal_bin(file, s->F_up, s->step_count, s->mesh->N);
    write_temporal_signal_bin(file, s->T_up, s->step_count, s->mesh->N);

    write_temporal_signal_bin(file, s->E_down, s->step_count, s->mesh->N);
    write_temporal_signal_bin(file, s->F_down, s->step_count, s->mesh->N);
    write_temporal_signal_bin(file, s->T_down, s->step_count, s->mesh->N);

    write_temporal_signal_bin(file, s->E_left, s->step_count, s->mesh->M);
    write_temporal_signal_bin(file, s->F_left, s->step_count, s->mesh->M);
    write_temporal_signal_bin(file, s->T_left, s->step_count, s->mesh->M);

    write_temporal_signal_bin(file, s->E_right, s->step_count, s->mesh->M);
    write_temporal_signal_bin(file, s->F_right, s->step_count, s->mesh->M);
    write_temporal_signal_bin(file, s->T_right, s->step_count, s->mesh->M);

    int8_t newline = 0xA;
    file.write((char *)&newline, sizeof(int8_t));
    file.close();
}
