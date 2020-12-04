#include <iostream>
#include <fstream>

#include "include/config.hpp"

using namespace std;


Config::Config(string file_path){
    file_name = file_path;

    values["x_min"] = "";
    values["x_max"] = "";
    values["y_min"] = "";
    values["y_max"] = "";
    values["N"] = "";
    values["M"] = "";

    values["c"] = "";
    values["a"] = "";
    values["C_v"] = "";

    values["CFL"] = "";
    values["precision"] = "";
    values["t_0"] = "";
    values["t_f"] = "";

    values["rho"] = "";
    values["sigma_a"] = "";
    values["sigma_c"] = "";

    values["E_exact"] = "";
    values["F_exact_x"] = "";
    values["F_exact_y"] = "";
    values["T_exact"] = "";

    values["E_0"] = "";
    values["F_0_x"] = "";
    values["F_0_y"] = "";
    values["T_0"] = "";

    values["E_u"] = "";
    values["F_u_x"] = "";
    values["F_u_y"] = "";
    values["T_u"] = "";

    values["E_d"] = "";
    values["F_d_x"] = "";
    values["F_d_y"] = "";
    values["T_d"] = "";

    values["E_l"] = "";
    values["F_l_x"] = "";
    values["F_l_y"] = "";
    values["T_l"] = "";

    values["E_r"] = "";
    values["F_r_x"] = "";
    values["F_r_y"] = "";
    values["T_r"] = "";

    values["export_file"] = "";
    values["export_mode"] = "";
    values["write_mode"] = "";
    values["simu_count"] = "";

    size = values.size();
}


void Config::read(){
    map<string, string> :: iterator it;     // parcours des valeurs

    map<string, int> count;                 // compte les occurence de chaque parametre dans le fichier
    for (it = values.begin(); it != values.end(); ++it)
        count[it->first] = 0;

    string name;                        // nom du parametre en cours de traitement

    bool read_unkown = false;           // lecture d'un parametre inconnue
    string unkown_name;                 // nom de l'inconnu

    ifstream file(file_name);
    if(file){
        while(!file.eof()){
            file >> name;
            it = values.find(name);
            if (it != values.end()){
                count[it->first] ++ ;
                file >> values[it->first];
                cout << "   -- " << name << " : " << values[it->first] << endl;
            } 
            else{
                read_unkown = true;
                unkown_name = name;
                break;
            }
        }
        
        file.close();
    
    } else
        throw string ("ERREUR: Erreur d'ouverture du fichier de configuration '" + file_name + "'");

    if (read_unkown == true)
        throw string ("ERREUR: Parametre inconnu '" + unkown_name + "' dans le fichier de configuration");

    n_param = 0;        // Décompte du nombre de paramètres lus
    for (it = values.begin(); it != values.end(); ++it){
        if (count[it->first] < 1 && !(it->first == "E_exact") && !(it->first == "F_exact_x") && !(it->first == "F_exact_y") && !(it->first == "T_exact") && !(it->first == "simu_count") && !(it->first == "y_min") && !(it->first == "y_max") && !(it->first == "N") && !(it->first == "M")){
            throw string ("ERREUR: Paramètre '"+it->first+"' manquant dans le fichier de configuration");}
        else if (count[it->first] > 1)
            throw string ("ERREUR: Paramètre '"+it->first+"' dupliqué dans le fichier de configuration");
        else if (count[it->first] == 1)
            n_param++;
    } 

}


void Config::check_numpy(){

    // Traitement du cas special numpy. M et N sont donnés, y_max s'obtient pour un maillage uniforme
    if  (values["rho"].find("numpy") != string::npos){

        cnpy::NpyArray npydata = cnpy::npy_load(values["rho"].substr(6));

        int N = npydata.shape[1];         // Nombre de divisions suivant l'abcisse
        int M = npydata.shape[0];         // Nombre de divisions suivant l'ordonée
        values["N"] = to_string(N);
        values["M"] = to_string(M);

        // Le problème est traité sur le cercle unité
        values["x_min"] = to_string(0);
        values["x_max"] = to_string(1);
        values["y_min"] = to_string(0);
        double x_min = atof(values["x_min"].c_str());
        double x_max = atof(values["x_max"].c_str());
        // double y_min = atof(values["y_min"].c_str());
        
        double y_max = M * (x_max - x_min) / N;
        values["y_min"] = to_string(0);
        values["y_max"] = to_string(y_max);
    }   

    // This is the normal case where the density is computed on the uniform grid.
    else{
        double x_min = atof(values["x_min"].c_str());
        double x_max = atof(values["x_max"].c_str());
        double y_min = atof(values["y_min"].c_str());
        double y_max = atof(values["y_max"].c_str());
        int N = atoi(values["N"].c_str());

        int M = int(N*(y_max-y_min) / (x_max-x_min));
        values["M"] = to_string(M);
    }

}


Config::~Config(){}