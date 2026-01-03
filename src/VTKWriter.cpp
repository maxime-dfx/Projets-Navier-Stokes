#include "VTKWriter.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem> // Nécessite C++17 (flag -std=c++17 dans le Makefile)

using namespace std;
namespace fs = std::filesystem; // Alias pour simplifier le code

VTKWriter::VTKWriter(DataFile* df, MACgrid* grid, string sim_name)
    : _df(df), _grid(grid), _sim_name(sim_name) 
{
    // --- 1. Gestion des Dossiers ---
    
    // Récupère le dossier racine défini dans le TOML (ex: "results")
    string root_dir = _df->Get_results();
    
    // Construit le chemin complet : results/MaSimu/VTK
    // L'opérateur "/" gère les séparateurs OS (Windows/Linux)
    fs::path output_path = fs::path(root_dir) / _sim_name / "VTK";

    // Crée les dossiers récursivement s'ils n'existent pas
    try {
        if (!fs::exists(output_path)) {
            fs::create_directories(output_path);
            cout << ">> VTKWriter : Dossier cree -> " << output_path << endl;
        }
    } catch (const fs::filesystem_error& e) {
        cerr << "Erreur Critique VTKWriter : Impossible de creer le dossier " << output_path << endl;
        cerr << e.what() << endl;
        exit(1);
    }
}

void VTKWriter::Write(int iteration, double t) {
    // --- 2. Construction du Nom de Fichier ---
    string root_dir = _df->Get_results();
    fs::path output_dir = fs::path(root_dir) / _sim_name / "VTK";
    
    stringstream ss;
    // Format : NomSimu_00123.vtk
    ss << output_dir.string() << "/" << _sim_name << "_" 
       << setfill('0') << setw(5) << iteration << ".vtk";
       
    string filename = ss.str();
    ofstream file(filename);

    if (!file.is_open()) {
        cerr << "Erreur VTKWriter : Impossible d'ecrire dans " << filename << endl;
        return;
    }

    // --- 3. Ecriture de l'En-tête VTK ---
    int Nx = _df->Get_Nx(); // Attention : Vérifier si GetNx() ou Get_Nx() dans MACgrid.h
    int Ny = _df->Get_Ny(); // On suppose ici que MACgrid a ces accesseurs

    double hx = _df->Get_hx();
    double hy = _df->Get_hy();
    double xmin = _df->Get_xmin();
    double ymin = _df->Get_ymin();

    file << "# vtk DataFile Version 3.0" << endl;
    file << "Navier-Stokes 2D - t = " << t << endl;
    file << "ASCII" << endl;
    file << "DATASET STRUCTURED_POINTS" << endl;
    file << "DIMENSIONS " << Nx << " " << Ny << " 1" << endl; // Grille centrée (Cell centers -> Points dans VTK Structured)
    file << "ORIGIN " << xmin + hx/2.0 << " " << ymin + hy/2.0 << " 0" << endl; // Décalage pour centrer les points
    file << "SPACING " << hx << " " << hy << " 1" << endl;
    file << "POINT_DATA " << Nx * Ny << endl;

    // --- 4. Ecriture de la Pression (Scalaire) ---
    file << "SCALARS pressure double" << endl;
    file << "LOOKUP_TABLE default" << endl;
    
    const Eigen::VectorXd& P = _grid->GetP();

    // Boucle sur les lignes (Y) puis colonnes (X)
    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            // P est définie au centre (i,j)
            file << P(_grid->GetPIndex(i, j)) << " ";
        }
        file << endl;
    }

    // --- 5. Ecriture de la Vitesse (Vecteur) ---
    file << "VECTORS velocity double" << endl;
    
    const Eigen::VectorXd& U = _grid->GetU();
    const Eigen::VectorXd& V = _grid->GetV();

    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            // Interpolation au centre de la maille (i,j)
            
            // U est sur les faces verticales
            // Index (i, j) = Face Gauche
            // Index (i, j+1) = Face Droite
            double u_val = 0.5 * ( U(_grid->GetUIndex(i, j)) + U(_grid->GetUIndex(i, j+1)) );

            // V est sur les faces horizontales
            // Index (i, j) = Face Bas
            // Index (i+1, j) = Face Haut
            double v_val = 0.5 * ( V(_grid->GetVIndex(i, j)) + V(_grid->GetVIndex(i+1, j)) );

            // Ecriture (u, v, w=0)
            file << u_val << " " << v_val << " 0.0" << endl;
        }
    }

    file.close();
}