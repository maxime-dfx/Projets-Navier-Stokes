// ====================================================================================
//                                  VTKWRITER.CPP
// ====================================================================================
// Implémentation de l'export VTK.
// Correction : Ajout des includes manquants pour DataFile et MACgrid.
// ====================================================================================

#include "VTKWriter.h"
#include "DataFile.h"
#include "MACgrid.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem> 
#include <algorithm> // pour std::max, std::min

using namespace std;
namespace fs = std::filesystem;

VTKWriter::VTKWriter(DataFile* df, MACgrid* grid, string sim_name)
    : _df(df), _grid(grid), _sim_name(sim_name) 
{
    // --- 1. Gestion robuste des dossiers ---
    string root_dir = _df->Get_results();
    fs::path output_path = fs::path(root_dir) / _sim_name / "VTK";

    try {
        if (!fs::exists(output_path)) {
            fs::create_directories(output_path);
            cout << ">> [VTK] Dossier cree : " << output_path << endl;
        }
    } catch (const fs::filesystem_error& e) {
        cerr << "ERREUR CRITIQUE : Impossible de creer le dossier " << output_path << endl;
        cerr << "Raison : " << e.what() << endl;
        exit(EXIT_FAILURE);
    }
}

void VTKWriter::Write(int iteration, double t) {
    // --- 2. Nom du fichier ---
    string root_dir = _df->Get_results();
    fs::path output_dir = fs::path(root_dir) / _sim_name / "VTK";
    
    stringstream ss;
    ss << output_dir.string() << "/" << _sim_name << "_" 
       << setfill('0') << setw(5) << iteration << ".vtk";
       
    ofstream file(ss.str());
    if (!file.is_open()) {
        cerr << "ERREUR : Impossible d'ecrire le fichier " << ss.str() << endl;
        return;
    }

    // --- 3. Paramètres Grille ---
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();
    double xmin = _df->Get_xmin();
    double ymin = _df->Get_ymin();

    // En-tête VTK (Legacy ASCII)
    file << "# vtk DataFile Version 3.0" << endl;
    file << "Navier-Stokes 2D - t=" << t << endl;
    file << "ASCII" << endl;
    file << "DATASET STRUCTURED_POINTS" << endl;
    file << "DIMENSIONS " << Nx << " " << Ny << " 1" << endl;
    // On centre les points VTK sur les centres de pression
    file << "ORIGIN " << xmin + hx/2.0 << " " << ymin + hy/2.0 << " 0" << endl;
    file << "SPACING " << hx << " " << hy << " 1" << endl;
    file << "POINT_DATA " << Nx * Ny << endl;

    // --- 4. CHAMP 1 : PRESSION (Scalaire) ---
    file << "SCALARS pressure double" << endl;
    file << "LOOKUP_TABLE default" << endl;
    const Eigen::VectorXd& P = _grid->GetP();
    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            file << P(_grid->GetPIndex(i, j)) << " ";
        }
        file << endl;
    }

    // --- 5. CHAMP 2 : VITESSE (Vecteur) ---
    // Interpolation au centre des mailles
    file << "VECTORS velocity double" << endl;
    const Eigen::VectorXd& U = _grid->GetU();
    const Eigen::VectorXd& V = _grid->GetV();

    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            // Moyenne des faces gauche/droite pour U
            double u_val = 0.5 * ( U(_grid->GetUIndex(i, j)) + U(_grid->GetUIndex(i, j+1)) );
            // Moyenne des faces bas/haut pour V
            double v_val = 0.5 * ( V(_grid->GetVIndex(i, j)) + V(_grid->GetVIndex(i+1, j)) );
            
            file << u_val << " " << v_val << " 0.0" << endl;
        }
    }

    // --- 6. CHAMP 3 : VORTICITÉ (Scalaire - AJOUTÉ) ---
    // w = dv/dx - du/dy
    file << "SCALARS vorticity double" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            // Calcul par différences finies centrées simples sur les valeurs interpolées
            // dv/dx
            double v_E = 0.5 * (V(_grid->GetVIndex(i, min(j+1, Nx-1))) + V(_grid->GetVIndex(i+1, min(j+1, Nx-1))));
            double v_W = 0.5 * (V(_grid->GetVIndex(i, max(j-1, 0)))    + V(_grid->GetVIndex(i+1, max(j-1, 0))));
            double dv_dx = (v_E - v_W) / (2.0 * hx);

            // du/dy
            double u_N = 0.5 * (U(_grid->GetUIndex(min(i+1, Ny-1), j)) + U(_grid->GetUIndex(min(i+1, Ny-1), j+1)));
            double u_S = 0.5 * (U(_grid->GetUIndex(max(i-1, 0), j))    + U(_grid->GetUIndex(max(i-1, 0), j+1)));
            double du_dy = (u_N - u_S) / (2.0 * hy);

            file << (dv_dx - du_dy) << " ";
        }
        file << endl;
    }

    file.close();
}