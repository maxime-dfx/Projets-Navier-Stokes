#include "ErrorAnalysis.h"
#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace Eigen;

ErrorAnalysis::ErrorAnalysis(DataFile* df, MACgrid* grid) : _df(df), _grid(grid) {}

ErrorAnalysis::~ErrorAnalysis() {}

double ErrorAnalysis::ComputePoiseuilleErrorL2() {
    double error_L2 = 0.0;
    
    // Récupération géométrie
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();
    double ymin = _df->Get_ymin();
    double ymax = _df->Get_ymax();
    double H = ymax - ymin; // Hauteur du canal

    // Paramètres de l'écoulement (suppose U_max = 1.0 ou défini par BC left)
    // Pour un Poiseuille standard avec BC gauche = 1.0 (profil parabolique max)
    double U_max = _df->Get_BC_Left_dir();
    const VectorXd& U = _grid->GetU();

    // On parcourt la grille fluide
    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j <= Nx; ++j) { // <--- MODIFICATION ICI (<= Nx)
            // Ignorer les obstacles solides
            if (_grid->IsSolidU(i, j)) continue;

            // Coordonnée Y au centre de la face U
            // Rappel : U est décalé en X mais centré en Y par rapport à la maille
            double y = ymin + (i + 0.5) * hy; 

            // Solution Exacte Poiseuille : U(y) = 4 * Umax * Y_adim * (1 - Y_adim)
            // Y_adim varie de 0 à 1
            double Y_adim = (y - ymin) / H;
            double u_exact = 4.0 * U_max * Y_adim * (1.0 - Y_adim);

            // Valeur Numérique
            double u_num = U(_grid->GetUIndex(i, j));

            // Somme quadratique pondérée par le volume (surface en 2D)
            error_L2 += pow(u_num - u_exact, 2) * hx * hy;
        }
    }

    return sqrt(error_L2);
}