// ====================================================================================
//                                ERROR_ANALYSIS.CPP
// ====================================================================================
// Implémentation du calcul d'erreur.
// ====================================================================================

#include "ErrorAnalysis.h"
#include <cmath>
#include <iostream>

ErrorAnalysis::ErrorAnalysis(DataFile* df, MACgrid* grid) 
    : _df(df), _grid(grid) 
{}

double ErrorAnalysis::ComputePoiseuilleErrorL2() {
    double error_L2_sq = 0.0;
    
    // Géométrie
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();
    double ymin = _df->Get_ymin();
    double ymax = _df->Get_ymax();
    double H = ymax - ymin;

    // Paramètre analytique
    double U_max = _df->Get_BC_Left_dir(); // Vitesse max au centre
    const Eigen::VectorXd& U = _grid->GetU();

    // Parcours du domaine Fluide
    // On ignore les ghost cells (i=0..Ny-1, j=1..Nx-1 pour U interne)
    // Mais pour l'erreur globale, on peut sommer sur tout le domaine interne.
    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j <= Nx; ++j) { // U est défini sur les faces verticales
            
            // Si c'est un solide (obstacle), on ne compte pas l'erreur
            if (_grid->IsSolidU(i, j)) continue;

            // Coordonnée Y du centre de la face U
            // U(i,j) est en x = xmin + j*hx, y = ymin + (i+0.5)*hy
            double y = ymin + (i + 0.5) * hy;

            // Solution Exacte (Parabole)
            double Y_adim = (y - ymin) / H;
            double u_exact = 4.0 * U_max * Y_adim * (1.0 - Y_adim);

            // Solution Numérique
            double u_num = U(_grid->GetUIndex(i, j));

            // Accumulation quadratique (Intégrale ~ Somme * surface)
            double diff = u_num - u_exact;
            error_L2_sq += diff * diff * (hx * hy);
        }
    }

    return std::sqrt(error_L2_sq);
}