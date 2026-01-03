#include "Laplacian.h"
#include <iostream>
#include <vector>

using namespace std;
using namespace Eigen;

Laplacian::Laplacian(Function* function, DataFile* data_file, MACgrid* grid) :
    _fct(function), _df(data_file), _grid(grid)
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    
    // Prédimensionnement de la matrice
    _H.resize(Nx * Ny, Nx * Ny);
}

void Laplacian::BuildMatrix()
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();
    
    // Coefficients de discrétisation 1/h^2
    double Cx = 1.0 / (hx * hx);
    double Cy = 1.0 / (hy * hy);

    // Utilisation de liste de triplets pour une construction efficace
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(5 * Nx * Ny); // Estimation : 5 coefficients par ligne (stencil à 5 points)

    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            int k = _grid->GetPIndex(i, j); // Index global courant
            double diag = 0.0; // Valeur de la diagonale

            // --- Voisins en X (Ouest / Est) ---
            if (j > 0) { // Voisin Gauche
                int k_W = _grid->GetPIndex(i, j - 1);
                triplets.push_back(Eigen::Triplet<double>(k, k_W, -Cx));
                diag += Cx;
            } 
            // else: Bord Gauche -> Neumann (dP/dn = 0) => pas de flux sortant

            if (j < Nx - 1) { // Voisin Droite
                int k_E = _grid->GetPIndex(i, j + 1);
                triplets.push_back(Eigen::Triplet<double>(k, k_E, -Cx));
                diag += Cx;
            }
            // else: Bord Droite -> Neumann

            // --- Voisins en Y (Sud / Nord) ---
            if (i > 0) { // Voisin Bas
                int k_S = _grid->GetPIndex(i - 1, j);
                triplets.push_back(Eigen::Triplet<double>(k, k_S, -Cy));
                diag += Cy;
            } 
            // else: Bord Bas -> Neumann

            if (i < Ny - 1) { // Voisin Haut
                int k_N = _grid->GetPIndex(i + 1, j);
                triplets.push_back(Eigen::Triplet<double>(k, k_N, -Cy));
                diag += Cy;
            }
            // else: Bord Haut -> Neumann

            // Remplissage de la diagonale
            // Note : Pour Neumann pur, la somme des coefs d'une ligne est nulle.
            triplets.push_back(Eigen::Triplet<double>(k, k, diag));
        }
    }
    
    // --- Fixation du mode hydrostatique ---
    // Avec des CL Neumann partout, la matrice est singulière (infinité de solutions à une constante près).
    // On impose arbitrairement P(0,0) = 0 en pénalisant fortement la diagonale du premier élément.
    int k_ref = _grid->GetPIndex(0, 0);
    triplets.push_back(Eigen::Triplet<double>(k_ref, k_ref, 1.0e9)); 

    // Assemblage final et factorisation
    _H.setFromTriplets(triplets.begin(), triplets.end());
    _solver.compute(_H);
    
    if(_solver.info() != Eigen::Success) {
        cerr << "Erreur: Factorisation de Cholesky a échoué." << endl;
        exit(1);
    }
}

Eigen::VectorXd Laplacian::ComputeDivergence(const Eigen::VectorXd& U, const Eigen::VectorXd& V)
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();

    Eigen::VectorXd div(Nx * Ny);
    
    // Parcours des cellules de pression
    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            int k = _grid->GetPIndex(i, j);

            // Indices des vitesses aux faces entourant la cellule (i,j)
            int k_u_E = _grid->GetUIndex(i, j + 1); // Est
            int k_u_W = _grid->GetUIndex(i, j);     // Ouest
            int k_v_N = _grid->GetVIndex(i + 1, j); // Nord
            int k_v_S = _grid->GetVIndex(i, j);     // Sud

            // Divergence discrète : (du/dx + dv/dy)
            double du_dx = (U(k_u_E) - U(k_u_W)) / hx;
            double dv_dy = (V(k_v_N) - V(k_v_S)) / hy;

            div(k) = du_dx + dv_dy;
        }
    }
    return div;
}

void Laplacian::Solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& p_sol)
{
    // Résolution du système : H * P = -RHS
    // (Le signe moins vient de l'équation de projection : Lap(P) = div(U*) / dt)
    // Ici on suppose que 'rhs' contient div(U*)/dt.
    
    p_sol = _solver.solve(-rhs);
    
    if(_solver.info() != Eigen::Success) {
        cerr << "Erreur: Résolution du système linéaire a échoué." << endl;
    }
}

void Laplacian::ComputeGradient(const Eigen::VectorXd& p, Eigen::VectorXd& gradPx, Eigen::VectorXd& gradPy)
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();

    // Redimensionnement
    gradPx.resize((Nx + 1) * Ny);
    gradPy.resize(Nx * (Ny + 1));
    gradPx.setZero();
    gradPy.setZero();

    // 1. Gradient X (défini aux faces U)
    // On ignore les faces de bord (j=0 et j=Nx) car on y impose souvent Dirichlet pour U
    for (int i = 0; i < Ny; ++i) {
        for (int j = 1; j < Nx; ++j) {
            int k_u = _grid->GetUIndex(i, j);
            int k_p_E = _grid->GetPIndex(i, j);     // Pression à droite
            int k_p_W = _grid->GetPIndex(i, j - 1); // Pression à gauche

            gradPx(k_u) = (p(k_p_E) - p(k_p_W)) / hx;
        }
    }

    // 2. Gradient Y (défini aux faces V)
    // On ignore les faces de bord (i=0 et i=Ny)
    for (int i = 1; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            int k_v = _grid->GetVIndex(i, j);
            int k_p_N = _grid->GetPIndex(i, j);     // Pression en haut
            int k_p_S = _grid->GetPIndex(i - 1, j); // Pression en bas

            gradPy(k_v) = (p(k_p_N) - p(k_p_S)) / hy;
        }
    }
}