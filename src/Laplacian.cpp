// ====================================================================================
//                                 LAPLACIAN.CPP
// ====================================================================================
// Implémentation du solveur de Poisson.
// CORRECTION : Redimensionnement explicite des vecteurs gradients pour éviter le crash.
// ====================================================================================

#include "Laplacian.h"
#include "Function.h"
#include "DataFile.h"
#include "MACgrid.h"

#include <iostream>
#include <vector>
#include <cmath>

// ============================================================================
// CONSTRUCTEUR
// ============================================================================
Laplacian::Laplacian(Function* fct, DataFile* df, MACgrid* grid)
    : _fct(fct), _df(df), _grid(grid) 
{
    // Pré-dimensionnement pour éviter les réallocations dynamiques
    int N = _df->Get_Nx() * _df->Get_Ny();
    _H.resize(N, N);
}

// ============================================================================
// CONSTRUCTION DE LA MATRICE (BUILD MATRIX)
// ============================================================================
void Laplacian::BuildMatrix() 
{
    std::cout << ">> [Laplacian] Construction de la matrice..." << std::endl;

    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();

    // Coefficients de discrétisation FD (1/h^2)
    double Cx = 1.0 / (hx * hx);
    double Cy = 1.0 / (hy * hy);

    // Liste de triplets pour l'insertion efficace dans Eigen::SparseMatrix
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(5 * Nx * Ny); 

    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            int k = _grid->GetPIndex(i, j); // Index global (Ligne de la matrice)
            double diag = 0.0; 

            // --- VOISIN OUEST (j-1) ---
            if (j > 0) { 
                int k_W = _grid->GetPIndex(i, j - 1);
                triplets.emplace_back(k, k_W, -Cx);
                diag += Cx;
            } 
            // else: Bord Ouest -> Neumann homogène (dP/dn = 0) => Flux nul

            // --- VOISIN EST (j+1) ---
            if (j < Nx - 1) { 
                int k_E = _grid->GetPIndex(i, j + 1);
                triplets.emplace_back(k, k_E, -Cx);
                diag += Cx;
            }

            // --- VOISIN SUD (i-1) ---
            if (i > 0) { 
                int k_S = _grid->GetPIndex(i - 1, j);
                triplets.emplace_back(k, k_S, -Cy);
                diag += Cy;
            }

            // --- VOISIN NORD (i+1) ---
            if (i < Ny - 1) { 
                int k_N = _grid->GetPIndex(i + 1, j);
                triplets.emplace_back(k, k_N, -Cy);
                diag += Cy;
            }

            // --- DIAGONALE ---
            triplets.emplace_back(k, k, diag);
        }
    }

    // --- GESTION DE LA SINGULARITÉ (Fixation du mode constant) ---
    int k_ref = _grid->GetPIndex(0, 0);
    triplets.emplace_back(k_ref, k_ref, 1.0e9);

    // Assemblage final
    _H.setFromTriplets(triplets.begin(), triplets.end());

    // Factorisation (LLT)
    _solver.compute(_H);

    if (_solver.info() != Eigen::Success) {
        std::cerr << "ERREUR CRITIQUE : Échec de la factorisation Cholesky !" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// ============================================================================
// CALCUL DE LA DIVERGENCE (COMPUTE DIVERGENCE)
// ============================================================================

void Laplacian::ComputeDivergence(const Eigen::VectorXd& U, const Eigen::VectorXd& V, Eigen::VectorXd& div_out) 
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    double inv_hx = 1.0 / _df->Get_hx();
    double inv_hy = 1.0 / _df->Get_hy();

    // On s'assure que le vecteur de sortie a la bonne taille (sécurité)
    // Comme on l'a pré-alloué dans TimeScheme, cela ne coûtera rien (pas de réallocation).
    if (div_out.size() != Nx * Ny) {
        div_out.resize(Nx * Ny);
    }

    // Parcours de toutes les cellules de pression
    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            int k = _grid->GetPIndex(i, j);

            // Récupération des indices de faces (Staggered Grid)
            int k_u_E = _grid->GetUIndex(i, j + 1);
            int k_u_W = _grid->GetUIndex(i, j);
            int k_v_N = _grid->GetVIndex(i + 1, j);
            int k_v_S = _grid->GetVIndex(i, j);

            // Calcul Divergence discrète : (du/dx + dv/dy)
            double du_dx = (U(k_u_E) - U(k_u_W)) * inv_hx;
            double dv_dy = (V(k_v_N) - V(k_v_S)) * inv_hy;

            // Écriture directe dans le buffer fourni
            div_out(k) = du_dx + dv_dy;
        }
    }
}

// ============================================================================
// RÉSOLUTION DU SYSTÈME (SOLVE)
// ============================================================================
void Laplacian::Solve(const Eigen::VectorXd& rhs_in, Eigen::VectorXd& p_sol) 
{
    // L'équation de projection est : Lap(P) = (rho/dt) * div(u*)
    Eigen::VectorXd rhs = -rhs_in; 

    // CONDITION DE FREDHOLM (Solvabilité pour Neumann)
    double mean_val = rhs.mean();
    rhs.array() -= mean_val;

    // Résolution
    p_sol = _solver.solve(rhs);

    // Centrage de la pression
    double p_mean = p_sol.mean();
    p_sol.array() -= p_mean;

    if (_solver.info() != Eigen::Success) {
        std::cerr << "ERREUR : Le solveur linéaire a échoué à converger." << std::endl;
    }
}

// ============================================================================
// CALCUL DU GRADIENT (COMPUTE GRADIENT)
// ============================================================================
void Laplacian::ComputeGradient(const Eigen::VectorXd& p, Eigen::VectorXd& gradPx, Eigen::VectorXd& gradPy) 
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    double inv_hx = 1.0 / _df->Get_hx();
    double inv_hy = 1.0 / _df->Get_hy();

    // On redimensionne les vecteurs avant de les remplir !
    // Sinon, l'accès par () cause un crash "Assertion failed".
    gradPx.setZero(_grid->GetU().size());
    gradPy.setZero(_grid->GetV().size());

    // --- GRADIENT X (Défini sur les faces U verticales) ---
    for (int i = 0; i < Ny; ++i) {
        for (int j = 1; j < Nx; ++j) {
            int k_u   = _grid->GetUIndex(i, j);
            int k_p_E = _grid->GetPIndex(i, j);     // Cellule Droite
            int k_p_W = _grid->GetPIndex(i, j - 1); // Cellule Gauche
            
            gradPx(k_u) = (p(k_p_E) - p(k_p_W)) * inv_hx;
        }
    }

    // --- GRADIENT Y (Défini sur les faces V horizontales) ---
    for (int i = 1; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            int k_v   = _grid->GetVIndex(i, j);
            int k_p_N = _grid->GetPIndex(i, j);     // Cellule Haut
            int k_p_S = _grid->GetPIndex(i - 1, j); // Cellule Bas

            gradPy(k_v) = (p(k_p_N) - p(k_p_S)) * inv_hy;
        }
    }
}