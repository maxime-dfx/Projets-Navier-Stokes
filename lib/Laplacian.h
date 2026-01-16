// ====================================================================================
//                                 LAPLACIAN.H
// ====================================================================================
// Description : Résolution de l'équation de Poisson pour la pression.
//               Div(Grad P) = RHS
//               Gère la matrice creuse et le solveur direct Cholesky.
// ====================================================================================

#ifndef _LAPLACIAN_H_
#define _LAPLACIAN_H_

#include <Eigen/Sparse>
#include <Eigen/Dense>

// Forward declarations pour limiter les dépendances dans le header
class Function;
class DataFile;
class MACgrid;

class Laplacian {
private:
    Function* _fct;
    DataFile* _df;
    MACgrid* _grid;

    // Matrice du Laplacien (Stockée pour éviter la reconstruction à chaque pas)
    Eigen::SparseMatrix<double> _H;

    // Solveur direct optimisé pour matrices Symétriques Définies Positives (SPD)
    // SimplicialLLT est plus rapide que LU pour ce type de problème
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> _solver;

public:
    // ============================================================================
    // CONSTRUCTEUR & DESTRUCTEUR
    // ============================================================================
    Laplacian(Function* fct, DataFile* df, MACgrid* grid);
    ~Laplacian() = default;

    // ============================================================================
    // MÉTHODES PRINCIPALES
    // ============================================================================
    
    // Construit la matrice du Laplacien (Stencil à 5 points)
    // Applique les CL de Neumann et gère la singularité via pénalisation.
    void BuildMatrix();

    // Calcule la divergence du champ de vitesse intermédiaire (u*, v*)
    // C'est le second membre de l'équation de Poisson (RHS)
    void ComputeDivergence(const Eigen::VectorXd& U, const Eigen::VectorXd& V, Eigen::VectorXd& div_out);
        
    // Résout le système linéaire H * P = RHS
    // Applique la correction de compatibilité (Condition de Fredholm)
    void Solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& p_sol);

    // Calcule le gradient de pression pour l'étape de correction
    void ComputeGradient(const Eigen::VectorXd& p, Eigen::VectorXd& gradPx, Eigen::VectorXd& gradPy);
};

#endif // _LAPLACIAN_H_