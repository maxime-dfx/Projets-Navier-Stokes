#ifndef _LAPLACIAN_H
#define _LAPLACIAN_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Function.h"
#include "MACgrid.h"

/**
 * @brief Classe gérant la résolution de l'équation de Poisson pour la pression.
 * Discrétisation par Différences Finies sur grille MAC décalée.
 */
class Laplacian
{
private:
    Function* _fct;
    DataFile* _df;
    MACgrid* _grid;
    
    // Matrice du Laplacien (Sparse)
    Eigen::SparseMatrix<double> _H; 
    
    // Solveur direct Cholesky (rapide pour matrices symétriques définies positives)
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> _solver;

public:
    // Constructeur : Prend le pointeur de grille par valeur (Correction appliquée)
    Laplacian(Function* function, DataFile* data_file, MACgrid* grid);

    // Construit la matrice du Laplacien avec les CL de Neumann
    void BuildMatrix();

    // Calcule la divergence du champ de vitesse (div u)
    Eigen::VectorXd ComputeDivergence(const Eigen::VectorXd& U, const Eigen::VectorXd& V);

    // Résout le système linéaire H * P = RHS
    void Solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& p_sol);

    // Calcule le gradient de pression pour la correction de vitesse
    void ComputeGradient(const Eigen::VectorXd& p, Eigen::VectorXd& gradPx, Eigen::VectorXd& gradPy);
    
    // Accesseur (pour debug ou validation)
    const Eigen::SparseMatrix<double>& Get_H() const { return _H; }
};

#endif