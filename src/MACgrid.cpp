// ====================================================================================
//                                  MACGRID.CPP
// ====================================================================================
// Implémentation de la grille MAC et de l'initialisation des champs.
// ====================================================================================

#include "MACgrid.h"
#include "Function.h"
#include "DataFile.h"
#include <iostream>
#include <cmath>

// ============================================================================
// CONSTRUCTEUR
// ============================================================================
MACgrid::MACgrid(Function* fct, DataFile* df) : _fct(fct), _df(df) 
{
    _Nx = _df->Get_Nx();
    _Ny = _df->Get_Ny();
    
    // --- ALLOCATION MÉMOIRE ---
    _p.resize(_Nx * _Ny);
    _U.resize((_Nx + 1) * _Ny);     
    _V.resize(_Nx * (_Ny + 1));     
    
    _is_solid_u.assign((_Nx + 1) * _Ny, 0);
    _is_solid_v.assign(_Nx * (_Ny + 1), 0);
    
    // 1. Définition de la géométrie (Obstacles)
    this->BuildObstacles(); 
    
    // 2. Initialisation des champs
    _p.setZero();

    // --- INITIALISATION U ---
    for (int i = 0; i < _Ny; ++i) {
        for (int j = 0; j <= _Nx; ++j) {
            int k = GetUIndex(i, j);
            if (_is_solid_u[k]) {
                _U(k) = 0.0;
            } else {
                Eigen::VectorXd c = GetUcoord(i, j);
                _U(k) = _fct->InitialConditionU(c(0), c(1));
            }
        }
    }

    // --- INITIALISATION V ---
    for (int i = 0; i <= _Ny; ++i) {
        for (int j = 0; j < _Nx; ++j) {
            int k = GetVIndex(i, j);
            if (_is_solid_v[k]) {
                _V(k) = 0.0;
            } else {
                Eigen::VectorXd c = GetVcoord(i, j);
                _V(k) = _fct->InitialConditionV(c(0), c(1));
            }
        }
    }
}

// ============================================================================
// COORDONNÉES SPATIALES
// ============================================================================
// Note : P est centré en (i+0.5, j+0.5)
//        U est décalé en (i+0.5, j)
//        V est décalé en (i, j+0.5)

Eigen::VectorXd MACgrid::GetPcoord(int i, int j) const {
    return Eigen::Vector2d(_df->Get_xmin() + (j + 0.5) * _df->Get_hx(),
                           _df->Get_ymin() + (i + 0.5) * _df->Get_hy());
}

Eigen::VectorXd MACgrid::GetUcoord(int i, int j) const {
    return Eigen::Vector2d(_df->Get_xmin() + j * _df->Get_hx(),
                           _df->Get_ymin() + (i + 0.5) * _df->Get_hy());
}

Eigen::VectorXd MACgrid::GetVcoord(int i, int j) const {
    return Eigen::Vector2d(_df->Get_xmin() + (j + 0.5) * _df->Get_hx(),
                           _df->Get_ymin() + i * _df->Get_hy());
}

// ============================================================================
// BUILD OBSTACLES
// ============================================================================
void MACgrid::BuildObstacles() {
    double cx = _df->Get_CylCx(); 
    double cy = _df->Get_CylCy(); 
    double R  = _df->Get_CylRadius(); 
    double R2 = R * R;

    // Si rayon nul ou négatif, pas d'obstacle
    if (R <= 0.0) return;

    // --- MARQUAGE OBSTACLE SUR U ---
    for (int i = 0; i < _Ny; ++i) {
        for (int j = 0; j <= _Nx; ++j) {
            Eigen::VectorXd c = GetUcoord(i, j);
            // Distance au carré
            if (std::pow(c(0)-cx, 2) + std::pow(c(1)-cy, 2) <= R2) {
                _is_solid_u[GetUIndex(i, j)] = 1;
                _U(GetUIndex(i, j)) = 0.0; // Vitesse nulle dans le solide
            }
        }
    }

    // --- MARQUAGE OBSTACLE SUR V ---
    for (int i = 0; i <= _Ny; ++i) {
        for (int j = 0; j < _Nx; ++j) {
            Eigen::VectorXd c = GetVcoord(i, j);
            if (std::pow(c(0)-cx, 2) + std::pow(c(1)-cy, 2) <= R2) {
                _is_solid_v[GetVIndex(i, j)] = 1;
                _V(GetVIndex(i, j)) = 0.0;
            }
        }
    }
    
    std::cout << ">> [MACgrid] Obstacle construit : Cylindre R=" << R 
              << " en (" << cx << "," << cy << ")" << std::endl;
}