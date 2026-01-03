#include "MACgrid.h"
#include <iostream>
#include <cassert> 

using namespace std;
using namespace Eigen;

MACgrid::MACgrid(Function* function, DataFile* data_file) :
_fct(function), _df(data_file)
{
    _Nx = _df->Get_Nx();
    _Ny = _df->Get_Ny();
    
    _p.resize(_Nx * _Ny);
    _U.resize((_Nx + 1) * _Ny);     
    _V.resize(_Nx * (_Ny + 1));     
    
    // [OPTIMISATION] 0 vaut 'false', 1 vaut 'true'
    _is_solid_u.resize((_Nx + 1) * _Ny, 0);
    _is_solid_v.resize(_Nx * (_Ny + 1), 0);
    
    // 1. Construire les obstacles AVANT d'initialiser U/V
    this->BuildObstacles(); 
    
    _p.setZero();

    // 2. Initialisation U
    for (int i = 0; i < _Ny; ++i) {
        for (int j = 0; j <= _Nx; ++j) {
            int k = GetUIndex(i, j);
            VectorXd coord = GetUcoord(i, j);
            
            // Si on est dans le solide, U=0, sinon condition initiale
            if (_is_solid_u[k]) _U(k) = 0.0;
            else _U(k) = _fct->InitialConditionU(coord(0), coord(1));
        }
    }

    // 3. Initialisation V
    for (int i = 0; i <= _Ny; ++i) {
        for (int j = 0; j < _Nx; ++j) {
            int k = GetVIndex(i, j);
            VectorXd coord = GetVcoord(i, j);
            
            if (_is_solid_v[k]) _V(k) = 0.0;
            else _V(k) = _fct->InitialConditionV(coord(0), coord(1));
        }
    }
}

int MACgrid::GetPIndex(int i, int j) const { return j + i * _Nx; }
int MACgrid::GetUIndex(int i, int j) const { return j + i * (_Nx + 1); }
int MACgrid::GetVIndex(int i, int j) const { return j + i * _Nx; }

VectorXd MACgrid::GetPcoord(int i, int j) const {
    VectorXd coord(2);
    coord(0) = _df->Get_xmin() + (j + 0.5) * _df->Get_hx(); 
    coord(1) = _df->Get_ymin() + (i + 0.5) * _df->Get_hy();
    return coord;
}
VectorXd MACgrid::GetUcoord(int i, int j) const {
    VectorXd coord(2);
    coord(0) = _df->Get_xmin() + j * _df->Get_hx(); 
    coord(1) = _df->Get_ymin() + (i + 0.5) * _df->Get_hy();
    return coord;
}
VectorXd MACgrid::GetVcoord(int i, int j) const {
    VectorXd coord(2);
    coord(0) = _df->Get_xmin() + (j + 0.5) * _df->Get_hx();
    coord(1) = _df->Get_ymin() + i * _df->Get_hy();
    return coord;
}

void MACgrid::BuildObstacles() {
    double x_min = _df->Get_xmin();
    double y_min = _df->Get_ymin();
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();

    double cx = _df->Get_CylCx(); 
    double cy = _df->Get_CylCy(); 
    double R  = _df->Get_CylRadius(); 
    double R2 = R * R;

    // Masque U
    for (int i = 0; i < _Ny; ++i) {
        for (int j = 0; j <= _Nx; ++j) {
            double x = x_min + j * hx;
            double y = y_min + (i + 0.5) * hy;
            if ( (x-cx)*(x-cx) + (y-cy)*(y-cy) <= R2 ) {
                _is_solid_u[GetUIndex(i, j)] = 1; // 1 pour true
                _U(GetUIndex(i, j)) = 0.0;
            }
        }
    }

    // Masque V
    for (int i = 0; i <= _Ny; ++i) {
        for (int j = 0; j < _Nx; ++j) {
            double x = x_min + (j + 0.5) * hx;
            double y = y_min + i * hy; 
            if ( (x-cx)*(x-cx) + (y-cy)*(y-cy) <= R2 ) {
                _is_solid_v[GetVIndex(i, j)] = 1; // 1 pour true
                _V(GetVIndex(i, j)) = 0.0;
            }
        }
    }
    
    if (R > 0.0)
        cout << "Obstacle construit : Cylindre en (" << cx << "," << cy << ") R=" << R << endl;
    else
        cout << "Pas d'obstacle defini ou Rayon nul." << endl;
}