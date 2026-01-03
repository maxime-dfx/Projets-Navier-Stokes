#ifndef _TIME_SCHEME_H_
#define _TIME_SCHEME_H_

#include "DataFile.h"
#include "Laplacian.h"
#include "MACgrid.h"
#include <string>
#include <vector>

class TimeScheme
{
protected:
    DataFile* _df;
    Laplacian* _lap;
    MACgrid* _grid;
    double _t;
    
    // [OPTIMISATION] Buffers pré-alloués pour les calculs intermédiaires
    // Utilisés par ComputeTendency et les schémas
    Eigen::VectorXd _du; 
    Eigen::VectorXd _dv;

    // Applique les conditions Dirichlet/Neumann aux frontières sur la grille actuelle
    void ApplyBoundaryConditions();

    // Calcule les variations (Diffusion - Advection)
    void ComputeTendency(const Eigen::VectorXd& u_in, const Eigen::VectorXd& v_in, 
                         Eigen::VectorXd& du, Eigen::VectorXd& dv);

public:
    TimeScheme(DataFile* data_file, Laplacian* lap, MACgrid* grid);
    virtual ~TimeScheme();
    
    virtual void Advance() = 0;
    
    double GetTime() const { return _t; }
    void SaveSolution(int n_iteration);
};

// Euler Explicite (Ordre 1)
class EulerScheme : public TimeScheme
{
public:
    EulerScheme(DataFile* data_file, Laplacian* lap, MACgrid* grid);
    void Advance() override;
};

// Runge-Kutta 2 (Ordre 2 - Point milieu)
class RungeKutta2Scheme : public TimeScheme
{
private:
    // Buffers spécifiques RK2
    Eigen::VectorXd _k1_u, _k1_v;
    Eigen::VectorXd _u_tmp, _v_tmp; 

public:
    RungeKutta2Scheme(DataFile* data_file, Laplacian* lap, MACgrid* grid);
    void Advance() override;
};

// Runge-Kutta 4 (Ordre 4 - Standard)
class RungeKutta4Scheme : public TimeScheme
{
private:
    // Buffers spécifiques RK4
    Eigen::VectorXd _k1_u, _k1_v;
    Eigen::VectorXd _k2_u, _k2_v;
    Eigen::VectorXd _k3_u, _k3_v;
    Eigen::VectorXd _k4_u, _k4_v;
    Eigen::VectorXd _u_tmp, _v_tmp;

public:
    RungeKutta4Scheme(DataFile* data_file, Laplacian* lap, MACgrid* grid);
    void Advance() override;
};

#endif