// ====================================================================================
//                                  TIME_SCHEME.H
// ====================================================================================
// Description : Classe abstraite et dérivées pour l'intégration temporelle.
//               Gère Euler, RK2 et RK4.
//               Contient la physique (Advection + Diffusion) dans ComputeTendency.
// ====================================================================================

#ifndef _TIME_SCHEME_H_
#define _TIME_SCHEME_H_

#include <Eigen/Dense>
#include <vector>

// Forward declarations
class DataFile;
class Laplacian;
class MACgrid;

// Classe de base abstraite
class TimeScheme {
protected:
    DataFile* _df;
    Laplacian* _lap;
    MACgrid* _grid;
    double _t;

    // Buffers pour stocker les tendances (du/dt, dv/dt)
    Eigen::VectorXd _du;
    Eigen::VectorXd _dv;

    // --- Méthodes Internes (Protected) ---
    void ApplyBoundaryConditions();
    void ComputeTendency(const Eigen::VectorXd& u_in, const Eigen::VectorXd& v_in, 
                         Eigen::VectorXd& du_out, Eigen::VectorXd& dv_out);

public:
    // Constructeur de base
    TimeScheme(DataFile* df, Laplacian* lap, MACgrid* grid);
    virtual ~TimeScheme() = default;

    // Fait avancer la simulation de t à t + dt
    virtual void Advance() = 0;

    // Sauvegarde manuelle (.dat)
    void SaveSolution(int n_iteration);

    double GetTime() const { return _t; }
};

// ============================================================================
// EULER EXPLICITE
// ============================================================================
class EulerScheme : public TimeScheme {
public:
    // Déclaration explicite du constructeur pour matcher le .cpp
    EulerScheme(DataFile* df, Laplacian* lap, MACgrid* grid);
    void Advance() override;
};

// ============================================================================
// RUNGE-KUTTA 2
// ============================================================================
class RungeKutta2Scheme : public TimeScheme {
private:
    Eigen::VectorXd _k1_u, _k1_v;
    Eigen::VectorXd _u_tmp, _v_tmp;
public:
    RungeKutta2Scheme(DataFile* df, Laplacian* lap, MACgrid* grid);
    void Advance() override;
};

// ============================================================================
// RUNGE-KUTTA 4
// ============================================================================
class RungeKutta4Scheme : public TimeScheme {
private:
    Eigen::VectorXd _k1_u, _k1_v;
    Eigen::VectorXd _k2_u, _k2_v;
    Eigen::VectorXd _k3_u, _k3_v;
    Eigen::VectorXd _k4_u, _k4_v;
    Eigen::VectorXd _u_tmp, _v_tmp;
public:
    RungeKutta4Scheme(DataFile* df, Laplacian* lap, MACgrid* grid);
    void Advance() override;
};

#endif // _TIME_SCHEME_H_

