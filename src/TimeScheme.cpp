// ====================================================================================
//                                  TIME_SCHEME.CPP
// ====================================================================================
// Implémentation des schémas temporels (Euler, RK2, RK4).
//
// AMELIORATIONS APPORTEES :
// 1. Gestion mémoire optimisée (plus d'allocations dynamiques dans les boucles).
// 2. Schéma d'advection HYBRIDE (Upwind + Centré) pour moins de diffusion.
// 3. Calculs de coordonnées optimisés (accès direct sans objets temporaires).
// ====================================================================================

#include "TimeScheme.h"
#include "DataFile.h"
#include "Laplacian.h"
#include "MACgrid.h"
#include "Function.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm> 
#include <locale>
#include <cmath>

using namespace Eigen;
using namespace std;

// =========================================================================
// HELPER : SCHÉMA D'ADVECTION HYBRIDE
// =========================================================================
// Mélange Upwind (Stable mais diffusif) et Centré (Précis mais instable).
// gamma = 0.0 -> Pur Upwind (Diffusion max, très stable)
// gamma = 1.0 -> Pur Centré (Diffusion nulle, oscillations possibles)
// gamma = 0.5 -> Bon compromis (Moins de diffusion, reste stable)
inline double AdvectionHybrid(double vel, double val_minus, double val_center, double val_plus, double inv_h) {
    // 1. Contribution Upwind (Ordre 1)
    double flux_upwind = 0.0;
    if (vel > 0) flux_upwind = vel * (val_center - val_minus) * inv_h;
    else         flux_upwind = vel * (val_plus - val_center) * inv_h;

    // 2. Contribution Centrée (Ordre 2)
    double flux_centered = vel * (val_plus - val_minus) * 0.5 * inv_h;

    // 3. Mélange (Facteur GAMMA)
    const double GAMMA = 0.5; 
    return (1.0 - GAMMA) * flux_upwind + GAMMA * flux_centered;
}

// =========================================================================
// CONSTRUCTEURS (AVEC PRE-ALLOCATION)
// =========================================================================

// Constructeur Mère
TimeScheme::TimeScheme(DataFile* df, Laplacian* lap, MACgrid* grid) 
    : _df(df), _lap(lap), _grid(grid), _t(df->Get_t0())
{
    long size_u = _grid->GetU().size();
    long size_v = _grid->GetV().size();
    long size_p = _grid->GetP().size();

    // Allocation unique des buffers de travail
    _du.resize(size_u);     _dv.resize(size_v);
    _u_star.resize(size_u); _v_star.resize(size_v);
    _u_next.resize(size_u); _v_next.resize(size_v);
    
    _p_next.resize(size_p);
    _div.resize(size_p);
    _rhs.resize(size_p);
    
    _gradPx.resize(size_u);
    _gradPy.resize(size_v);

    // Initialisation à zéro
    _du.setZero();      _dv.setZero();
}

// Constructeur Euler
EulerScheme::EulerScheme(DataFile* df, Laplacian* lap, MACgrid* grid) 
    : TimeScheme(df, lap, grid) 
{}

// Constructeur RK2
RungeKutta2Scheme::RungeKutta2Scheme(DataFile* df, Laplacian* lap, MACgrid* grid) 
    : TimeScheme(df, lap, grid) 
{
    long size_u = _grid->GetU().size();
    long size_v = _grid->GetV().size();
    _k1_u.resize(size_u); _k1_v.resize(size_v);
    _u_tmp.resize(size_u); _v_tmp.resize(size_v);
}

// Constructeur RK4
RungeKutta4Scheme::RungeKutta4Scheme(DataFile* df, Laplacian* lap, MACgrid* grid) 
    : TimeScheme(df, lap, grid) 
{
    long size_u = _grid->GetU().size();
    long size_v = _grid->GetV().size();
    _k1_u.resize(size_u); _k1_v.resize(size_v);
    _k2_u.resize(size_u); _k2_v.resize(size_v);
    _k3_u.resize(size_u); _k3_v.resize(size_v);
    _k4_u.resize(size_u); _k4_v.resize(size_v);
    _u_tmp.resize(size_u); _v_tmp.resize(size_v);
}

// =========================================================================
// OUTILS & BC
// =========================================================================

void TimeScheme::ApplyBoundaryConditions()
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    VectorXd& U = const_cast<VectorXd&>(_grid->GetU()); 
    VectorXd& V = const_cast<VectorXd&>(_grid->GetV());
    Function* fct = _grid->GetFunction();
    double ymin = _df->Get_ymin();
    double hy = _df->Get_hy();
    double xmin = _df->Get_xmin();
    double hx = _df->Get_hx();

    // Murs Verticaux
    for (int i = 0; i < Ny; ++i) {
        int k_left = _grid->GetUIndex(i, 0);
        // Optimisation coordonnée : calcul direct
        double y_left = ymin + (i + 0.5) * hy;

        if (fct->IsDirichletLeft()) U(k_left) = fct->GetLeftU_Normal(y_left);
        else U(k_left) = U(_grid->GetUIndex(i, 1)); 

        int k_right = _grid->GetUIndex(i, Nx);
        double y_right = ymin + (i + 0.5) * hy;

        if (fct->IsDirichletRight()) U(k_right) = fct->GetRightU_Normal(y_right);
        else U(k_right) = U(_grid->GetUIndex(i, Nx - 1));
    }

    // Murs Horizontaux
    for (int j = 0; j < Nx; ++j) {
        int k_bott = _grid->GetVIndex(0, j);
        double x_bott = xmin + (j + 0.5) * hx;

        if (fct->IsDirichletBottom()) V(k_bott) = fct->GetBottomV_Normal(x_bott); 
        else V(k_bott) = V(_grid->GetVIndex(1, j));

        int k_top = _grid->GetVIndex(Ny, j);
        double x_top = xmin + (j + 0.5) * hx;

        if (fct->IsDirichletTop()) V(k_top) = fct->GetTopV_Normal(x_top); 
        else V(k_top) = V(_grid->GetVIndex(Ny - 1, j));
    }
}

// =========================================================================
// COMPUTE TENDENCY (C'est ici que la Physique est calculée !)
// =========================================================================
void TimeScheme::ComputeTendency(const VectorXd& u_in, const VectorXd& v_in, VectorXd& du, VectorXd& dv)
{
    du.setZero();
    dv.setZero();

    double nu = _df->Get_nu();
    double hx = _df->Get_hx(); double hy = _df->Get_hy();
    int Nx = _df->Get_Nx(); int Ny = _df->Get_Ny();
    double odx = 1.0 / hx; double ody = 1.0 / hy;
    double odx2 = 1.0 / (hx*hx); double ody2 = 1.0 / (hy*hy);
    double xmin = _df->Get_xmin();
    double ymin = _df->Get_ymin();

    Function* fct = _grid->GetFunction();

    // --- TENDANCE U ---
    for (int i = 0; i < Ny; ++i) {
        for (int j = 1; j < Nx; ++j) { 
            if (_grid->IsSolidU(i,j)) continue; 

            int k = _grid->GetUIndex(i, j);
            double u_curr = u_in(k);
            double x_curr = xmin + j * hx; // Coordonnée U
            
            double u_E = u_in(_grid->GetUIndex(i, j + 1));
            double u_W = u_in(_grid->GetUIndex(i, j - 1));
            double u_N, u_S;
            
            if (i > 0) u_S = u_in(_grid->GetUIndex(i - 1, j));
            else { 
                if (fct->IsDirichletBottom()) u_S = 2.0 * fct->GetBottomU_Tangent(x_curr) - u_curr;
                else u_S = u_curr; 
            }
            if (i < Ny - 1) u_N = u_in(_grid->GetUIndex(i + 1, j));
            else {
                if (fct->IsDirichletTop()) u_N = 2.0 * fct->GetTopU_Tangent(x_curr) - u_curr;
                else u_N = u_curr;
            }

            double diffusion = nu * ((u_E - 2*u_curr + u_W)*odx2 + (u_N - 2*u_curr + u_S)*ody2);
            
            // Advection X (Hybride)
            double adv_x = AdvectionHybrid(u_curr, u_W, u_curr, u_E, odx);

            // Advection Y (Hybride)
            double v_avg = 0.25 * (
                (i < Ny ? v_in(_grid->GetVIndex(i + 1, j)) : 0.0) +      
                (i < Ny ? v_in(_grid->GetVIndex(i + 1, j - 1)) : 0.0) + 
                v_in(_grid->GetVIndex(i, j)) +                          
                v_in(_grid->GetVIndex(i, j - 1))                         
            );
            double adv_y = AdvectionHybrid(v_avg, u_S, u_curr, u_N, ody);

            du(k) = diffusion - (adv_x + adv_y);
        }
    }

    // --- TENDANCE V ---
    for (int i = 1; i < Ny; ++i) { 
        for (int j = 0; j < Nx; ++j) {
            if (_grid->IsSolidV(i,j)) continue; 

            int k = _grid->GetVIndex(i, j);
            double v_curr = v_in(k);
            double y_curr = ymin + i * hy; // Coordonnée V

            double v_N = v_in(_grid->GetVIndex(i + 1, j));
            double v_S = v_in(_grid->GetVIndex(i - 1, j));
            double v_E, v_W;

            if (j > 0) v_W = v_in(_grid->GetVIndex(i, j - 1));
            else {
                if (fct->IsDirichletLeft()) v_W = 2.0 * fct->GetLeftV_Tangent(y_curr) - v_curr;
                else v_W = v_curr;
            }
            if (j < Nx - 1) v_E = v_in(_grid->GetVIndex(i, j + 1));
            else {
                if (fct->IsDirichletRight()) v_E = 2.0 * fct->GetRightV_Tangent(y_curr) - v_curr;
                else v_E = v_curr;
            }

            double diffusion = nu * ((v_E - 2*v_curr + v_W)*odx2 + (v_N - 2*v_curr + v_S)*ody2);
            
            // Advection Y (Hybride)
            double adv_y = AdvectionHybrid(v_curr, v_S, v_curr, v_N, ody);
            
            // Advection X (Hybride)
            double u_avg = 0.25 * (
                (j < Nx ? u_in(_grid->GetUIndex(i, j + 1)) : 0.0) +      
                u_in(_grid->GetUIndex(i, j)) +                          
                (j < Nx ? u_in(_grid->GetUIndex(i - 1, j + 1)) : 0.0) + 
                u_in(_grid->GetUIndex(i - 1, j))                         
            );
            double adv_x = AdvectionHybrid(u_avg, v_W, v_curr, v_E, odx);

            dv(k) = diffusion - (adv_x + adv_y);
        }
    }
}

// =========================================================================
// 1. EULER EXPLICITE
// =========================================================================
void EulerScheme::Advance()
{
    ApplyBoundaryConditions(); 

    const VectorXd& u_n = _grid->GetU();
    const VectorXd& v_n = _grid->GetV();

    ComputeTendency(u_n, v_n, _du, _dv);

    double dt = _df->Get_dt();
    
    // Optimisation : Utilisation des buffers membres
    _u_star = u_n + dt * _du;
    _v_star = v_n + dt * _dv;

    double rho = _df->Get_rho();
    
    // Optimisation : Passage par référence
    _lap->ComputeDivergence(_u_star, _v_star, _div);
    
    _rhs = (rho / dt) * _div;
    
    _lap->Solve(_rhs, _p_next);
    _lap->ComputeGradient(_p_next, _gradPx, _gradPy);

    _u_next = _u_star - (dt / rho) * _gradPx;
    _v_next = _v_star - (dt / rho) * _gradPy;

    // Pénalisation
    int Nx = _df->Get_Nx(); int Ny = _df->Get_Ny();
    for (int i=0; i<Ny; ++i) for (int j=0; j<=Nx; ++j) if (_grid->IsSolidU(i,j)) _u_next(_grid->GetUIndex(i,j))=0.;
    for (int i=0; i<=Ny; ++i) for (int j=0; j<Nx; ++j) if (_grid->IsSolidV(i,j)) _v_next(_grid->GetVIndex(i,j))=0.;

    _grid->SetU(_u_next); _grid->SetV(_v_next); _grid->SetP(_p_next);
    _t += dt;
}

// =========================================================================
// 2. RUNGE-KUTTA 2
// =========================================================================
void RungeKutta2Scheme::Advance()
{
    double dt = _df->Get_dt();
    const VectorXd& u_old = _grid->GetU();
    const VectorXd& v_old = _grid->GetV();
    
    ApplyBoundaryConditions(); 

    // K1
    ComputeTendency(u_old, v_old, _k1_u, _k1_v);

    // K2 (Point milieu)
    _u_tmp = u_old + 0.5 * dt * _k1_u;
    _v_tmp = v_old + 0.5 * dt * _k1_v;
    _grid->SetU(_u_tmp); _grid->SetV(_v_tmp);
    ApplyBoundaryConditions(); 
    ComputeTendency(_grid->GetU(), _grid->GetV(), _du, _dv); // _du sert de k2

    // Prediction
    _u_star = u_old + dt * _du;
    _v_star = v_old + dt * _dv;

    // Projection
    double rho = _df->Get_rho();
    _lap->ComputeDivergence(_u_star, _v_star, _div); // Optimisation
    
    _rhs = (rho / dt) * _div;
    
    _lap->Solve(_rhs, _p_next);
    _lap->ComputeGradient(_p_next, _gradPx, _gradPy);

    _u_next = _u_star - (dt / rho) * _gradPx;
    _v_next = _v_star - (dt / rho) * _gradPy;

    // Pénalisation
    int Nx = _df->Get_Nx(); int Ny = _df->Get_Ny();
    for (int i=0; i<Ny; ++i) for (int j=0; j<=Nx; ++j) if (_grid->IsSolidU(i,j)) _u_next(_grid->GetUIndex(i,j))=0.;
    for (int i=0; i<=Ny; ++i) for (int j=0; j<Nx; ++j) if (_grid->IsSolidV(i,j)) _v_next(_grid->GetVIndex(i,j))=0.;

    _grid->SetU(_u_next); _grid->SetV(_v_next); _grid->SetP(_p_next);
    _t += dt;
}

// =========================================================================
// 3. RUNGE-KUTTA 4
// =========================================================================
void RungeKutta4Scheme::Advance()
{
    double dt = _df->Get_dt();
    const VectorXd& u_old = _grid->GetU();
    const VectorXd& v_old = _grid->GetV();

    ApplyBoundaryConditions(); 

    // K1
    ComputeTendency(u_old, v_old, _k1_u, _k1_v);

    // K2
    _u_tmp = u_old + 0.5 * dt * _k1_u;
    _v_tmp = v_old + 0.5 * dt * _k1_v;
    _grid->SetU(_u_tmp); _grid->SetV(_v_tmp); ApplyBoundaryConditions();
    ComputeTendency(_grid->GetU(), _grid->GetV(), _k2_u, _k2_v);

    // K3
    _u_tmp = u_old + 0.5 * dt * _k2_u;
    _v_tmp = v_old + 0.5 * dt * _k2_v;
    _grid->SetU(_u_tmp); _grid->SetV(_v_tmp); ApplyBoundaryConditions();
    ComputeTendency(_grid->GetU(), _grid->GetV(), _k3_u, _k3_v);

    // K4
    _u_tmp = u_old + dt * _k3_u;
    _v_tmp = v_old + dt * _k3_v;
    _grid->SetU(_u_tmp); _grid->SetV(_v_tmp); ApplyBoundaryConditions();
    ComputeTendency(_grid->GetU(), _grid->GetV(), _k4_u, _k4_v);

    // Combinaison RK4
    _u_star = u_old + (dt / 6.0) * (_k1_u + 2.0*_k2_u + 2.0*_k3_u + _k4_u);
    _v_star = v_old + (dt / 6.0) * (_k1_v + 2.0*_k2_v + 2.0*_k3_v + _k4_v);

    // Projection
    double rho = _df->Get_rho();
    _lap->ComputeDivergence(_u_star, _v_star, _div); // Optimisation
    
    _rhs = (rho / dt) * _div;
    
    _lap->Solve(_rhs, _p_next);
    _lap->ComputeGradient(_p_next, _gradPx, _gradPy);

    _u_next = _u_star - (dt / rho) * _gradPx;
    _v_next = _v_star - (dt / rho) * _gradPy;

    // Pénalisation
    int Nx = _df->Get_Nx(); int Ny = _df->Get_Ny();
    for (int i=0; i<Ny; ++i) for (int j=0; j<=Nx; ++j) if (_grid->IsSolidU(i,j)) _u_next(_grid->GetUIndex(i,j))=0.;
    for (int i=0; i<=Ny; ++i) for (int j=0; j<Nx; ++j) if (_grid->IsSolidV(i,j)) _v_next(_grid->GetVIndex(i,j))=0.;

    _grid->SetU(_u_next); _grid->SetV(_v_next); _grid->SetP(_p_next);
    _t += dt;
}

// =========================================================================
// SAUVEGARDE MANUELLE
// =========================================================================
void TimeScheme::SaveSolution(int n_iteration) {
    stringstream ss_dat;
    string resultsPath = _df->Get_results();
    ss_dat << resultsPath << "/sol_" << n_iteration << ".dat";
    ofstream dat(ss_dat.str());

    if (dat.is_open()) {
        dat.imbue(std::locale("C"));
        int Nx = _df->Get_Nx(); int Ny = _df->Get_Ny();
        double hx = _df->Get_hx(); double hy = _df->Get_hy();
        double xmin = _df->Get_xmin();
        double ymin = _df->Get_ymin();
        const VectorXd& P = _grid->GetP();
        const VectorXd& U = _grid->GetU();
        const VectorXd& V = _grid->GetV();

        auto get_u = [&](int r, int c) { return 0.5*(U(_grid->GetUIndex(r,c)) + U(_grid->GetUIndex(r,c+1))); };
        auto get_v = [&](int r, int c) { return 0.5*(V(_grid->GetVIndex(r,c)) + V(_grid->GetVIndex(r+1,c))); };

        for (int i = 0; i < Ny; ++i) {
            for (int j = 0; j < Nx; ++j) {
                double x = xmin + (j + 0.5) * hx;
                double y = ymin + (i + 0.5) * hy;
                double p_val = P(_grid->GetPIndex(i, j));
                double u_val = get_u(i, j);
                double v_val = get_v(i, j);
                
                // Calcul Vorticité
                double dv_dx = (j<Nx-1 && j>0) ? (get_v(i,j+1)-get_v(i,j-1))/(2*hx) : 0; 
                double du_dy = (i<Ny-1 && i>0) ? (get_u(i+1,j)-get_u(i-1,j))/(2*hy) : 0;
                double omega = dv_dx - du_dy;

                dat << x << " " << y << " " << p_val << " " << u_val << " " << v_val << " " << omega << "\n";
            }
            dat << "\n";
        }
    }
    dat.close();
}