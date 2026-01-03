#include "TimeScheme.h"
#include "Function.h" 
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm> 
#include <locale>

using namespace Eigen;
using namespace std;

// =========================================================================
// CONSTRUCTEURS
// =========================================================================

TimeScheme::TimeScheme(DataFile* data_file, Laplacian* lap, MACgrid* grid) :
    _df(data_file), _lap(lap), _grid(grid), _t(data_file->Get_t0())
{
    long size_u = _grid->GetU().size();
    long size_v = _grid->GetV().size();
    _du.resize(size_u);
    _dv.resize(size_v);
}
TimeScheme::~TimeScheme() {}

EulerScheme::EulerScheme(DataFile* data_file, Laplacian* lap, MACgrid* grid) :
    TimeScheme(data_file, lap, grid) {}

RungeKutta2Scheme::RungeKutta2Scheme(DataFile* data_file, Laplacian* lap, MACgrid* grid) :
    TimeScheme(data_file, lap, grid) 
{
    long size_u = _grid->GetU().size();
    long size_v = _grid->GetV().size();
    _k1_u.resize(size_u); _k1_v.resize(size_v);
    _u_tmp.resize(size_u); _v_tmp.resize(size_v);
}

RungeKutta4Scheme::RungeKutta4Scheme(DataFile* data_file, Laplacian* lap, MACgrid* grid) :
    TimeScheme(data_file, lap, grid) 
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

inline double Upwind(double vel, double val_minus, double val_center, double val_plus, double inv_h) {
    if (vel > 0) return vel * (val_center - val_minus) * inv_h;
    else return vel * (val_plus - val_center) * inv_h;
}

void TimeScheme::ApplyBoundaryConditions()
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    VectorXd& U = const_cast<VectorXd&>(_grid->GetU()); 
    VectorXd& V = const_cast<VectorXd&>(_grid->GetV());
    Function* fct = _grid->GetFunction();

    // Murs Verticaux
    for (int i = 0; i < Ny; ++i) {
        int k_left = _grid->GetUIndex(i, 0);
        double y_left = _grid->GetUcoord(i, 0)(1);
        if (fct->IsDirichletLeft()) U(k_left) = fct->GetLeftU_Normal(y_left);
        else U(k_left) = U(_grid->GetUIndex(i, 1)); 

        int k_right = _grid->GetUIndex(i, Nx);
        double y_right = _grid->GetUcoord(i, Nx)(1);
        if (fct->IsDirichletRight()) U(k_right) = fct->GetRightU_Normal(y_right);
        else U(k_right) = U(_grid->GetUIndex(i, Nx - 1));
    }

    // Murs Horizontaux
    for (int j = 0; j < Nx; ++j) {
        int k_bott = _grid->GetVIndex(0, j);
        double x_bott = _grid->GetVcoord(0, j)(0);
        if (fct->IsDirichletBottom()) V(k_bott) = fct->GetBottomV_Normal(x_bott); 
        else V(k_bott) = V(_grid->GetVIndex(1, j));

        int k_top = _grid->GetVIndex(Ny, j);
        double x_top = _grid->GetVcoord(Ny, j)(0);
        if (fct->IsDirichletTop()) V(k_top) = fct->GetTopV_Normal(x_top); 
        else V(k_top) = V(_grid->GetVIndex(Ny - 1, j));
    }
}

// =========================================================================
// COMPUTE TENDENCY (Inchangé)
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

    Function* fct = _grid->GetFunction();

    // --- TENDANCE U ---
    for (int i = 0; i < Ny; ++i) {
        for (int j = 1; j < Nx; ++j) { 
            if (_grid->IsSolidU(i,j)) continue; 

            int k = _grid->GetUIndex(i, j);
            double u_curr = u_in(k);
            double x_curr = _grid->GetUcoord(i,j)(0);
            
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
            double adv_x = Upwind(u_curr, u_W, u_curr, u_E, odx);

            double v_avg = 0.25 * (
                (i < Ny ? v_in(_grid->GetVIndex(i + 1, j)) : 0.0) +     
                (i < Ny ? v_in(_grid->GetVIndex(i + 1, j - 1)) : 0.0) + 
                v_in(_grid->GetVIndex(i, j)) +                          
                v_in(_grid->GetVIndex(i, j - 1))                        
            );
            double adv_y = Upwind(v_avg, u_S, u_curr, u_N, ody);

            du(k) = diffusion - (adv_x + adv_y);
        }
    }

    // --- TENDANCE V ---
    for (int i = 1; i < Ny; ++i) { 
        for (int j = 0; j < Nx; ++j) {
            if (_grid->IsSolidV(i,j)) continue; 

            int k = _grid->GetVIndex(i, j);
            double v_curr = v_in(k);
            double y_curr = _grid->GetVcoord(i, j)(1);

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
            double adv_y = Upwind(v_curr, v_S, v_curr, v_N, ody);
            
            double u_avg = 0.25 * (
                (j < Nx ? u_in(_grid->GetUIndex(i, j + 1)) : 0.0) +     
                u_in(_grid->GetUIndex(i, j)) +                          
                (j < Nx ? u_in(_grid->GetUIndex(i - 1, j + 1)) : 0.0) + 
                u_in(_grid->GetUIndex(i - 1, j))                        
            );
            double adv_x = Upwind(u_avg, v_W, v_curr, v_E, odx);

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

    // Pour Euler, pas besoin de copie car on ne fait qu'une étape
    const VectorXd& u_n = _grid->GetU();
    const VectorXd& v_n = _grid->GetV();

    ComputeTendency(u_n, v_n, _du, _dv);

    double dt = _df->Get_dt();
    VectorXd u_star = u_n + dt * _du;
    VectorXd v_star = v_n + dt * _dv;

    // --- PROJECTION ---
    double rho = _df->Get_rho();
    VectorXd div = _lap->ComputeDivergence(u_star, v_star);
    VectorXd rhs = (rho / dt) * div;
    VectorXd p_next;
    _lap->Solve(rhs, p_next);
    VectorXd gradPx, gradPy;
    _lap->ComputeGradient(p_next, gradPx, gradPy);

    VectorXd u_next = u_star - (dt / rho) * gradPx;
    VectorXd v_next = v_star - (dt / rho) * gradPy;

    // Pénalisation
    int Nx = _df->Get_Nx(); int Ny = _df->Get_Ny();
    for (int i=0; i<Ny; ++i) for (int j=0; j<=Nx; ++j) if (_grid->IsSolidU(i,j)) u_next(_grid->GetUIndex(i,j))=0.;
    for (int i=0; i<=Ny; ++i) for (int j=0; j<Nx; ++j) if (_grid->IsSolidV(i,j)) v_next(_grid->GetVIndex(i,j))=0.;

    _grid->SetU(u_next); _grid->SetV(v_next); _grid->SetP(p_next);
    _t += dt;
}

// =========================================================================
// 2. RUNGE-KUTTA 2 (CORRIGÉ)
// =========================================================================
void RungeKutta2Scheme::Advance()
{
    double dt = _df->Get_dt();
    
    // 1. SAUVEGARDE DE L'ETAT INITIAL (CRUCIAL !)
    // On copie les vecteurs car _grid->SetU va écraser la mémoire de la grille
    VectorXd u_old = _grid->GetU();
    VectorXd v_old = _grid->GetV();
    
    ApplyBoundaryConditions(); // Applique BC sur u_old/v_old implicitement car grille pas encore modifiée

    // --- K1 ---
    // Calcul de k1 à partir de l'état initial
    ComputeTendency(u_old, v_old, _k1_u, _k1_v);

    // --- K2 (Point milieu) ---
    // u_tmp = u_old + 0.5 * dt * k1
    // ATTENTION : On utilise bien u_old ici !
    _u_tmp = u_old + 0.5 * dt * _k1_u;
    _v_tmp = v_old + 0.5 * dt * _k1_v;
    
    // Mise à jour grille pour appliquer BC sur l'état intermédiaire
    _grid->SetU(_u_tmp); 
    _grid->SetV(_v_tmp);
    ApplyBoundaryConditions(); 
    
    // Calcul de k2 à partir de l'état intermédiaire
    // _du sert de buffer pour k2
    ComputeTendency(_grid->GetU(), _grid->GetV(), _du, _dv);

    // --- PREDICTION FINALE ---
    // u* = u_old + dt * k2
    // ATTENTION : On réutilise u_old (l'état initial), pas la grille modifiée !
    VectorXd u_star = u_old + dt * _du;
    VectorXd v_star = v_old + dt * _dv;

    // --- PROJECTION ---
    double rho = _df->Get_rho();
    VectorXd div = _lap->ComputeDivergence(u_star, v_star);
    VectorXd rhs = (rho / dt) * div;
    VectorXd p_next;
    _lap->Solve(rhs, p_next);
    VectorXd gradPx, gradPy;
    _lap->ComputeGradient(p_next, gradPx, gradPy);

    VectorXd u_next = u_star - (dt / rho) * gradPx;
    VectorXd v_next = v_star - (dt / rho) * gradPy;

    // Pénalisation
    int Nx = _df->Get_Nx(); int Ny = _df->Get_Ny();
    for (int i=0; i<Ny; ++i) for (int j=0; j<=Nx; ++j) if (_grid->IsSolidU(i,j)) u_next(_grid->GetUIndex(i,j))=0.;
    for (int i=0; i<=Ny; ++i) for (int j=0; j<Nx; ++j) if (_grid->IsSolidV(i,j)) v_next(_grid->GetVIndex(i,j))=0.;

    _grid->SetU(u_next); _grid->SetV(v_next); _grid->SetP(p_next);
    _t += dt;
}

// =========================================================================
// 3. RUNGE-KUTTA 4 (CORRIGÉ)
// =========================================================================
void RungeKutta4Scheme::Advance()
{
    double dt = _df->Get_dt();
    
    // 1. SAUVEGARDE DE L'ETAT INITIAL
    VectorXd u_old = _grid->GetU();
    VectorXd v_old = _grid->GetV();

    ApplyBoundaryConditions(); 

    // --- K1 ---
    ComputeTendency(u_old, v_old, _k1_u, _k1_v);

    // --- K2 ---
    // u_tmp = u_old + 0.5 * dt * k1
    _u_tmp = u_old + 0.5 * dt * _k1_u;
    _v_tmp = v_old + 0.5 * dt * _k1_v;
    _grid->SetU(_u_tmp); _grid->SetV(_v_tmp); ApplyBoundaryConditions();
    
    ComputeTendency(_grid->GetU(), _grid->GetV(), _k2_u, _k2_v);

    // --- K3 ---
    // u_tmp = u_old + 0.5 * dt * k2
    _u_tmp = u_old + 0.5 * dt * _k2_u;
    _v_tmp = v_old + 0.5 * dt * _k2_v;
    _grid->SetU(_u_tmp); _grid->SetV(_v_tmp); ApplyBoundaryConditions();

    ComputeTendency(_grid->GetU(), _grid->GetV(), _k3_u, _k3_v);

    // --- K4 ---
    // u_tmp = u_old + dt * k3
    _u_tmp = u_old + dt * _k3_u;
    _v_tmp = v_old + dt * _k3_v;
    _grid->SetU(_u_tmp); _grid->SetV(_v_tmp); ApplyBoundaryConditions();

    ComputeTendency(_grid->GetU(), _grid->GetV(), _k4_u, _k4_v);

    // --- COMBINAISON RK4 ---
    // u* = u_old + dt/6 * (k1 + 2k2 + 2k3 + k4)
    // Ici aussi, on part de u_old !
    VectorXd u_star = u_old + (dt / 6.0) * (_k1_u + 2.0*_k2_u + 2.0*_k3_u + _k4_u);
    VectorXd v_star = v_old + (dt / 6.0) * (_k1_v + 2.0*_k2_v + 2.0*_k3_v + _k4_v);

    // --- PROJECTION ---
    double rho = _df->Get_rho();
    VectorXd div = _lap->ComputeDivergence(u_star, v_star);
    VectorXd rhs = (rho / dt) * div;
    VectorXd p_next;
    _lap->Solve(rhs, p_next);
    VectorXd gradPx, gradPy;
    _lap->ComputeGradient(p_next, gradPx, gradPy);

    VectorXd u_next = u_star - (dt / rho) * gradPx;
    VectorXd v_next = v_star - (dt / rho) * gradPy;

    // Pénalisation
    int Nx = _df->Get_Nx(); int Ny = _df->Get_Ny();
    for (int i=0; i<Ny; ++i) for (int j=0; j<=Nx; ++j) if (_grid->IsSolidU(i,j)) u_next(_grid->GetUIndex(i,j))=0.;
    for (int i=0; i<=Ny; ++i) for (int j=0; j<Nx; ++j) if (_grid->IsSolidV(i,j)) v_next(_grid->GetVIndex(i,j))=0.;

    _grid->SetU(u_next); _grid->SetV(v_next); _grid->SetP(p_next);
    _t += dt;
}

void TimeScheme::SaveSolution(int n_iteration) {
    // ... (Reste inchangé)
    stringstream ss_dat;
    string resultsPath = _df->Get_results();
    ss_dat << resultsPath << "/sol_" << n_iteration << ".dat";
    ofstream dat(ss_dat.str());

    if (dat.is_open()) {
        dat.imbue(std::locale("C"));
        int Nx = _df->Get_Nx(); int Ny = _df->Get_Ny();
        double hx = _df->Get_hx(); double hy = _df->Get_hy();
        const VectorXd& P = _grid->GetP();
        const VectorXd& U = _grid->GetU();
        const VectorXd& V = _grid->GetV();

        auto get_u = [&](int r, int c) { return 0.5*(U(_grid->GetUIndex(r,c)) + U(_grid->GetUIndex(r,c+1))); };
        auto get_v = [&](int r, int c) { return 0.5*(V(_grid->GetVIndex(r,c)) + V(_grid->GetVIndex(r+1,c))); };

        for (int i = 0; i < Ny; ++i) {
            for (int j = 0; j < Nx; ++j) {
                double x = _df->Get_xmin() + (j + 0.5) * hx;
                double y = _df->Get_ymin() + (i + 0.5) * hy;
                double p_val = P(_grid->GetPIndex(i, j));
                double u_val = get_u(i, j);
                double v_val = get_v(i, j);
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