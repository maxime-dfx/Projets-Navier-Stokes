// ====================================================================================
//                                  DATAFILE.H
// ====================================================================================
// Description : Classe de configuration (Singleton de fait).
//               Lit le fichier .toml et expose les paramètres de simulation
//               en lecture seule pour tous les autres modules.
// ====================================================================================

#ifndef _DATA_FILE_H
#define _DATA_FILE_H

#include <string>
#include <vector>
#include <iostream>

class DataFile {
private:
    // --- MAILLAGE (Space) ---
    double _xmin, _xmax;
    double _ymin, _ymax;
    double _hx, _hy;
    int _Nx, _Ny;

    // --- TEMPS & SCHÉMA (Time) ---
    double _t0, _tfinal, _dt;
    std::string _scheme; // "Euler", "RK2", "RK4"

    // --- PHYSIQUE (Physics) ---
    double _rho;
    double _nu; // Viscosité cinématique

    // --- CONDITIONS AUX LIMITES (Boundary) ---
    std::string _bc_left, _bc_right, _bc_bottom, _bc_top;
    double _bc_left_val, _bc_right_val, _bc_bottom_val, _bc_top_val;
    std::string _inlet_profile; // "Constant" ou "Poiseuille"

    // --- OBSTACLE (Cylindre) ---
    double _cyl_cx, _cyl_cy, _cyl_r;

    // --- SORTIES (Output) ---
    std::string _results_dir;
    std::string _sim_name;
    const std::string _file_name;

    // --- PARAMÈTRES ÉTUDE CONVERGENCE ---
    bool _has_convergence_section;
    std::vector<int> _conv_resolutions;
    std::vector<std::string> _conv_schemes;
    double _conv_L, _conv_CFL, _conv_Nu;

public:
    // Constructeur : Parse le fichier TOML
    explicit DataFile(std::string file_name);

    // ========================================================================
    // GETTERS (Accesseurs Lecture Seule)
    // ========================================================================
    
    // Espace
    double Get_xmin() const { return _xmin; }
    double Get_xmax() const { return _xmax; }
    double Get_ymin() const { return _ymin; }
    double Get_ymax() const { return _ymax; }
    double Get_hx()   const { return _hx; }
    double Get_hy()   const { return _hy; }
    int    Get_Nx()   const { return _Nx; }
    int    Get_Ny()   const { return _Ny; }

    // Temps & Physique
    double Get_t0()     const { return _t0; }
    double Get_tfinal() const { return _tfinal; }
    double Get_dt()     const { return _dt; }
    double Get_rho()    const { return _rho; }
    double Get_nu()     const { return _nu; }
    std::string Get_scheme() const { return _scheme; }

    // Obstacle
    double Get_CylCx()     const { return _cyl_cx; }
    double Get_CylCy()     const { return _cyl_cy; }
    double Get_CylRadius() const { return _cyl_r; }

    // Conditions Limites (BC)
    // Retourne le type ("Dirichlet", "Neumann") ou la valeur
    std::string Get_BC_Left()   const { return _bc_left; }
    std::string Get_BC_Right()  const { return _bc_right; }
    std::string Get_BC_Bottom() const { return _bc_bottom; }
    std::string Get_BC_Top()    const { return _bc_top; }

    double Get_BC_Left_dir()   const { return _bc_left_val; }
    double Get_BC_Right_dir()  const { return _bc_right_val; }
    double Get_BC_Bottom_dir() const { return _bc_bottom_val; }
    double Get_BC_Top_dir()    const { return _bc_top_val; }
    
    std::string Get_InletProfile() const { return _inlet_profile; }

    // Output
    std::string Get_results() const { return _results_dir; }
    std::string Get_SimName() const { return _sim_name; }

    // Convergence (Pour le main)
    bool HasConvergenceSection() const { return _has_convergence_section; }
    const std::vector<int>& Get_Conv_Resolutions() const { return _conv_resolutions; }
    const std::vector<std::string>& Get_Conv_Schemes() const { return _conv_schemes; }
    double Get_Conv_L()   const { return _conv_L; }
    double Get_Conv_CFL() const { return _conv_CFL; }
    double Get_Conv_Nu()  const { return _conv_Nu; }

    // ========================================================================
    // SETTERS (Uniquement pour l'étude de convergence dynamique)
    // ========================================================================
    void Set_Nx(int nx) { _Nx = nx; }
    void Set_Ny(int ny) { _Ny = ny; }
    void Set_hx(double h) { _hx = h; }
    void Set_hy(double h) { _hy = h; }
    void Set_dt(double dt) { _dt = dt; }
    void Set_Scheme(std::string s) { _scheme = s; }
};

#endif // _DATA_FILE_H