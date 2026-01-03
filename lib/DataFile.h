#ifndef _DATA_FILE_H
#define _DATA_FILE_H

#include <string>
#include <iostream>
#include <vector>

class DataFile {
private:
   // --- MAILLAGE ---
   double _xmin, _xmax, _ymin, _ymax, _hx, _hy;
   int _Nx, _Ny;

   // --- PHYSIQUE & TEMPS ---
   double _nu, _rho;
   double _t0, _tfinal, _dt;
   std::string _scheme;

   // --- SORTIES ---
   std::string _results;
   std::string _sim_name; // Nom de la simulation (ex: "Stokes")
   const std::string _file_name;

   // --- CONDITIONS AUX LIMITES ---
   std::string _bc_left;
   std::string _bc_right;
   std::string _bc_bottom;
   std::string _bc_top;

   double _bc_left_dir;
   double _bc_right_dir;
   double _bc_bottom_dir;
   double _bc_top_dir;

   // Type de profil d'entrée ("Constant" ou "Poiseuille")
   std::string _inlet_profile; 

   // --- OBSTACLE (CYLINDRE) ---
   double _cyl_cx;
   double _cyl_cy;
   double _cyl_r;

   // --- PARAMETRES ETUDE CONVERGENCE ---
   bool _has_convergence;       
   std::vector<int> _conv_res;  
   double _conv_L;              
   double _conv_CFL;            
   double _conv_Nu;
   std::vector<std::string> _conv_schemes;
             

public:
   DataFile(std::string file_name);

   // Getters Maillage
   double Get_xmin() const { return _xmin; }
   double Get_xmax() const { return _xmax; }
   double Get_ymin() const { return _ymin; }
   double Get_ymax() const { return _ymax; }
   double Get_hx() const { return _hx; }
   double Get_hy() const { return _hy; }
   int Get_Nx() const { return _Nx; }
   int Get_Ny() const { return _Ny; }

   // Getters Physique
   double Get_nu() const { return _nu; }
   double Get_rho() const { return _rho; }
   double Get_t0() const { return _t0; }
   double Get_tfinal() const { return _tfinal; }
   double Get_dt() const { return _dt; }
   std::string Get_scheme() const { return _scheme; }

   // Getters Sortie
   std::string Get_results() const { return _results; }
   std::string Get_SimName() const { return _sim_name; }

   // Getters BC
   std::string Get_BC_Left()   const { return _bc_left; }
   std::string Get_BC_Right()  const { return _bc_right; }
   std::string Get_BC_Bottom() const { return _bc_bottom; }
   std::string Get_BC_Top()    const { return _bc_top; }
   
   double Get_BC_Left_dir() const { return _bc_left_dir; }
   double Get_BC_Right_dir() const { return _bc_right_dir; }
   double Get_BC_Bottom_dir() const { return _bc_bottom_dir; }
   double Get_BC_Top_dir() const { return _bc_top_dir; }

   std::string Get_InletProfile() const { return _inlet_profile; }

   // Getters Obstacle
   double Get_CylCx() const { return _cyl_cx; }
   double Get_CylCy() const { return _cyl_cy; }
   double Get_CylRadius() const { return _cyl_r; }

   // --- SETTERS POUR L'ETUDE DE CONVERGENCE ---
   void Set_Nx(int nx) { _Nx = nx; }
   void Set_Ny(int ny) { _Ny = ny; }
   void Set_hx(double hx) { _hx = hx; }
   void Set_hy(double hy) { _hy = hy; }
   void Set_dt(double dt) { _dt = dt; }
   void Set_SimName(std::string name) { _sim_name = name; }
   const std::vector<std::string>& Get_Conv_Schemes() const { return _conv_schemes; }
   void Set_Scheme(std::string s) { _scheme = s; }

   // --- GETTERS POUR L'ETUDE DE CONVERGENCE ---
   bool HasConvergenceSection() const { return _has_convergence; }
   const std::vector<int>& Get_Conv_Resolutions() const { return _conv_res; }
   double Get_Conv_L() const { return _conv_L; }
   double Get_Conv_CFL() const { return _conv_CFL; }
   double Get_Conv_Nu() const { return _conv_Nu; }
};

#endif