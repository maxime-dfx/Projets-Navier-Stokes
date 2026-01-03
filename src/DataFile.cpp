#include "DataFile.h"
#include <fstream>
#include <iostream>
#include "toml.hpp" // Assure-toi que ce fichier est accessible

using namespace std;

DataFile::DataFile(std::string file_name) : _file_name(file_name), 
    // Initialisation par défaut pour éviter les valeurs aléatoires
    _has_convergence(false), _conv_L(1.0), _conv_CFL(0.5), _conv_Nu(1.0)
{
    // Lecture du fichier TOML
    auto config = toml::parse(file_name);

    // --- 1. Lecture des blocs standards (Space, Time, Physics...) ---
    // (Je reprends tes blocs existants, assure-toi qu'ils sont là)
    const auto& space = toml::find(config, "space");
    _xmin = toml::find<double>(space, "xmin");
    _xmax = toml::find<double>(space, "xmax");
    _ymin = toml::find<double>(space, "ymin");
    _ymax = toml::find<double>(space, "ymax");
    _hx   = toml::find<double>(space, "hx");
    _hy   = toml::find<double>(space, "hy");
    _Nx   = (_xmax - _xmin) / _hx;
    _Ny   = (_ymax - _ymin) / _hy;

    const auto& time = toml::find(config, "time");
    _t0      = toml::find<double>(time, "t0");
    _tfinal  = toml::find<double>(time, "tfinal");
    _dt      = toml::find<double>(time, "dt");
    _scheme  = toml::find_or<std::string>(time, "scheme", "Euler");

    const auto& physics = toml::find(config, "physics");
    _rho = toml::find<double>(physics, "rho");
    _nu  = toml::find<double>(physics, "nu");

    const auto& boundary = toml::find(config, "boundary");
    _bc_left      = toml::find<std::string>(boundary, "left");
    _bc_left_dir  = toml::find_or<double>(boundary, "left_val", 0.0);
    _bc_right     = toml::find<std::string>(boundary, "right");
    _bc_right_dir = toml::find_or<double>(boundary, "right_val", 0.0);
    _bc_bottom    = toml::find<std::string>(boundary, "bottom");
    _bc_bottom_dir= toml::find_or<double>(boundary, "bottom_val", 0.0);
    _bc_top       = toml::find<std::string>(boundary, "top");
    _bc_top_dir   = toml::find_or<double>(boundary, "top_val", 0.0);
    _inlet_profile= toml::find_or<std::string>(boundary, "inlet_profile", "Constant");

    const auto& output = toml::find(config, "output");
    _results  = toml::find<std::string>(output, "results");
    _sim_name = toml::find<std::string>(output, "name");

    if (config.contains("obstacle")) {
        const auto& obs = toml::find(config, "obstacle");
        _cyl_cx = toml::find_or<double>(obs, "cx", -100.0);
        _cyl_cy = toml::find_or<double>(obs, "cy", 0.0);
        _cyl_r  = toml::find_or<double>(obs, "radius", 0.0);
    } else {
        _cyl_cx = -100.0; _cyl_cy = 0.0; _cyl_r = 0.0;
    }

    // --- 2. Lecture du bloc CONVERGENCE (CRUCIAL POUR EVITER NaN) ---
    if (config.contains("convergence")) {
        const auto& conv = toml::find(config, "convergence");
        if (conv.contains("schemes")) {
            _conv_schemes = toml::find<std::vector<std::string>>(conv, "schemes");
        } else {
            _conv_schemes = {_scheme}; // Par défaut, celui du bloc [time]
        }
        _has_convergence = true;
        // Lecture du vecteur d'entiers
        _conv_res = toml::find<std::vector<int>>(conv, "resolutions");
        
        // Lecture avec valeurs par défaut de sécurité
        _conv_L   = toml::find_or<double>(conv, "L", 1.0);
        _conv_CFL = toml::find_or<double>(conv, "CFL", 0.5);
        _conv_Nu  = toml::find_or<double>(conv, "nu", _nu); // Utilise _nu physique par défaut
        
        cout << ">> DataFile: Mode Convergence active avec " << _conv_res.size() << " resolutions." << endl;
    } else {
        _has_convergence = false;
    }
}