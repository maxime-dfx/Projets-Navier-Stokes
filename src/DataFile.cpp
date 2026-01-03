// ====================================================================================
//                                  DATAFILE.CPP
// ====================================================================================
// Parsing robuste via toml11.
// Correction : Stockage par valeur (auto) au lieu de référence (auto&) pour éviter
// les dangling references qui corrompent les données d'initialisation.
// ====================================================================================

#include "DataFile.h"
#include "toml.hpp" // Librairie header-only standard
#include <iostream>
#include <filesystem>

DataFile::DataFile(std::string file_name) 
    : _file_name(file_name), 
      _has_convergence_section(false)
{
    // Lecture du fichier TOML
    auto config = toml::parse(file_name);

    // --- 1. MAILLAGE ---
    // CORRECTION ICI : "const auto" au lieu de "const auto&"
    const auto space = toml::find(config, "space");
    
    _xmin = toml::find<double>(space, "xmin");
    _xmax = toml::find<double>(space, "xmax");
    _ymin = toml::find<double>(space, "ymin");
    _ymax = toml::find<double>(space, "ymax");
    _hx   = toml::find<double>(space, "hx");
    _hy   = toml::find<double>(space, "hy");
    
    // Calcul dérivé
    _Nx = static_cast<int>((_xmax - _xmin) / _hx);
    _Ny = static_cast<int>((_ymax - _ymin) / _hy);

    // --- 2. TEMPS ---
    const auto time = toml::find(config, "time"); // Pas de &
    _t0     = toml::find<double>(time, "t0");
    _tfinal = toml::find<double>(time, "tfinal");
    _dt     = toml::find<double>(time, "dt");
    _scheme = toml::find_or<std::string>(time, "scheme", "Euler");

    // --- 3. PHYSIQUE ---
    const auto phys = toml::find(config, "physics"); // Pas de &
    _rho = toml::find<double>(phys, "rho");
    _nu  = toml::find<double>(phys, "nu");

    // --- 4. CONDITIONS AUX LIMITES ---
    const auto bc = toml::find(config, "boundary"); // Pas de &
    _bc_left   = toml::find<std::string>(bc, "left");
    _bc_right  = toml::find<std::string>(bc, "right");
    _bc_bottom = toml::find<std::string>(bc, "bottom");
    _bc_top    = toml::find<std::string>(bc, "top");

    // Valeurs associées (ex: Vitesse d'entrée)
    _bc_left_val   = toml::find_or<double>(bc, "left_val", 0.0);
    _bc_right_val  = toml::find_or<double>(bc, "right_val", 0.0);
    _bc_bottom_val = toml::find_or<double>(bc, "bottom_val", 0.0);
    _bc_top_val    = toml::find_or<double>(bc, "top_val", 0.0);

    // Profil (Constant vs Poiseuille)
    _inlet_profile = toml::find_or<std::string>(bc, "inlet_profile", "Constant");

    // --- 5. OBSTACLE (Optionnel) ---
    if (config.contains("obstacle")) {
        const auto obs = toml::find(config, "obstacle"); // Pas de &
        _cyl_cx = toml::find_or<double>(obs, "cx", -100.0); 
        _cyl_cy = toml::find_or<double>(obs, "cy", 0.0);
        _cyl_r  = toml::find_or<double>(obs, "radius", 0.0);
    } else {
        _cyl_cx = -100.0; _cyl_cy = 0.0; _cyl_r = 0.0;
    }

    // --- 6. SORTIES ---
    const auto out = toml::find(config, "output"); // Pas de &
    _results_dir = toml::find<std::string>(out, "results");
    _sim_name    = toml::find<std::string>(out, "name");

    // --- 7. MODE CONVERGENCE (Optionnel) ---
    if (config.contains("convergence")) {
        _has_convergence_section = true;
        const auto conv = toml::find(config, "convergence"); // Pas de &
        
        _conv_resolutions = toml::find<std::vector<int>>(conv, "resolutions");
        _conv_schemes     = toml::find_or<std::vector<std::string>>(conv, "schemes", {_scheme});
        _conv_L           = toml::find_or<double>(conv, "L", 1.0);
        _conv_CFL         = toml::find_or<double>(conv, "CFL", 0.5);
        _conv_Nu          = toml::find_or<double>(conv, "nu", _nu);

        std::cout << ">> [DataFile] Mode Convergence detecte (" 
                  << _conv_resolutions.size() << " resolutions)." << std::endl;
    }
}