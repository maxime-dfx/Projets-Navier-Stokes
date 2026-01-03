#include "ParticleSystem.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;
using namespace Eigen;

// [CORRECTION] Constructeur avec DataFile
ParticleSystem::ParticleSystem(MACgrid* grid, DataFile* df) : _grid(grid), _df(df) {
}

ParticleSystem::~ParticleSystem() {
    _particles.clear();
}

void ParticleSystem::InitParticles(int n_particles) {
    _particles.clear();
    
    // On adapte la position de départ aux dimensions réelles
    double x_start = _df->Get_xmin() + 0.05; 
    double y_start = _df->Get_ymin() + 0.1;
    double y_end   = _df->Get_ymax() - 0.1;
    
    cout << "--- Initialisation de " << n_particles << " particules ---" << endl;
    
    for(int i = 0; i < n_particles; ++i) {
        Particle p;
        p.x = x_start;
        p.y = y_start + (double)i / (n_particles - 1) * (y_end - y_start);
        p.id = i;
        _particles.push_back(p);
    }
}

// =========================================================================
// INTERPOLATION BILINÉAIRE DYNAMIQUE
// =========================================================================

double ParticleSystem::InterpolateU(double x, double y) {
    // Récupération dynamique des pas d'espace
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();
    double x_min = _df->Get_xmin();
    double y_min = _df->Get_ymin();
    
    // Indexation relative à l'origine
    double j_idx = (x - x_min) / hx;
    double i_idx = (y - y_min) / hy - 0.5; // Décalage U (staggered)

    int j = static_cast<int>(floor(j_idx));
    int i = static_cast<int>(floor(i_idx));
    
    double tx = j_idx - j;
    double ty = i_idx - i;

    // [CORRECTION] Dimensions réelles
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    
    // Clamp sécurisé
    if (i < 0) i = 0; 
    if (i >= Ny - 1) i = Ny - 2;
    
    if (j < 0) j = 0; 
    if (j >= Nx) j = Nx - 1; // U a Nx+1 colonnes, donc j peut aller jusqu'à Nx

    const VectorXd& U = _grid->GetU();
    
    // Sécurité accès mémoire (Assertion fix)
    int idx00 = _grid->GetUIndex(i, j);
    int idx11 = _grid->GetUIndex(i + 1, j + 1);
    if (idx00 < 0 || idx11 >= U.size()) return 0.0; // Sortie silencieuse

    double u00 = U(idx00);
    double u10 = U(_grid->GetUIndex(i, j + 1));
    double u01 = U(_grid->GetUIndex(i + 1, j));
    double u11 = U(idx11);

    double lerp_y1 = (1.0 - tx) * u00 + tx * u10;
    double lerp_y2 = (1.0 - tx) * u01 + tx * u11;
    
    return (1.0 - ty) * lerp_y1 + ty * lerp_y2;
}

double ParticleSystem::InterpolateV(double x, double y) {
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();
    double x_min = _df->Get_xmin();
    double y_min = _df->Get_ymin();

    double j_idx = (x - x_min) / hx - 0.5; // Décalage V
    double i_idx = (y - y_min) / hy;

    int j = static_cast<int>(floor(j_idx));
    int i = static_cast<int>(floor(i_idx));

    double tx = j_idx - j;
    double ty = i_idx - i;
    
    // [CORRECTION] Dimensions réelles
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    
    // Clamp sécurisé
    if (i < 0) i = 0; 
    if (i >= Ny) i = Ny - 1; // V a Ny+1 lignes
    
    if (j < 0) j = 0; 
    if (j >= Nx - 1) j = Nx - 2;

    const VectorXd& V = _grid->GetV();

    int idx00 = _grid->GetVIndex(i, j);
    int idx11 = _grid->GetVIndex(i + 1, j + 1);
    if (idx00 < 0 || idx11 >= V.size()) return 0.0;

    double v00 = V(idx00);
    double v10 = V(_grid->GetVIndex(i, j + 1));
    double v01 = V(_grid->GetVIndex(i + 1, j));
    double v11 = V(idx11);

    double lerp_y1 = (1.0 - tx) * v00 + tx * v10;
    double lerp_y2 = (1.0 - tx) * v01 + tx * v11;

    return (1.0 - ty) * lerp_y1 + ty * lerp_y2;
}

void ParticleSystem::Advance(double dt) 
{
    double x_max = _df->Get_xmax();
    double y_min = _df->Get_ymin();
    double y_max = _df->Get_ymax();
    double x_min = _df->Get_xmin(); // Important pour le reset

    for (auto& p : _particles) {
        double u_part = InterpolateU(p.x, p.y);
        double v_part = InterpolateV(p.x, p.y);
        
        p.x += u_part * dt;
        p.y += v_part * dt;

        // 1. Protection contre les NaNs (Explosion numérique)
        if (std::isnan(p.x) || std::isnan(p.y)) {
            // On réinitialise la particule à l'entrée
            p.x = x_min + 0.05; 
            p.y = y_min + 0.5; 
        }

        // 2. Recyclage normal (Sortie de l'écran)
        if (p.x > x_max) { 
            p.x = x_min + 0.05; 
            p.y = y_min + 0.1 + ((double)rand() / RAND_MAX) * (y_max - y_min - 0.2); 
        }
        
        // 3. Murs glissants
        if (p.y < y_min) p.y = y_min + 0.01;
        if (p.y > y_max) p.y = y_max - 0.01;
    }
}

void ParticleSystem::Save(int iter, const string& path) {
    // On crée un sous-dossier pour ne pas polluer
    string vtk_path = path + "/Particles";
    // Création du dossier si inexistant (C++17)
    // #include <filesystem> est nécessaire en haut si vous utilisez fs::create_directory
    // Sinon, assurez-vous juste que le dossier existe via le main
    
    // Nom du fichier : particles_00100.vtk
    string filename = path + "/particles_" + to_string(iter) + ".vtk";
    ofstream file(filename);
    
    if (file.is_open()) {
        int n_part = _particles.size();

        // En-tête VTK PolyData
        file << "# vtk DataFile Version 3.0" << endl;
        file << "Lagrangian Particles" << endl;
        file << "ASCII" << endl;
        file << "DATASET POLYDATA" << endl;
        file << "POINTS " << n_part << " float" << endl;

        // Écriture des coordonnées (x y z)
        for (const auto& p : _particles) {
            // On ajoute 0.0 pour la coordonnée Z
            file << p.x << " " << p.y << " 0.0" << endl;
        }

        // Optionnel : Ajouter des IDs ou d'autres données
        file << "POINT_DATA " << n_part << endl;
        file << "SCALARS ParticleID int 1" << endl;
        file << "LOOKUP_TABLE default" << endl;
        for (const auto& p : _particles) {
            file << p.id << endl;
        }
    }
    file.close();
}