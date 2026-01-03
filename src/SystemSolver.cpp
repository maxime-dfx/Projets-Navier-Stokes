// ====================================================================================
//                                  SYSTEM_SOLVER.CPP
// ====================================================================================
// Implémentation de la boucle principale.
// Ajout de la sonde (Probe) pour le calcul de Strouhal.
// ====================================================================================

#include "SystemSolver.h"
#include "Function.h"
#include "MACgrid.h"
#include "Laplacian.h"
#include "TimeScheme.h"
#include "ErrorAnalysis.h"
#include "VTKWriter.h"

#include <iostream>
#include <fstream>
#include <memory>
#include <algorithm>
#include <cmath>

double SystemSolver::Run(DataFile* df, bool write_vtk) 
{
    // ========================================================================
    // 1. INSTANCIATION (Smart Pointers pour gestion mémoire auto)
    // ========================================================================
    std::cout << ">> [SystemSolver] Initialisation..." << std::endl;
    
    auto fct  = std::make_unique<Function>(df);
    auto grid = std::make_unique<MACgrid>(fct.get(), df);
    // Le Laplacien a besoin des autres composants
    auto lap  = std::make_unique<Laplacian>(fct.get(), df, grid.get());

    // ========================================================================
    // 2. PRÉ-CALCUL
    // ========================================================================
    // Construction de la matrice (coûteux, fait une seule fois)
    lap->BuildMatrix();

    // ========================================================================
    // 3. SÉLECTION DU SCHÉMA TEMPOREL
    // ========================================================================
    std::unique_ptr<TimeScheme> time_scheme;
    std::string s = df->Get_scheme();

    if (s == "RK4" || s == "RungeKutta4") {
        time_scheme = std::make_unique<RungeKutta4Scheme>(df, lap.get(), grid.get());
    } 
    else if (s == "RK2" || s == "RungeKutta2") {
        time_scheme = std::make_unique<RungeKutta2Scheme>(df, lap.get(), grid.get());
    } 
    else {
        time_scheme = std::make_unique<EulerScheme>(df, lap.get(), grid.get());
    }
    std::cout << ">> [SystemSolver] Schema choisi : " << s << std::endl;

    // ========================================================================
    // 4. CONFIGURATION DES SORTIES (SONDE & VTK)
    // ========================================================================
    std::unique_ptr<VTKWriter> vtk = nullptr;
    std::ofstream probe_file;
    int i_probe = -1, j_probe = -1;

    if (write_vtk) {
        // --- VTK Writer ---
        vtk = std::make_unique<VTKWriter>(df, grid.get(), df->Get_SimName());
        
        // --- SONDE (PROBE) pour VON KARMAN ---
        // On place la sonde dans le sillage (Wake)
        double cx = df->Get_CylCx();
        double r  = std::max(df->Get_CylRadius(), 0.01);
        
        // Position cible : 3 diamètres derrière, légèrement décalé en Y (Instabilité)
        double target_x = (cx > 0) ? cx + 3.0 * 2.0 * r : df->Get_xmax() * 0.75;
        double target_y = (cx > 0) ? df->Get_CylCy() + 0.5 * r : df->Get_ymax() * 0.55;

        // Conversion en indices de grille
        j_probe = static_cast<int>((target_x - df->Get_xmin()) / df->Get_hx());
        i_probe = static_cast<int>((target_y - df->Get_ymin()) / df->Get_hy());
        
        // Sécurité bornes
        if (j_probe < 0) j_probe = 0; 
        if (j_probe >= df->Get_Nx()) j_probe = df->Get_Nx() - 1;
        if (i_probe < 0) i_probe = 0; 
        if (i_probe >= df->Get_Ny()) i_probe = df->Get_Ny() - 1;

        // Fichier de sortie texte simple pour Gnuplot/Python
        std::string probe_name = df->Get_results() + "/" + df->Get_SimName() + "_probe.dat";
        probe_file.open(probe_name);
        probe_file << "# Time V_velocity" << std::endl;
        
        std::cout << ">> [SystemSolver] Sonde active en (" << target_x << ", " << target_y 
                  << ") -> Indices [" << i_probe << "," << j_probe << "]" << std::endl;
    }

    // ========================================================================
    // 5. BOUCLE TEMPORELLE
    // ========================================================================
    double t = df->Get_t0();
    double t_end = df->Get_tfinal();
    double dt = df->Get_dt();
    int iter = 0;

    // Fréquence de sauvegarde VTK (on vise ~50-100 frames total pour la vidéo)
    int total_steps = static_cast<int>((t_end - t) / dt);
    int save_freq = std::max(1, total_steps / 100);

    // État initial
    if (write_vtk && vtk) vtk->Write(0, t);

    std::cout << ">> [SystemSolver] Debut de simulation (t_final=" << t_end << ")..." << std::endl;

    while (t < t_end) {
        // --- Avance ---
        time_scheme->Advance();
        t = time_scheme->GetTime();
        iter++;

        // --- Sortie Sonde (Haute fréquence : chaque pas de temps) ---
        if (probe_file.is_open()) {
            double v_val = grid->GetV()(grid->GetVIndex(i_probe, j_probe));
            probe_file << t << " " << v_val << std::endl;
        }

        // --- Sortie VTK (Basse fréquence) ---
        if (write_vtk && vtk && (iter % save_freq == 0)) {
            vtk->Write(iter, t);
            // Petit feedback visuel
            double progress = 100.0 * t / t_end;
            std::cout << "\r   Progress: " << int(progress) << "% (t=" << t << ")" << std::flush;
        }
    }
    std::cout << std::endl << ">> [SystemSolver] Simulation terminee." << std::endl;

    if (probe_file.is_open()) probe_file.close();

    // ========================================================================
    // 6. ANALYSE D'ERREUR (POST-MORTEM)
    // ========================================================================
    // Calculé uniquement si demandé ou si le cas est Poiseuille
    bool check_poiseuille = (df->Get_SimName().find("Poiseuille") != std::string::npos);
    
    if (df->HasConvergenceSection() || check_poiseuille) {
        ErrorAnalysis err(df, grid.get());
        double e = err.ComputePoiseuilleErrorL2();
        std::cout << ">> [SystemSolver] Erreur L2 Finale : " << e << std::endl;
        return e;
    }

    return 0.0;
}