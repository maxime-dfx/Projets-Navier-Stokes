// ====================================================================================
//                                  SYSTEM_SOLVER.CPP
// ====================================================================================
// Implémentation de la boucle principale.
// Intègre :
//  - Initialisation des objets (Grille, Laplacien, Schéma Temporel)
//  - Boucle en temps
//  - Sécurité CFL (Arrêt d'urgence si instabilité)
//  - Sorties VTK (avec barre de progression) et Sonde ponctuelle
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
#include <iomanip> // Pour std::setw

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

    // Fréquence de sauvegarde VTK (on vise ~100 frames total pour la vidéo)
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

        // ====================================================================
        // CHECK CFL & NAN (Sécurité)
        // ====================================================================
        // Vérification périodique pour ne pas ralentir le calcul
        if (iter % 10 == 0) {
            double max_U = grid->GetU().cwiseAbs().maxCoeff();
            double max_V = grid->GetV().cwiseAbs().maxCoeff();
            double h_min = std::min(df->Get_hx(), df->Get_hy());
            
            // Calcul de la CFL courante : (V_max * dt) / h
            double current_cfl = std::max(max_U, max_V) * dt / h_min;

            if (current_cfl > 1.5) { // Marge tolérante à 1.5 (Upwind est stable mais faut pas exagérer)
                std::cerr << "\n\n[!!! ALERTE !!!] EXPLOSION CFL DETECTEE" << std::endl;
                std::cerr << "  -> CFL atteinte : " << current_cfl << " (Max recommandé ~1.0)" << std::endl;
                std::cerr << "  -> Vitesse Max  : " << std::max(max_U, max_V) << std::endl;
                std::cerr << "  -> Arrêt d'urgence à t=" << t << std::endl;
                break; // On sort proprement pour sauvegarder
            }
            
            if (std::isnan(max_U) || std::isinf(max_U)) {
                 std::cerr << "\n\n[!!! ERREUR !!!] La simulation contient des NaNs (Not A Number)." << std::endl;
                 break;
            }
        }

        // --- Sortie Sonde (Haute fréquence) ---
        if (probe_file.is_open()) {
            double v_val = grid->GetV()(grid->GetVIndex(i_probe, j_probe));
            probe_file << t << " " << v_val << std::endl;
        }

        // --- Sortie VTK & Feedback (Basse fréquence) ---
        if (write_vtk && vtk && (iter % save_freq == 0)) {
            vtk->Write(iter, t);
            
            // Feedback Console Amélioré (Barre de chargement)
            double progress = t / t_end;
            int barWidth = 40;
            
            std::cout << "\r   [";
            int pos = barWidth * progress;
            for (int i = 0; i < barWidth; ++i) {
                if (i < pos) std::cout << "=";
                else if (i == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0) << " % (t=" 
                      << std::fixed << std::setprecision(3) << t << "s)" << std::flush;
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