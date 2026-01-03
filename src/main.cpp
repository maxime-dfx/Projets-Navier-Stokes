#include <iostream>
#include <string>
#include <memory>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <filesystem> 

#include "DataFile.h"
#include "SystemSolver.h"

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    // 1. Vérification des arguments (Juste le nom du fichier suffit maintenant)
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input_file.toml" << std::endl;
        return EXIT_FAILURE;
    }

    std::string toml_path = argv[1];

    try {
        // 2. Chargement de la configuration
        auto df = std::make_unique<DataFile>(toml_path);

        // 3. AIGUILLAGE AUTOMATIQUE
        // Si le fichier contient une section [convergence], on lance l'étude.
        if (df->HasConvergenceSection()) {
            std::cout << "=========================================" << std::endl;
            std::cout << "      ETUDE DE CONVERGENCE (AUTO)        " << std::endl;
            std::cout << "      Detectee via le fichier .toml      " << std::endl;
            std::cout << "=========================================" << std::endl;

            // Création du dossier de sortie
            std::string res_dir = df->Get_results();
            if (!fs::exists(res_dir)) fs::create_directories(res_dir);

            const auto& resolutions = df->Get_Conv_Resolutions();
            const auto& schemes = df->Get_Conv_Schemes();
            double L = df->Get_Conv_L();
            double CFL = df->Get_Conv_CFL();
            double Nu = df->Get_Conv_Nu();

            for (const auto& scheme : schemes) {
                df->Set_Scheme(scheme);
                std::cout << "\n>>> Test Schema : " << scheme << " <<<" << std::endl;

                std::string out_file = res_dir + "/convergence_" + scheme + ".dat";
                std::ofstream file(out_file);
                
                if (!file.is_open()) {
                    std::cerr << "ERREUR : Impossible de créer " << out_file << std::endl;
                    continue;
                }

                file << "# h  ErreurL2  dt  Time(s)" << std::endl;

                for (int N : resolutions) {
                    double h = L / static_cast<double>(N);
                    double dt = CFL * (h * h) / Nu;

                    df->Set_Nx(N); df->Set_Ny(N);
                    df->Set_hx(h); df->Set_hy(h);
                    df->Set_dt(dt);

                    auto start = std::chrono::high_resolution_clock::now();
                    
                    // Lancement silencieux (false)
                    double error = SystemSolver::Run(df.get(), false);
                    
                    auto end = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> diff = end - start;

                    std::cout << "   N=" << std::setw(3) << N 
                              << " | h=" << std::scientific << std::setprecision(2) << h
                              << " | Err=" << error 
                              << " | CPU=" << std::fixed << std::setprecision(4) << diff.count() << "s" << std::endl;

                    file << h << " " << error << " " << dt << " " << diff.count() << std::endl;
                }
                file.close();
            }
            std::cout << "\n>> Fin de l'etude. Données sauvegradées dans " << res_dir << std::endl;
        } 
        // 4. Mode SIMULATION STANDARD (Pas de section convergence)
        else {
            std::cout << "=========================================" << std::endl;
            std::cout << "      SIMULATION STANDARD (VTK)          " << std::endl;
            std::cout << "=========================================" << std::endl;
            
            auto start = std::chrono::high_resolution_clock::now();
            double err = SystemSolver::Run(df.get(), true); // true = écrit les VTK
            auto end = std::chrono::high_resolution_clock::now();
            
            std::chrono::duration<double> diff = end - start;
            std::cout << ">> Temps Execution Total : " << diff.count() << " s" << std::endl;
            if (err > 0.0) std::cout << ">> Erreur Finale (L2) : " << err << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "ERREUR FATALE : " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}