#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <iomanip>
#include <chrono> // Pour le temps d'exécution

#include "DataFile.h"
#include "SystemSolver.h"

using namespace std;

int main(int argc, char** argv) {
    // 1. Vérification arguments
    if (argc < 2) {
        cerr << "Usage: ./ns_solver_2d input/file.toml" << endl;
        return 1;
    }

    // 2. Chargement Configuration
    auto df = make_unique<DataFile>(argv[1]);

    // 3. Détection du Mode
    if (df->HasConvergenceSection()) {
        cout << "=== ETUDE COMPARATIVE (Euler / RK2 / RK4) ===" << endl;

        const auto& Resolutions = df->Get_Conv_Resolutions();
        const auto& Schemes = df->Get_Conv_Schemes();
        double L = df->Get_Conv_L();
        double CFL = df->Get_Conv_CFL();
        double Nu = df->Get_Conv_Nu();
        string base_name = df->Get_SimName();

        // 1. Boucle sur les SCHEMAS
        for (const string& scheme : Schemes) {
            df->Set_Scheme(scheme); // Changement dynamique
            cout << "\n>>> Test Schema : " << scheme << " <<<" << endl;

            string out_file = "results/convergence_" + scheme + ".dat";
            ofstream file(out_file);
            file << "# h  ErreurL2  dt  Time(s)" << endl;

            // 2. Boucle sur les RESOLUTIONS
            for (int N : Resolutions) {
                double h = L / N;
                double dt = CFL * (h * h) / Nu;

                df->Set_Nx(N); df->Set_Ny(N);
                df->Set_hx(h); df->Set_hy(h);
                df->Set_dt(dt);

                // Mode silencieux (false) pour aller vite
                auto start = chrono::high_resolution_clock::now();
                double error = SystemSolver::Run(df.get(), false); 
                auto end = chrono::high_resolution_clock::now();

                chrono::duration<double> diff = end - start;
                cout << "   N=" << setw(3) << N << " | Err=" << scientific << error << endl;
                file << h << " " << error << " " << dt << " " << diff.count() << endl;
            }
        }
    } 
    else {
        // --- MODE SIMULATION STANDARD ---
        cout << "\n=== LANCEMENT SIMULATION STANDARD ===" << endl;
        
        auto start = chrono::high_resolution_clock::now();
        
        // Lancement (write_vtk = true pour visualiser)
        double err = SystemSolver::Run(df.get(), true);
        
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> diff = end - start;
        
        cout << "Simulation terminee en " << diff.count() << "s." << endl;
        if (err > 0.0) cout << "Erreur L2 finale : " << err << endl;
    }

    return 0;
}