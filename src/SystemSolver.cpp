#include "SystemSolver.h"
#include "Function.h"
#include "MACgrid.h"
#include "Laplacian.h"
#include "TimeScheme.h"
#include "ErrorAnalysis.h"
#include "VTKWriter.h"

#include <iostream>
#include <memory>
#include <algorithm> // std::max

using namespace std;

double SystemSolver::Run(DataFile* df, bool write_vtk) {
    // 1. Instanciation
    auto fct = make_unique<Function>(df);
    auto grid = make_unique<MACgrid>(fct.get(), df);
    auto lap = make_unique<Laplacian>(fct.get(), df, grid.get()); 

    // 2. Initialisation
    lap->BuildMatrix();

    // 3. Schéma Temporel
    unique_ptr<TimeScheme> time_scheme = nullptr;
    string scheme_name = df->Get_scheme();

    if (scheme_name == "RungeKutta4" || scheme_name == "RK4") {
        time_scheme = make_unique<RungeKutta4Scheme>(df, lap.get(), grid.get());
    } 
    else if (scheme_name == "RungeKutta2" || scheme_name == "RK2") {
        time_scheme = make_unique<RungeKutta2Scheme>(df, lap.get(), grid.get());
    } 
    else {
        time_scheme = make_unique<EulerScheme>(df, lap.get(), grid.get());
    }
    // 4. Sorties (Conditionnelles)
    unique_ptr<VTKWriter> vtk_writer = nullptr;
    if (write_vtk) {
        vtk_writer = make_unique<VTKWriter>(df, grid.get(), df->Get_SimName());
    }
    
    // 5. Boucle Temporelle
    double t = df->Get_t0();
    double t_final = df->Get_tfinal();
    double dt = df->Get_dt();
    int iteration = 0;
    int total_iters = static_cast<int>((t_final - t) / dt);
    int save_freq = std::max(1, total_iters / 20); 

    // Sauvegarde initiale
    if (write_vtk && vtk_writer) vtk_writer->Write(0, t);

    while (t < t_final) {
        time_scheme->Advance();
        t = time_scheme->GetTime();
        iteration++;

        if (write_vtk && vtk_writer && (iteration % save_freq == 0)) {
           vtk_writer->Write(iteration, t); 
        }
    }

    // 6. Calcul d'Erreur (si applicable)
    // On ne calcule l'erreur que si on est en mode convergence ou test Poiseuille
    double error = 0.0;
    if (df->HasConvergenceSection() || df->Get_SimName().find("Poiseuille") != string::npos) {
        ErrorAnalysis error_ana(df, grid.get());
        error = error_ana.ComputePoiseuilleErrorL2();
    }

    return error;
}