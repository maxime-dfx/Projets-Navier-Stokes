// ====================================================================================
//                                  SYSTEM_SOLVER.H
// ====================================================================================
// Description : Chef d'orchestre de la simulation.
//               Initialise les objets, lance la boucle temporelle,
//               gère les sorties (VTK) et la sonde (Probe).
// ====================================================================================

#ifndef _SYSTEM_SOLVER_H_
#define _SYSTEM_SOLVER_H_

#include "DataFile.h"

class SystemSolver {
public:
    // Méthode statique principale.
    // write_vtk : Si true, génère les fichiers lourds .vtk pour Paraview.
    //             Si false (mode convergence), ne sort rien ou juste les erreurs.
    static double Run(DataFile* df, bool write_vtk = true);
};

#endif // _SYSTEM_SOLVER_H_