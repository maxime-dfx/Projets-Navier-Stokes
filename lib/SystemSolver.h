#ifndef _SYSTEM_SOLVER_H_
#define _SYSTEM_SOLVER_H_

#include "DataFile.h"

class SystemSolver {
public:
    // write_vtk = true  : Mode Simulation (écrit les fichiers, plus lent)
    // write_vtk = false : Mode Convergence (calcul pur, rapide)
    static double Run(DataFile* df, bool write_vtk = true);
};

#endif