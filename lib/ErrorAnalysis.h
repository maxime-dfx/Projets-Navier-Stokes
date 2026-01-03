#ifndef _ERROR_ANALYSIS_H_
#define _ERROR_ANALYSIS_H_

#include "DataFile.h"
#include "MACgrid.h"
#include <string>

class ErrorAnalysis {
private:
    DataFile* _df;
    MACgrid* _grid;

public:
    ErrorAnalysis(DataFile* df, MACgrid* grid);
    ~ErrorAnalysis();

    // Calcule l'erreur L2 pour le cas Poiseuille (Canal plan)
    // Retourne la norme de l'erreur
    double ComputePoiseuilleErrorL2();
};

#endif