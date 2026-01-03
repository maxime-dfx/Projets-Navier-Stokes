// ====================================================================================
//                                ERROR_ANALYSIS.H
// ====================================================================================
// Description : Outil de validation. Calcule l'erreur L2 entre la solution numérique
//               et la solution analytique de Poiseuille.
// ====================================================================================

#ifndef _ERROR_ANALYSIS_H_
#define _ERROR_ANALYSIS_H_

#include "DataFile.h"
#include "MACgrid.h"

class ErrorAnalysis {
private:
    DataFile* _df;
    MACgrid* _grid;

public:
    ErrorAnalysis(DataFile* df, MACgrid* grid);
    ~ErrorAnalysis() = default;

    // Calcule la norme L2 de l'erreur sur la vitesse U
    // E = sqrt( sum( (U_num - U_exact)^2 * dV ) )
    double ComputePoiseuilleErrorL2();
};

#endif // _ERROR_ANALYSIS_H_