// ====================================================================================
//                                  FUNCTION.CPP
// ====================================================================================
// Implémentation des profils physiques.
// ====================================================================================

#include "Function.h"
#include <cmath>

Function::Function(DataFile* df) : _df(df) {}

// ============================================================================
// PROFIL D'ENTRÉE (U à Gauche)
// ============================================================================
double Function::GetLeftU_Normal(double y) const {
    double U_val = _df->Get_BC_Left_dir();

    // Si le fichier demande un profil de Poiseuille (Parabolique)
    if (_df->Get_InletProfile() == "Poiseuille") {
        double y_min = _df->Get_ymin();
        double y_max = _df->Get_ymax();
        double H = y_max - y_min;
        
        // Formule Poiseuille Plan : U(y) = 4 * Umax * Y * (1-Y)
        // avec Y = (y-ymin)/H coordonnée adimensionnée [0,1]
        double Y = (y - y_min) / H;
        return 4.0 * U_val * Y * (1.0 - Y);
    }
    
    // Sinon profil Constant (Piston)
    return U_val;
}

// ============================================================================
// CONDITIONS INITIALES
// ============================================================================
double Function::InitialConditionU(double x, double y) {
    // Astuce : Si on a Dirichlet à l'entrée, on peut initialiser le champ
    // avec la valeur au bord pour accélérer la convergence (écoulement établi).
    // Ici on reste neutre (0.0).
    return 0.0;
}

double Function::InitialConditionV(double x, double y) {
    // Fluide au repos initialement
    return 0.0;
}