#include "Function.h"
#include <cmath>

Function::Function(DataFile* data_file) : _df(data_file) {}

// [NOUVEAU] Gestion du profil d'entrée
double Function::GetLeftU_Normal(double y) const {
    double U_val = _df->Get_BC_Left_dir();

    // Si le fichier input demande "Poiseuille"
    if (_df->Get_InletProfile() == "Poiseuille") {
        double y_min = _df->Get_ymin();
        double y_max = _df->Get_ymax();
        double H = y_max - y_min;
        
        // Formule parabolique : 4 * Umax * (y - ymin) * (ymax - y) / H^2
        return 4.0 * U_val * (y - y_min) * (y_max - y) / (H * H);
    }
    
    // Sinon, par défaut : Piston (Constant)
    return U_val;
}

// --- Conditions Initiales U ---
double Function::InitialConditionU(double x, double y) {
    if (_df->Get_BC_Left() == "Dirichlet") {
        // On appelle la fonction ci-dessus pour initialiser cohérent
        // (Parabole ou Constant selon le choix)
        return GetLeftU_Normal(y);
    }
    return 0.0;
}

// --- Conditions Initiales V ---
double Function::InitialConditionV(double x, double y) {
    return 0.0;
}