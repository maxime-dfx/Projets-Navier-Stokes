// ====================================================================================
//                                  FUNCTION.H
// ====================================================================================
// Description : Gestionnaire des conditions initiales et aux limites.
//               Centralise la logique "métier" (ex: Profil parabolique).
// ====================================================================================

#ifndef _FUNCTION_H_
#define _FUNCTION_H_

#include "DataFile.h"
#include <string>

class Function {
private:
    DataFile* _df;

public:
    explicit Function(DataFile* df);
    ~Function() = default;

    // --- Conditions Initiales (t=0) ---
    double InitialConditionU(double x, double y);
    double InitialConditionV(double x, double y);

    // --- Helpers pour Conditions aux Limites ---
    // Renvoient true si la paroi est de type Dirichlet (Vitesse imposée)
    bool IsDirichletLeft()   const { return _df->Get_BC_Left()   == "Dirichlet"; }
    bool IsDirichletRight()  const { return _df->Get_BC_Right()  == "Dirichlet"; }
    bool IsDirichletBottom() const { return _df->Get_BC_Bottom() == "Dirichlet"; }
    bool IsDirichletTop()    const { return _df->Get_BC_Top()    == "Dirichlet"; }

    // --- Valeurs aux Bords (Composantes Normales) ---
    // Utilisé pour la pénétration (U gauche/droite, V haut/bas)
    
    // Celui-ci est implémenté dans le .cpp (car il peut être complexe), on garde le nom 'y'
    double GetLeftU_Normal(double y)   const; 

    // Pour ceux-ci, on commente l'argument pour éviter le warning "unused parameter"
    double GetRightU_Normal(double /*y*/)  const { return _df->Get_BC_Right_dir(); }
    double GetBottomV_Normal(double /*x*/) const { return _df->Get_BC_Bottom_dir(); }
    double GetTopV_Normal(double /*x*/)    const { return _df->Get_BC_Top_dir(); }

    // --- Valeurs aux Bords (Composantes Tangentielles) ---
    // Utilisé pour le cisaillement (V gauche/droite, U haut/bas)
    
    double GetLeftV_Tangent(double /*y*/)   const { return 0.0; }
    double GetRightV_Tangent(double /*y*/)  const { return 0.0; }
    double GetBottomU_Tangent(double /*x*/) const { return _df->Get_BC_Bottom_dir(); } // Ex: Cavité entraînée
    double GetTopU_Tangent(double /*x*/)    const { return _df->Get_BC_Top_dir(); }
};

#endif // _FUNCTION_H_