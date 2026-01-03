#ifndef _FUNCTION_H_
#define _FUNCTION_H_

#include "DataFile.h" 
#include <string>
#include <cmath>

class Function {
private:
    DataFile* _df;

public:
    Function(DataFile* data_file); 

    // Conditions Initiales (t=0)
    double InitialConditionU(double x, double y);
    double InitialConditionV(double x, double y);

    // --- LOGIQUE CONDITIONNELLE ---
    bool IsDirichletLeft()   const { return _df->Get_BC_Left() == "Dirichlet"; } 
    bool IsDirichletRight()  const { return _df->Get_BC_Right() == "Dirichlet"; }
    bool IsDirichletBottom() const { return _df->Get_BC_Bottom() == "Dirichlet"; }
    bool IsDirichletTop()    const { return _df->Get_BC_Top() == "Dirichlet"; } 

    // --- VITESSES NORMALES (Traversée du mur) ---
    // [MODIF] Déclaration seule ici, implémentation dans .cpp
    double GetLeftU_Normal(double y) const;
    
    double GetRightU_Normal(double y)  const { return 0.0; } 
    double GetBottomV_Normal(double x) const { return 0.0; } 
    double GetTopV_Normal(double x)    const { return 0.0; } 

    // --- VITESSES TANGENTIELLES (Glissement du mur) ---
    double GetLeftV_Tangent(double y)   const { return 0.0; }
    double GetRightV_Tangent(double y)  const { return 0.0; } 
    double GetBottomU_Tangent(double x) const { return _df->Get_BC_Bottom_dir(); } 
    double GetTopU_Tangent(double x)    const { return _df->Get_BC_Top_dir(); }    
};

#endif