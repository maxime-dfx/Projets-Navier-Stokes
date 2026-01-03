// ====================================================================================
//                                  MACGRID.H
// ====================================================================================
// Description : Gestion de la grille décalée (Staggered Grid) type Marker-and-Cell.
//               Stockage des champs U, V, P et gestion des obstacles (Solid Mask).
// ====================================================================================

#ifndef _MACGRID_H_
#define _MACGRID_H_

#include <Eigen/Dense>
#include <vector>

// Forward declarations
class Function;
class DataFile;

class MACgrid {
private:
    Function* _fct;
    DataFile* _df;
    
    // Dimensions
    int _Nx, _Ny;

    // Champs Physiques (1D Vector mapping 2D)
    Eigen::VectorXd _p; // Pression aux centres
    Eigen::VectorXd _U; // Vitesse X aux faces verticales
    Eigen::VectorXd _V; // Vitesse Y aux faces horizontales
    
    // Masques Solides (Pour les obstacles arbitraires)
    // uint8_t est plus efficace que bool ou int pour les tableaux de masques
    std::vector<uint8_t> _is_solid_u; 
    std::vector<uint8_t> _is_solid_v; 

public:
    // ============================================================================
    // CONSTRUCTEUR
    // ============================================================================
    MACgrid(Function* fct, DataFile* df);
    
    // ============================================================================
    // ACCESSEURS (GETTERS / SETTERS)
    // ============================================================================
    Function* GetFunction() const { return _fct; }
    
    const Eigen::VectorXd& GetP() const { return _p; }
    const Eigen::VectorXd& GetU() const { return _U; }
    const Eigen::VectorXd& GetV() const { return _V; }

    // Utilisation de références const pour éviter la copie
    void SetP(const Eigen::VectorXd& p) { _p = p; }
    void SetU(const Eigen::VectorXd& u) { _U = u; }
    void SetV(const Eigen::VectorXd& v) { _V = v; }
    
    // ============================================================================
    // INDEXATION (MAPPING 2D -> 1D)
    // ============================================================================
    // P est de taille Nx * Ny
    int GetPIndex(int i, int j) const { return j + i * _Nx; }
    // U est de taille (Nx+1) * Ny
    int GetUIndex(int i, int j) const { return j + i * (_Nx + 1); }
    // V est de taille Nx * (Ny+1)
    int GetVIndex(int i, int j) const { return j + i * _Nx; }

    // ============================================================================
    // COORDONNÉES SPATIALES
    // ============================================================================
    Eigen::VectorXd GetPcoord(int i, int j) const;
    Eigen::VectorXd GetUcoord(int i, int j) const;
    Eigen::VectorXd GetVcoord(int i, int j) const;

    // ============================================================================
    // GESTION DES OBSTACLES
    // ============================================================================
    // Génère les masques solides basés sur DataFile
    void BuildObstacles();
    
    // Vérification rapide de l'état solide
    bool IsSolidU(int i, int j) const { return _is_solid_u[GetUIndex(i,j)]; }
    bool IsSolidV(int i, int j) const { return _is_solid_v[GetVIndex(i,j)]; }
};

#endif // _MACGRID_H_