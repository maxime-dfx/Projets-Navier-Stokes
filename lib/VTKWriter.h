#ifndef _VTK_WRITER_H_
#define _VTK_WRITER_H_

#include "DataFile.h"
#include "MACgrid.h"
#include <string>

/**
 * @class VTKWriter
 * @brief Classe responsable de l'export des résultats au format VTK (Paraview).
 * Gère la création automatique de l'arborescence de fichiers.
 */
class VTKWriter {
private:
    DataFile* _df;
    MACgrid* _grid;
    std::string _sim_name;

public:
    /**
     * @brief Constructeur
     * @param df Pointeur vers les paramètres de simulation
     * @param grid Pointeur vers le maillage (contient U, V, P)
     * @param sim_name Nom spécifique de la simulation (pour le sous-dossier)
     */
    VTKWriter(DataFile* df, MACgrid* grid, std::string sim_name);

    /**
     * @brief Écrit un fichier .vtk pour l'itération courante.
     * @param iteration Numéro de l'itération (pour le nom de fichier)
     * @param t Temps physique actuel (pour info, pas utilisé par VTK Legacy)
     */
    void Write(int iteration, double t);
};

#endif