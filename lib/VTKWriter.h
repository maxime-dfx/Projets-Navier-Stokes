// ====================================================================================
//                                  VTKWRITER.H
// ====================================================================================
// Description : Export des résultats au format VTK (Legacy ASCII).
//               Calcule automatiquement la Vorticité pour l'analyse des tourbillons.
// ====================================================================================

#ifndef _VTK_WRITER_H_
#define _VTK_WRITER_H_

#include <string>

// Forward declarations
class DataFile;
class MACgrid;

class VTKWriter {
private:
    DataFile* _df;
    MACgrid* _grid;
    std::string _sim_name;
    std::string _output_dir;

public:
    VTKWriter(DataFile* df, MACgrid* grid, std::string sim_name);
    ~VTKWriter() = default;

    // Écrit le fichier .vtk pour l'itération donnée
    // Champs : Pression (Scalaire), Vitesse (Vecteur), Vorticité (Scalaire)
    void Write(int iteration, double t);
};

#endif // _VTK_WRITER_H_