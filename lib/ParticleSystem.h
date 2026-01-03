#ifndef _PARTICLE_SYSTEM_H_
#define _PARTICLE_SYSTEM_H_

#include "MACgrid.h"
#include "DataFile.h" // [AJOUT]
#include <vector>
#include <string>

struct Particle {
    double x, y;
    int id;
};

class ParticleSystem {
private:
    MACgrid* _grid;
    DataFile* _df; // [AJOUT]
    std::vector<Particle> _particles;
    
    double InterpolateU(double x, double y);
    double InterpolateV(double x, double y);

public:
    // On passe DataFile au constructeur
    ParticleSystem(MACgrid* grid, DataFile* df); 
    ~ParticleSystem();

    void InitParticles(int n_particles);
    void Advance(double dt);
    void Save(int iter, const std::string& path);
};

#endif