#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <omp.h>

// Struct zur Beschreibung eines Teilchens
struct Particle {
  float x, y, z;    // Koordinaten des Teilchens
  float vx, vy, vz; // Geschwindigkeiten des Teilchens
};

// 5c)
// Alle x, y, z, vx, vy, vz haben ihr eigenes array. Hierdurch koennen unit
// strided access erzeugt werden.
struct ParticleSoA {
  float *__restrict__ x;
  float *__restrict__ y;
  float *__restrict__ z;
  float *__restrict__ vx;
  float *__restrict__ vy;
  float *__restrict__ vz;
};

// Prototypen
void MoveParticles(const int nr_Particles, Particle *const partikel,
                   const float dt);
void MoveParticlesOpt(const int nr_Particles, ParticleSoA particles,
                      const float dt);