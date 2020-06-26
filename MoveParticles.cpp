#include "main.h"

#include <cstring>

// 5d) dies sind unsere Abstraktionen um dem assume_aligned und assume directive
#if __has_builtin(__builtin_assume_aligned)
#define ASSUME_ALIGNED(x, alignment)                                           \
  x = static_cast<decltype(x)>(__builtin_assume_aligned(x, alignment))
#elif defined(__INTEL_COMPILER)
#define ASSUME_ALIGNED(x, alignment) __assume_aligned(x, alignment)
#else
#define ASSUME_ALIGNED(x, alignment) (void)0
#endif

#if defined(__INTEL_COMPILER)
#define ASSUME(x) __assume(x)
#elif __has_builtin(__builtin_unreachable)
#define ASSUME(x)                                                              \
  do {                                                                         \
    if (!(x)) {                                                                \
      __builtin_unreachable();                                                 \
    }                                                                          \
  } while (0)
#else
#define ASSUME(x) (void)0
#endif

void MoveParticles(const int nr_Particles, Particle *const partikel,
                   const float dt) {

  // Schleife �ber alle Partikel
  for (int i = 0; i < nr_Particles; i++) {

    // Kraftkomponenten (x,y,z) der Kraft auf aktuellen Partikel (i)
    float Fx = 0, Fy = 0, Fz = 0;

    // Schleife �ber die anderen Partikel die Kraft auf Partikel i aus�ben
    for (int j = 0; j < nr_Particles; j++) {

      // Abschw�chung als zus�tzlicher Abstand, um Singularit�t und
      // Selbst-Interaktion zu vermeiden
      const float softening = 1e-20;

      // Gravitationsgesetz
      // Berechne Abstand der Partikel i und j
      const float dx = partikel[j].x - partikel[i].x;
      const float dy = partikel[j].y - partikel[i].y;
      const float dz = partikel[j].z - partikel[i].z;
      const float drSquared = dx * dx + dy * dy + dz * dz + softening;
      const float drPower32 = pow(drSquared, 3.0 / 2.0);

      // Addiere Kraftkomponenten zur Netto-Kraft
      Fx += dx / drPower32;
      Fy += dy / drPower32;
      Fz += dz / drPower32;
    }

    // Berechne �nderung der Geschwindigkeit des Partikel i durch einwirkende
    // Kraft
    partikel[i].vx += dt * Fx;
    partikel[i].vy += dt * Fy;
    partikel[i].vz += dt * Fz;
  }

  // Bewege Partikel entsprechend der aktuellen Geschwindigkeit
  for (int i = 0; i < nr_Particles; i++) {
    partikel[i].x += partikel[i].vx * dt;
    partikel[i].y += partikel[i].vy * dt;
    partikel[i].z += partikel[i].vz * dt;
  }
}

// 5c) Diese Funktion wurde angepasst um mit der Structure of Arrays
// datenstruktur klarzukommen
void MoveParticlesOpt(const int nr_Particles, ParticleSoA particles,
                      const float dt) {
  // 5d) Hiermit geben wir an das die anzahl der partikel immer Teilbar durch 8
  // ist. Somit garantieren wir das keine Peel loops noetig sind.
  ASSUME(nr_Particles % 8 == 0);

  // 5d) Hiermit geben wir an das unsere pointer alle 32 byte aligned sind
  ASSUME_ALIGNED(particles.x, 32);
  ASSUME_ALIGNED(particles.y, 32);
  ASSUME_ALIGNED(particles.z, 32);
  ASSUME_ALIGNED(particles.vx, 32);
  ASSUME_ALIGNED(particles.vy, 32);
  ASSUME_ALIGNED(particles.vz, 32);

  // 7c) groesse des tiles
  // Werte ueber 16 geben uns Fehler aus wegen einem bug im Intel compiler
  constexpr int tile_size = 16;

// Schleife �ber alle Partikel
// 7c) Hierbei muessen wir die schedule groesse anpassen mit der tile_size
#pragma omp parallel for schedule(static, 512 / tile_size)
  for (int ii = 0; ii < nr_Particles; ii += tile_size) {

    // Kraftkomponenten (x,y,z) der Kraft auf aktuellen Partikel (i)
    // 7c) erstelle arrays in der groesse unserer tile size fuer Fx, Fy und Fz
    alignas(32) float Fx[tile_size] = {};
    alignas(32) float Fy[tile_size] = {};
    alignas(32) float Fz[tile_size] = {};
    // Schleife �ber die anderen Partikel die Kraft auf Partikel i aus�ben
#pragma omp simd simdlen(8)
    for (int j = 0; j < nr_Particles; j++) {
      for (int i = ii; i != (ii + tile_size); i++) {

        // Abschw�chung als zus�tzlicher Abstand, um Singularit�t und
        // Selbst-Interaktion zu vermeiden
        constexpr float softening = 1e-20;

        // Gravitationsgesetz
        // Berechne Abstand der Partikel i und j
        const float dx = particles.x[j] - particles.x[i];
        const float dy = particles.y[j] - particles.y[i];
        const float dz = particles.z[j] - particles.z[i];
        const float drSquared = dx * dx + dy * dy + dz * dz + softening;

        // 3a) Strength reduction, sqrt is gunstiger zu berechnen als pow
        float drPower32 = std::sqrt(drSquared) * drSquared;
        // 3a) einmaliges vorberechnen der inversen
        float invDrPower32 = 1.f / drPower32;

        // Addiere Kraftkomponenten zur Netto-Kraft
        // 3a) Hier wird die inverse jetzt benutzt um die kosten der Division zu
        // sparen
        Fx[i - ii] += dx * invDrPower32;
        Fy[i - ii] += dy * invDrPower32;
        Fz[i - ii] += dz * invDrPower32;
      }
    }

    // Berechne �nderung der Geschwindigkeit des Partikel i durch einwirkende
    // Kraft

    for (int i = ii; i != (ii + tile_size); ++i) {
      particles.vx[i] += dt * Fx[i - ii];
      particles.vy[i] += dt * Fy[i - ii];
      particles.vz[i] += dt * Fz[i - ii];
    }
  }

  // Bewege Partikel entsprechend der aktuellen Geschwindigkeit
  for (int i = 0; i < nr_Particles; i++) {
    particles.x[i] += particles.vx[i] * dt;
    particles.y[i] += particles.vy[i] * dt;
    particles.z[i] += particles.vz[i] * dt;
  }
}