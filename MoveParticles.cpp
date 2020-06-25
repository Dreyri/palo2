#include "main.h"

#include <cstring>

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

  // Schleife �ber alle Partikel
  for (int i = 0; i < nr_Particles; i++) {

    // Kraftkomponenten (x,y,z) der Kraft auf aktuellen Partikel (i)
    float Fx = 0, Fy = 0, Fz = 0;

    // Schleife �ber die anderen Partikel die Kraft auf Partikel i aus�ben
    // 5a) Erfordere eine Vektorisierung und gebe an das unser partikel 32 byte
    // aligned ist
#pragma omp simd simdlen(8)
    for (int j = 0; j < nr_Particles; j++) {

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
      Fx += dx * invDrPower32;
      Fy += dy * invDrPower32;
      Fz += dz * invDrPower32;
    }

    // Berechne �nderung der Geschwindigkeit des Partikel i durch einwirkende
    // Kraft
    particles.vx[i] += dt * Fx;
    particles.vy[i] += dt * Fy;
    particles.vz[i] += dt * Fz;
  }

  // Bewege Partikel entsprechend der aktuellen Geschwindigkeit
  for (int i = 0; i < nr_Particles; i++) {
    particles.x[i] += particles.vx[i] * dt;
    particles.y[i] += particles.vy[i] * dt;
    particles.z[i] += particles.vz[i] * dt;
  }
}