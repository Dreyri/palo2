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

struct vector_results {
  float dx;
  float dy;
  float dz;
  float invDrPower32;
};

void MoveParticlesOpt(const int nr_Particles, Particle *const partikel,
                      const float dt) {

  constexpr int strip_size = 8;

  // Schleife �ber alle Partikel
  for (int i = 0; i < nr_Particles; i++) {

    // Kraftkomponenten (x,y,z) der Kraft auf aktuellen Partikel (i)
    float Fx = 0, Fy = 0, Fz = 0;

    // Schleife �ber die anderen Partikel die Kraft auf Partikel i aus�ben
    for (int jj = 0; jj < nr_Particles; jj += strip_size) {
      vector_results results[strip_size];

      // 4a) Dieser Teil enthaelt den vektorisierbaren Teil
      for (int strip_index = 0; strip_index != strip_size; ++strip_index) {
        // 4a) Berechnung des wirklichen j Wert weil unser jetziges j um
        // strip_size springt
        int j = jj + strip_index;

        // Abschw�chung als zus�tzlicher Abstand, um Singularit�t und
        // Selbst-Interaktion zu vermeiden
        constexpr float softening = 1e-20;

        // Gravitationsgesetz
        // Berechne Abstand der Partikel i und j
        float dx_ = partikel[j].x - partikel[i].x;
        float dy_ = partikel[j].y - partikel[i].y;
        float dz_ = partikel[j].z - partikel[i].z;
        const float drSquared = dx_ * dx_ + dy_ * dy_ + dz_ * dz_ + softening;

        // 3a) Strength reduction, sqrt is gunstiger zu berechnen als pow
        float drPower32 = std::sqrt(drSquared) * drSquared;
        // 3a) einmaliges vorberechnen der inversen
        float invDrPower32_ = 1.f / drPower32;

        results[strip_index].dx = dx_;
        results[strip_index].dy = dy_;
        results[strip_index].dz = dz_;
        results[strip_index].invDrPower32 = invDrPower32_;
      }

      // 4a) dieser Teil ist der nicht vektorisierbare Teil wegen der
      // Datenabhaengigkeit von F{x,y,z}
      for (int k = 0; k != strip_size; ++k) {
        // Addiere Kraftkomponenten zur Netto-Kraft
        // 3a) Hier wird die inverse jetzt benutzt um die kosten der Division zu
        // sparen
        auto &r = results[k];
        float invDrPower32_ = r.invDrPower32;

        Fx += r.dx * invDrPower32_;
        Fy += r.dy * invDrPower32_;
        Fz += r.dz * invDrPower32_;
      }
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