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

template <typename T> struct vec3 {
  T x;
  T y;
  T z;
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
      vec3<float> results[strip_size];

      // 4a) Dieser Teil enthaelt den vektorisierbaren Teil
      for (int strip_index = 0, j = jj; strip_index != strip_size;
           ++strip_index, ++j) {

        // Gravitationsgesetz
        // Berechne Abstand der Partikel i und j
        // 4a) Dieser teil wird nicht vektorisiert im precise fp-modell
        // auch nachdem nur noch diese sehr deutlich vektorisierbare operation
        // vorhanden ist
        results[strip_index].x = partikel[j].x - partikel[i].x;
        results[strip_index].y = partikel[j].y - partikel[i].y;
        results[strip_index].z = partikel[j].z - partikel[i].z;
      }

      // 4a) dieser Teil ist der nicht vektorisierbare Teil wegen der
      // Datenabhaengigkeit von F{x,y,z}
      for (int k = 0; k != strip_size; ++k) {

        // Abschw�chung als zus�tzlicher Abstand, um Singularit�t und
        // Selbst-Interaktion zu vermeiden
        constexpr float softening = 1e-20;

        auto &r = results[k];
        float drSquared = r.x * r.x + r.y * r.y + r.z * r.z + softening;
        float drPower32 = std::sqrt(drSquared) * drSquared;

        float invDrPower32 = 1.f / drPower32;

        // Addiere Kraftkomponenten zur Netto-Kraft
        // 3a) Hier wird die inverse jetzt benutzt um die kosten der Division zu
        // sparen
        Fx += r.x * invDrPower32;
        Fy += r.y * invDrPower32;
        Fz += r.z * invDrPower32;
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