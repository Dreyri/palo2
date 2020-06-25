#include "main.h"

#include <mm_malloc.h>

/// 2a) Diese Funktion nimmt die startwerte, die ergebnisse der
/// Referenzimplementation und unserer optimierten Ergebnisse.
/// Es wird getested ob der Unterschied zwischen den Startzustand und unseren
/// neuen Zustand innerhalb eines 0.1% Fehler zu den alten Ergebnissen haelt.
bool test_result(Particle *partikel_start, Particle *partikel,
                 Particle *optimized, int count,
                 const float thresh = 0.001f) noexcept {

  auto equal_epsilon = [](float start, float expected, float optimized,
                          float thresh) {
    float delta_expected = expected - start;
    float delta_optimized = optimized - start;
    float results_delta = std::fabs(delta_optimized - delta_expected);
    bool acceptable_error =
        results_delta < std::max(thresh, std::fabs(thresh * delta_expected));
    return acceptable_error;
  };

  const char *member;
  float expect;
  float found;

  int i = 0;

  for (; i != count; ++i) {
    if (!equal_epsilon(partikel_start[i].x, partikel[i].x, optimized[i].x,
                       thresh)) {
      member = "x";
      expect = partikel[i].x;
      found = optimized[i].x;

      goto fail;
    }

    if (!equal_epsilon(partikel_start[i].y, partikel[i].y, optimized[i].y,
                       thresh)) {
      member = "y";
      expect = partikel[i].y;
      found = optimized[i].y;

      goto fail;
    }
    if (!equal_epsilon(partikel_start[i].z, partikel[i].z, optimized[i].z,
                       thresh)) {
      member = "z";
      expect = partikel[i].z;
      found = optimized[i].z;

      goto fail;
    }
    if (!equal_epsilon(partikel_start[i].vx, partikel[i].vx, optimized[i].vx,
                       thresh)) {
      member = "vx";
      expect = partikel[i].vx;
      found = optimized[i].vx;

      goto fail;
    }
    if (!equal_epsilon(partikel_start[i].vy, partikel[i].vy, optimized[i].vy,
                       thresh)) {
      member = "vy";
      expect = partikel[i].vy;
      found = optimized[i].vy;

      goto fail;
    }
    if (!equal_epsilon(partikel_start[i].vz, partikel[i].vz, optimized[i].vz,
                       thresh)) {
      member = "vz";
      expect = partikel[i].vz;
      found = optimized[i].vz;

      goto fail;
    }
  }

  return true;

fail:
  printf("Error in field %s at index %i\n", member, i);
  printf("expected: %f, found: %f\n", expect, found);
  return false;
}

void initParticles(Particle *const partikel, const int nr_Particles) {
  srand(0);
  for (int i = 0; i < nr_Particles; i++) {
    partikel[i].x = float(rand()) / RAND_MAX;
    partikel[i].y = float(rand()) / RAND_MAX;
    partikel[i].z = float(rand()) / RAND_MAX;
    partikel[i].vx = float(rand()) / RAND_MAX;
    partikel[i].vy = float(rand()) / RAND_MAX;
    partikel[i].vz = float(rand()) / RAND_MAX;
  }
}

void copyParticles(Particle *const partikel_src, Particle *const partikel_dst,
                   const int nr_Particles) {
  for (int i = 0; i < nr_Particles; i++) {
    partikel_dst[i].x = partikel_src[i].x;
    partikel_dst[i].y = partikel_src[i].y;
    partikel_dst[i].z = partikel_src[i].z;
    partikel_dst[i].vx = partikel_src[i].vx;
    partikel_dst[i].vy = partikel_src[i].vy;
    partikel_dst[i].vz = partikel_src[i].vz;
  }
}

int main() {
  // Problemgr��e und Anzahl und Gr��e der Zeitschritte definieren
  constexpr int nrOfParticles = 16384;
  constexpr int nrRuns =
      10; // Anzahl der L�ufe und der Zeitschritte der Simulation
  constexpr int skipRuns =
      3; // Anzahl der Messungen, die nicht in Mittelwert ber�cksichtigt werden
  constexpr float dt = 0.01f; // L�nge eines Zeitschrittes

  // 5a) Alloziere die Particle mit korrektem alignment
  Particle *partikel_start =
      static_cast<Particle *>(_mm_malloc(sizeof(Particle) * nrOfParticles, 32));
  Particle *partikel =
      static_cast<Particle *>(_mm_malloc(sizeof(Particle) * nrOfParticles, 32));
  // Particle *partikel_start = new Particle[nrOfParticles];
  // Particle *partikel = new Particle[nrOfParticles];
  copyParticles(partikel_start, partikel, nrOfParticles);

  // Initiaslisierung der Partikel mit Zufallswerten
  initParticles(partikel_start, nrOfParticles);

  // Messen der Performance
  double runtimeStep[nrRuns] = {0.}; // Sammlung der Laufzeiten der Steps
  double GFlopsStep[nrRuns] = {0.};  // Sammlung der Leistungen der Steps
  double meanRuntime = 0.;
  double stdRuntime = 0.;
  double meanGFlops = 0.;
  double stdGFlops = 0.;

  // Berechnung der Anzahl an GFLOPs der Berechnung
  const float NrOfGFLOPs =
      20.0 * 1e-9 * float(nrOfParticles) * float(nrOfParticles - 1);
  printf("#### Runtime Measurements Particle Simulation  ###\n");

  for (int run = 0; run < nrRuns; run++) {
    copyParticles(partikel_start, partikel, nrOfParticles);

    const double tStart = omp_get_wtime(); // Start der Zeitmessung
    MoveParticlesOpt(nrOfParticles, partikel,
                     dt);                // Funktion, die optimiert werden soll
    const double tEnd = omp_get_wtime(); // Ende der Zeitmessung

    if (run == 0) {
      Particle *reference_particles = new Particle[nrOfParticles];
      copyParticles(partikel_start, reference_particles, nrOfParticles);
      MoveParticles(nrOfParticles, reference_particles, dt);

      if (!test_result(partikel_start, reference_particles, partikel,
                       nrOfParticles, 0.001f)) {
        std::terminate();
      } else {
        std::printf("Test passed\n");
      }
    }

    runtimeStep[run] = tEnd - tStart;
    GFlopsStep[run] = NrOfGFLOPs / runtimeStep[run];
    if (run >= skipRuns) { // Berechnung Mittelwerte
      meanRuntime += runtimeStep[run];
      meanGFlops += GFlopsStep[run];
    }

    printf("Run %d: Runtime: %f03,\t GFLOPS/s: %f01, \t %s\n", run,
           runtimeStep[run], GFlopsStep[run],
           (run < skipRuns ? "Not in Average" : ""));
    fflush(stdout); // Ausgabebuffer leeren
  }
  // Berechnung der Mittelwerte
  double nrRunsInStatistics = (double)(nrRuns - skipRuns);
  meanRuntime /= nrRunsInStatistics;
  meanGFlops /= nrRunsInStatistics;

  // Berechnung der Mittelwertfehler
  for (int i = skipRuns; i < nrRuns; i++) {
    stdRuntime +=
        (runtimeStep[i] - meanRuntime) * (runtimeStep[i] - meanRuntime);
    stdGFlops += (GFlopsStep[i] - meanGFlops) * (GFlopsStep[i] - meanGFlops);
  }
  stdRuntime =
      sqrt(stdRuntime / (nrRunsInStatistics * (nrRunsInStatistics - 1)));
  stdGFlops = sqrt(stdGFlops / (nrRunsInStatistics * (nrRunsInStatistics - 1)));

  // Ausgabe der Ergebnisse
  printf("\n\n####### Average Performance #########\n");
  printf("Average Runtime: %f03 +- %f Seconds \n", meanRuntime, stdRuntime);
  printf("Average Performance: %f03 +- %f03 GFLOPS/s \n", meanGFlops,
         stdGFlops);
  printf("#####################################\n");

  _mm_free(partikel);
  _mm_free(partikel_start);
}