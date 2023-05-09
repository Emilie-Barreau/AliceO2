#ifndef HISTOGRAMMER_H
#define HISTOGRAMMER_H

#include "Math/Vector4D.h"
#include <vector>

class TH1;

class Histogrammer {

public:

Histogrammer();

void fillSingleParticleHisto(const ROOT::Math::PxPyPzMVector& p);

void save(const char* filename);

private:
  std::vector<TH1*> mHistos;
};

#endif