#ifndef HISTOGRAMMER_H
#define HISTOGRAMMER_H

#include "Math/Vector4D.h"
#include <vector>

//Creation of histograms for classes KineReader and MinvTask
//Giving one or two LorentzVector, fill histograms to compare them

class TH1;
class TH2;

class Histogrammer {

public:
Histogrammer();  

void fillSingleParticleHistos(const ROOT::Math::PxPyPzMVector& lor);

void fillDoubleParticleHistos(const ROOT::Math::PxPyPzMVector& lor1, const ROOT::Math::PxPyPzMVector& lor2);

void save(const char* filename);

private:
  std::vector<TH1*> mHistos;
  std::vector<TH2*> mHistos2;
};

#endif