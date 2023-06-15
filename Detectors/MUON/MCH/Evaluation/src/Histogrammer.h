#ifndef HISTOGRAMMER_H
#define HISTOGRAMMER_H

#include "Math/Vector4D.h"
#include "MCHEvaluation/ExtendedTrack.h"
#include <vector>

// Creation of histograms for classes KineReader and MinvTask
//
// Giving one or two LorentzVector, fill histograms to compare them

class TH1;
class TH2;
class ExtendedTrack;

class Histogrammer
{

 public:
  Histogrammer();                                                                                              // constructor

  void fillSingleParticleHistos(const ROOT::Math::PxPyPzMVector& lor);                                         // filling histos 1D

  void fillDoubleParticleHistos(const ROOT::Math::PxPyPzMVector& lor1, const ROOT::Math::PxPyPzMVector& lor2); // filling histos 2D

  void fillTest(const o2::mch::eval::ExtendedTrack& E1, const o2::mch::eval::ExtendedTrack& E2);

  void DEtest(const int& DE);

  void DSclust(const int& DS_Index);

  void save(const char* filename);                                                                             // saving and writing in a root file

 private:
  std::vector<TH1*> mHistos;
  std::vector<TH2*> mHistos2;
  std::vector<TH1*> mMinv;
};

#endif