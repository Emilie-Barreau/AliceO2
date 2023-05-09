#include "Histogrammer.h"
#include "TH1.h"
#include "TH2.h"
#include <iostream>
#include <TFile.h>
#include <Math/GenVector/VectorUtil.h>

Histogrammer::Histogrammer()
{
  mHistos.emplace_back(new TH1F("pT", "pT", 100, 0., 16.));
  mHistos.emplace_back(new TH1F("Eta", "Eta", 100, -5., -2.));
  mHistos.emplace_back(new TH1F("Rapidity", "Rapidity", 100, -5., -2.));
  mHistos.emplace_back(new TH1F("Minv", "Minv", 80., 2., 5.));

  mHistos2.emplace_back(new TH2F("Minv pT", "M_{inv} depending of p_{T}", 200, 0., 16., 200, 0., 5.));
  mHistos2.emplace_back(new TH2F("Minv rap", "M_{inv} depending of rapidity", 8, -4.5, -2., 200, 0., 5.));
}

void Histogrammer::save(const char* filename)
{
  TFile myFile(filename, "RECREATE");

  myFile.WriteObject(mHistos[0], "pT");
  myFile.WriteObject(mHistos[1], "eta");
  myFile.WriteObject(mHistos[2], "y");
  myFile.WriteObject(mHistos[3], "minv");

  myFile.WriteObject(mHistos2[0], "minv with pT");
  myFile.WriteObject(mHistos2[1], "minv with y");
}

void Histogrammer::fillSingleParticleHistos(const ROOT::Math::PxPyPzMVector& lor)
{
  mHistos[0]->Fill(lor.Pt());
  mHistos[1]->Fill(lor.Eta());
  mHistos[2]->Fill(lor.Rapidity());
}

void Histogrammer::fillDoubleParticleHistos(const ROOT::Math::PxPyPzMVector& lor1, const ROOT::Math::PxPyPzMVector& lor2)
{
  auto lor12 = lor1 + lor2;
  //minv = GetInvariantMass(lor1, lor2);
  mHistos[3]->Fill(lor12.M());
  mHistos2[0]->Fill(lor12.Pt(), lor12.M());
  mHistos2[1]->Fill(lor12.Rapidity(), lor12.M());
}
