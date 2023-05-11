#include "Histogrammer.h"
#include "TH1.h"
#include "TH2.h"
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <Math/GenVector/VectorUtil.h>

Histogrammer::Histogrammer()
{
  // initializing bins and limits for histograms
  mHistos.emplace_back(new TH1F("pT", "pT", 100, 0., 10.));
  mHistos.emplace_back(new TH1F("Rapidity", "Rapidity", 100, -5., -2.));
  mHistos.emplace_back(new TH1F("Eta", "Eta", 100, -5., -2.));
  mHistos.emplace_back(new TH1F("Minv", "Minv", 100., 2., 5.));

  mHistos2.emplace_back(new TH2F("Minv pT", "M_{inv} depending of p_{T}", 100, 0., 16., 100, 0., 5.));
  mHistos2.emplace_back(new TH2F("Minv rap", "M_{inv} depending of rapidity", 100, -4.5, -2., 100, 0., 5.));
}

void Histogrammer::save(const char* filename)
{
  // write over the file if already created
  TFile myFile(filename, "RECREATE");

  // new vector indentation to know which histogram is created
  const int pt = 0;
  const int y = 1;
  const int eta = 2;
  const int minv = 3;

  // write the file
  myFile.WriteObject(mHistos[pt], "pT");
  myFile.WriteObject(mHistos[y], "y");
  myFile.WriteObject(mHistos[eta], "eta");
  myFile.WriteObject(mHistos[minv], "minv");

  myFile.WriteObject(mHistos2[pt], "minv with pT");
  myFile.WriteObject(mHistos2[y], "minv with y");
}

void Histogrammer::fillSingleParticleHistos(const ROOT::Math::PxPyPzMVector& lor)
{
  // fill 1D
  //mHistos[0]->Fill(lor.Pt());
  mHistos[2]->Fill(lor.Eta());
  //mHistos[2]->Fill(lor.Rapidity());
}

void Histogrammer::fillDoubleParticleHistos(const ROOT::Math::PxPyPzMVector& lor1, const ROOT::Math::PxPyPzMVector& lor2)
{
  // fill 2D
  auto lor12 = lor1 + lor2;
  mHistos[0]->Fill(lor12.Pt());
  mHistos[1]->Fill(lor12.Rapidity());
  mHistos[3]->Fill(lor12.M());
  mHistos2[0]->Fill(lor12.Pt(), lor12.M());
  mHistos2[1]->Fill(lor12.Rapidity(), lor12.M());
}