#include "Histogrammer.h"
#include "TH1.h"
#include "TH2.h"
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <Math/GenVector/VectorUtil.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>

Histogrammer::Histogrammer()
{
  // initializing bins and limits for histograms
  mHistos.emplace_back(new TH1F("pT", "pT", 100, 0., 16.));
  mHistos.emplace_back(new TH1F("Rapidity", "Rapidity", 100, -4.5, -2.));
  mHistos.emplace_back(new TH1F("Eta", "Eta", 100, -4.5, -2.));
  mHistos.emplace_back(new TH1F("Minv", "Minv", 100, 2., 3.5));

  for (int i = 0; i < 16; i++){
    mMinv.emplace_back(new TH1F(Form("minv_%d",i), Form("minv_%d",i), 100, 2., 3.5));
  }

  mHistos.emplace_back(new TH1F("pT Integral", "pT Integral", 100, 0., 16.));

  mHistos2.emplace_back(new TH2F("Minv pT", "M_{inv} depending of p_{T}", 100, 0., 16., 100, 0., 3.5));
  mHistos2.emplace_back(new TH2F("Minv rap", "M_{inv} depending of rapidity", 100, -4.5, -2., 100, 0., 3.5));
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

  myFile.WriteObject(mHistos[4], "test pt");

  myFile.WriteObject(mHistos2[pt], "minv with pT");
  myFile.WriteObject(mHistos2[y], "minv with y");

  for (int i = 0; i < 16; i++){
    myFile.WriteObject(mMinv[i], "Minv");
  }
}

void Histogrammer::fillSingleParticleHistos(const ROOT::Math::PxPyPzMVector& lor)
{
  // fill 1D
  // mHistos[0]->Fill(lor.Pt());
  mHistos[2]->Fill(lor.Eta());
  // mHistos[2]->Fill(lor.Rapidity());
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

  std::vector<int> ptmin{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  std::vector<int> ptmax{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};

  for (int i=0; i < 16; i++) {
    int bin1 = mHistos2[0]->GetXaxis()->FindBin(ptmin[i]);
    int bin2 = mHistos2[0]->GetXaxis()->FindBin(ptmax[i]);
    mMinv[i] = mHistos2[0]->ProjectionY(Form("projecY_%d",i), bin1, bin2);
    double Integrale = mMinv[i]->Integral();
    mHistos[4]->SetBinContent(bin1, Integrale);
  }
}
