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
  mHistos.emplace_back(new TH1F("pT", "pT", 32, 0., 16.));
  mHistos.emplace_back(new TH1F("Rapidity", "Rapidity", 15, -4., -2.5));
  mHistos.emplace_back(new TH1F("Eta", "Eta", 100, -4.5, -2.));
  mHistos.emplace_back(new TH1F("Minv", "Minv", 175, 2., 3.5));
  mHistos.emplace_back(new TH1F("pT Integral", "pT Integral", 32, 0., 16.));
  mHistos.emplace_back(new TH1F("y Integral", "y Integral", 15, -4., -2.5));

  for (int i = 0; i < 32; i++){
    mMinv.emplace_back(new TH1F(Form("minv_%d",i), Form("minv_%d",i), 175, 2., 3.5));
  }

  for (int i = 0; i < 15; i++){
    mMinv.emplace_back(new TH1F(Form("minvbis_%d",i), Form("minvbis_%d",i), 175, 2., 3.5));
  }


  mHistos2.emplace_back(new TH2F("Minv pT", "M_{inv} depending of p_{T}", 160, 0., 16., 175, 0., 3.5));
  mHistos2.emplace_back(new TH2F("Minv rap", "M_{inv} depending of rapidity", 100, -4.5, -2., 175, 0., 3.5));
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

  myFile.WriteObject(mHistos[4], "pT Integral");
  myFile.WriteObject(mHistos[5], "y Integral");

  myFile.WriteObject(mHistos2[pt], "minv with pT");
  myFile.WriteObject(mHistos2[y], "minv with y");

  for (int i = 0; i < 32; i++){
    myFile.WriteObject(mMinv[i],Form("minv_%d",i));
  }

  for (int i = 0; i < 15; i++){
    myFile.WriteObject(mMinv[i],Form("minvbis_%d",i));
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

  //std::vector<int> ptmin{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  //std::vector<int> ptmax{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
  std::vector<double> ptmin{0., 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5,
                            6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5,
                            11., 11.5, 12., 12.5, 13., 13.5, 14., 14.5, 15., 15.5};
  std::vector<double> ptmax{0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5,
                            6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5,
                            11., 11.5, 12., 12.5, 13., 13.5, 14., 14.5, 15., 15.5, 16.};

  std::vector<double> ymin{-4., -3.9, -3.8, -3.7, -3.6, -3.5, -3.4, -3.3, -3.2, -3.1, -3., -2.9, -2.8, -2.7, -2.6};
  std::vector<double> ymax{-3.9, -3.8, -3.7, -3.6, -3.5, -3.4, -3.3, -3.2, -3.1, -3., -2.9, -2.8, -2.7, -2.6, -2.5};

  for (int i=0; i < 32; i++) {
    int bin1 = mHistos2[0]->GetXaxis()->FindBin(ptmin[i]);
    int bin2 = mHistos2[0]->GetXaxis()->FindBin(ptmax[i]);
    mMinv[i] = mHistos2[0]->ProjectionY(Form("projecY_%d",i), bin1, bin2);
    double Integrale = mMinv[i]->Integral();
    mHistos[4]->SetBinContent(i+1, Integrale);
  }

  for (int i=0; i < 15; i++){
    int bin1 = mHistos2[1]->GetXaxis()->FindBin(ymin[i]);
    int bin2 = mHistos2[1]->GetXaxis()->FindBin(ymax[i]);
    mMinv[i] = mHistos2[1]->ProjectionY(Form("projecYbis_%d",i), bin1, bin2);
    double Integrale = mMinv[i]->Integral();
    mHistos[5]->SetBinContent(i+1, Integrale);
  }
}
