#include "TFile.h"
#include <TH1.h>
#include <TH2.h>
#include <iostream>

//Macro with previous histograms
//
//Compares and gives the acceptence/efficiency

void acceff(const char* fgen = "/Users/emiliebarreau/alice/TEST_50000evt/Histos.root",
            const char* freco = "/Users/emiliebarreau/alice/TEST_50000evt/Histos_bis.root") 
{
  /*TFile f1(fgen);
  TFile f2(freco);
  f1.ls();
  f2.ls();*/

  TFile* f_gen = TFile::Open(fgen);
  TFile* f_reco = TFile::Open(freco);

  // For Invariant mass
  TH1* hm_gen = static_cast<TH1*>(f_gen->Get("minv"));
  TH1* hm_reco = static_cast<TH1*>(f_reco->Get("minv"));
  double acceff_minv = hm_reco->GetEntries() / hm_gen->GetEntries();
  std::cout << "=====================================================" << std::endl;
  std::cout << "Entries gen et reco : " << hm_gen->GetEntries() << " et " << hm_reco->GetEntries() << std::endl;
  std::cout << "acceff minv : " << acceff_minv * 100 << " \%" << std::endl;
  hm_reco->Divide(hm_gen);
  //std::cout << hm_reco->Integral() << std::endl;

  // For Invariant mass depending of pT
  TH2* hmpt_gen = static_cast<TH2*>(f_gen->Get("minv with pT"));
  TH2* hmpt_reco = static_cast<TH2*>(f_reco->Get("minv with pT"));
  double acceff_minvpt = hmpt_reco->GetEntries() / hmpt_gen->GetEntries();
  std::cout << "=====================================================" << std::endl;
  std::cout << "Entries gen et reco : " << hmpt_gen->GetEntries() << " et " << hmpt_reco->GetEntries() << std::endl;
  std::cout << "acceff pt : " << acceff_minvpt * 100 << " \%" << std::endl;
  hmpt_reco->Divide(hmpt_gen);

  // For Invariant mass depending of y
  TH2* hmy_gen = static_cast<TH2*>(f_gen->Get("minv with y"));
  TH2* hmy_reco = static_cast<TH2*>(f_reco->Get("minv with y"));
  double acceff_minvy = hmy_reco->GetEntries() / hmy_gen->GetEntries();
  std::cout << "=====================================================" << std::endl;
  std::cout << "Entries gen et reco : " << hmy_gen->GetEntries() << " et " << hmy_reco->GetEntries() << std::endl;
  std::cout << "acceff y : " << acceff_minvy * 100 << " \%" << std::endl;
  hmy_reco->Divide(hmy_gen);

  hm_reco->Draw();
  hmpt_reco->Draw("COLZ");
  hmy_reco->Draw("COLZ");
}