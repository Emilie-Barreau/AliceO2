#include "TFile.h"
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <iostream>

//Macro with previous histograms
//
//Compares and gives the acceptence/efficiency

void acceff(const char* fgen = "/Users/emiliebarreau/alice/TEST_50000evt/Histos_gen.root",
            const char* freco = "/Users/emiliebarreau/alice/TEST_50000evt/Histos_reco.root") 
{
  /*TFile f1(fgen);
  TFile f2(freco);
  f1.ls();
  f2.ls();*/

  TFile* f_gen = TFile::Open(fgen);
  TFile* f_reco = TFile::Open(freco);

  TCanvas *c = new TCanvas();
  c->Divide(2,2);

  // For pT
  TH1* hpt_gen = static_cast<TH1*>(f_gen->Get("pT"));
  TH1* hpt_reco = static_cast<TH1*>(f_reco->Get("pT"));
  hpt_gen->Sumw2();
  hpt_reco->Sumw2();
  double acceff_pt = hpt_reco->GetEntries() / hpt_gen->GetEntries();
  std::cout << "=====================================================" << std::endl;
  std::cout << "Entries gen et reco : " << hpt_gen->GetEntries() << " et " << hpt_reco->GetEntries() << std::endl;
  std::cout << "acceff pT : " << acceff_pt * 100 << " \%" << std::endl;
  std::cout << "Test Integral : " << hpt_gen->Integral() << std::endl;
  hpt_reco->Divide(hpt_gen);

  // For y
  TH1* hy_gen = static_cast<TH1*>(f_gen->Get("y"));
  TH1* hy_reco = static_cast<TH1*>(f_reco->Get("y"));
  hy_gen->Sumw2();
  hy_reco->Sumw2();
  double acceff_y = hy_reco->GetEntries() / hy_gen->GetEntries();
  std::cout << "=====================================================" << std::endl;
  std::cout << "Entries gen et reco : " << hy_gen->GetEntries() << " et " << hy_reco->GetEntries() << std::endl;
  std::cout << "acceff y : " << acceff_y * 100 << " \%" << std::endl;
  hy_reco->Divide(hy_gen);

  // For Invariant mass
  TH1* hm_gen = static_cast<TH1*>(f_gen->Get("minv"));
  TH1* hm_reco = static_cast<TH1*>(f_reco->Get("minv"));
  hm_gen->Sumw2();
  hm_reco->Sumw2();
  double acceff_minv = hm_reco->GetEntries() / hm_gen->GetEntries();
  std::cout << "=====================================================" << std::endl;
  std::cout << "Entries gen et reco : " << hm_gen->GetEntries() << " et " << hm_reco->GetEntries() << std::endl;
  std::cout << "acceff minv : " << acceff_minv * 100 << " \%" << std::endl;
  /*double bin1 = hm_reco->FindBin(3.);
  double bin2 = hm_reco->FindBin(3.25);
  double integ = hm_reco->Integral(bin1, bin2);
  std::cout << "Integrale : " << integ << std::endl;*/
  hm_reco->Divide(hm_gen);
  //hm_reco->Scale(1./hm_reco->Integral());

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

  // TEST INTEGRALE
  TH2* htestpt_gen = static_cast<TH2*>(f_gen->Get("pT Integral"));
  TH2* htestpt_reco = static_cast<TH2*>(f_reco->Get("pT Integral"));
  htestpt_gen->Sumw2();
  htestpt_reco->Sumw2();
  double acceff_testpt = htestpt_reco->GetEntries() / htestpt_gen->GetEntries();
  std::cout << "=====================================================" << std::endl;
  std::cout << "Entries gen et reco : " << htestpt_gen->GetEntries() << " et " << htestpt_reco->GetEntries() << std::endl;
  std::cout << "acceff test : " << acceff_testpt * 100 << " \%" << std::endl;
  htestpt_reco->Divide(htestpt_gen);

  TH2* htesty_gen = static_cast<TH2*>(f_gen->Get("y Integral"));
  TH2* htesty_reco = static_cast<TH2*>(f_reco->Get("y Integral"));
  htesty_gen->Sumw2();
  htesty_reco->Sumw2();
  double acceff_testy = htesty_reco->GetEntries() / htesty_gen->GetEntries();
  std::cout << "=====================================================" << std::endl;
  std::cout << "Entries gen et reco : " << htesty_gen->GetEntries() << " et " << htesty_reco->GetEntries() << std::endl;
  std::cout << "acceff test : " << acceff_testy * 100 << " \%" << std::endl;
  htesty_reco->Divide(htesty_gen);

  c->cd(1);
  hpt_reco->SetXTitle("p_{T} (GeV/c^{2})");
  hpt_reco->SetYTitle("A.e");
  hpt_reco->Draw("HIST E");

  c->cd(2);
  hy_reco->SetXTitle("y");
  hy_reco->SetYTitle("A.e");
  hy_reco->Draw("HIST E");

  /*c->cd(3);
  hm_reco->SetXTitle("invariant mass (GeV/c)");
  hm_reco->SetYTitle("A.e");
  hm_reco->Draw("HIST E");*/

  c->cd(3);
  htestpt_reco->SetXTitle("p_{T} (GeV/c^{2})");
  htestpt_reco->SetYTitle("A.e");
  htestpt_reco->Draw("HIST E");

  c->cd(4);
  htesty_reco->SetXTitle("y");
  htesty_reco->SetYTitle("A.e");
  htesty_reco->Draw("HIST E");

  /*c->cd(6);
  hmpt_reco->SetXTitle("p_{T} (GeV/c^{2})");
  hmpt_reco->SetYTitle("invariant mass (GeV/c)");
  hmpt_reco->Draw("COLZ");

  c->cd(7);
  hmy_reco->SetXTitle("y");
  hmy_reco->SetYTitle("invariant mass (GeV/c)");
  hmy_reco->Draw("COLZ");*/
}