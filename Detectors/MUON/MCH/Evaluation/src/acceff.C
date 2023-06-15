#include "TFile.h"
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <iostream>
#include <TLegend.h>

// Macro with previous histograms
//
// Compares and gives the acceptence/efficiency

void acceff(const char* fgen = "/Users/emiliebarreau/alice/TEST_50000evt/Histos_gen.root",
            const char* freco = "/Users/emiliebarreau/alice/TEST_50000evt/Histos_reco.root",
            const char* freco_cut = "/Users/emiliebarreau/alice/TEST_50000evt/Histos_reco_cut.root")
{
  /*TFile f1(fgen);
  TFile f2(freco);
  f1.ls();
  f2.ls();*/

  TFile* f_gen = TFile::Open(fgen);
  TFile* f_reco = TFile::Open(freco);
  TFile* f_reco_cut = TFile::Open(freco_cut);

  TCanvas* c = new TCanvas();
  c->Divide(2, 1);
  c->SetWindowSize(1500, 700);

  // For pT
  TH1* hpt_gen = static_cast<TH1*>(f_gen->Get("pT"));
  TH1* hpt_reco = static_cast<TH1*>(f_reco->Get("pT"));
  hpt_gen->Sumw2();
  hpt_reco->Sumw2();
  double acceff_pt = hpt_reco->GetEntries() / hpt_gen->GetEntries();
  std::cout << "acceff pT : " << acceff_pt * 100 << " \%" << std::endl;
  hpt_reco->Divide(hpt_gen);

  TH1* hpt_reco_cut = static_cast<TH1*>(f_reco_cut->Get("pT"));
  hpt_reco_cut->Sumw2();
  double acceff_pt_cut = hpt_reco_cut->GetEntries() / hpt_gen->GetEntries();
  std::cout << "acceff pT cut : " << acceff_pt_cut * 100 << " \%" << std::endl;
  hpt_reco_cut->Divide(hpt_gen);

  // For y
  TH1* hy_gen = static_cast<TH1*>(f_gen->Get("y"));
  TH1* hy_reco = static_cast<TH1*>(f_reco->Get("y"));
  hy_gen->Sumw2();
  hy_reco->Sumw2();
  double acceff_y = hy_reco->GetEntries() / hy_gen->GetEntries();
  std::cout << "acceff y : " << acceff_y * 100 << " \%" << std::endl;
  hy_reco->Divide(hy_gen);

  TH1* hy_reco_cut = static_cast<TH1*>(f_reco_cut->Get("y"));
  hy_reco_cut->Sumw2();
  double acceff_y_cut = hy_reco_cut->GetEntries() / hy_gen->GetEntries();
  std::cout << "acceff y : " << acceff_y_cut * 100 << " \%" << std::endl;
  hy_reco_cut->Divide(hy_gen);

  // For Invariant mass
  /*TH1* hm_gen = static_cast<TH1*>(f_gen->Get("minv"));
  TH1* hm_reco = static_cast<TH1*>(f_reco->Get("minv"));
  hm_gen->Sumw2();
  hm_reco->Sumw2();
  hm_reco->Divide(hm_gen);

  // For Invariant mass depending of pT
  TH2* hmpt_gen = static_cast<TH2*>(f_gen->Get("minv with pT"));
  TH2* hmpt_reco = static_cast<TH2*>(f_reco->Get("minv with pT"));
  hmpt_reco->Divide(hmpt_gen);

  // For Invariant mass depending of y
  TH2* hmy_gen = static_cast<TH2*>(f_gen->Get("minv with y"));
  TH2* hmy_reco = static_cast<TH2*>(f_reco->Get("minv with y"));
  hmy_reco->Divide(hmy_gen);

  // TEST INTEGRALE
  TH2* htestpt_gen = static_cast<TH2*>(f_gen->Get("pT Integral"));
  TH2* htestpt_reco = static_cast<TH2*>(f_reco->Get("pT Integral"));
  htestpt_gen->Sumw2();
  htestpt_reco->Sumw2();
  double acceff_testpt = htestpt_reco->GetEntries() / htestpt_gen->GetEntries();
  htestpt_reco->Divide(htestpt_gen);

  TH2* htesty_gen = static_cast<TH2*>(f_gen->Get("y Integral"));
  TH2* htesty_reco = static_cast<TH2*>(f_reco->Get("y Integral"));
  htesty_gen->Sumw2();
  htesty_reco->Sumw2();
  double acceff_testy = htesty_reco->GetEntries() / htesty_gen->GetEntries();
  htesty_reco->Divide(htesty_gen);*/

  /*c->cd(1);
  hpt_reco->SetXTitle("p_{T} (GeV/c^{2})");
  hpt_reco->SetYTitle("A.e");
  hpt_reco->Draw("HIST E");
  hpt_reco_cut->Draw("SAME HIST E");
  hpt_reco_cut->SetLineColor(2);
  hpt_reco->GetXaxis()->SetLabelSize(0.04);
  hpt_reco->GetXaxis()->SetTitleSize(0.04);
  hpt_reco->GetYaxis()->SetLabelSize(0.04);
  hpt_reco->GetYaxis()->SetTitleSize(0.04);
  auto legend1 = new TLegend(0.1, 0.9, 0.6, 0.8); // position x, position y, width x, width y
  legend1->SetHeader(Form("A.e difference : %.3f (percent)", (acceff_pt - acceff_pt_cut) * 100));
  legend1->AddEntry(hpt_reco, "Perfect detector");
  legend1->AddEntry(hpt_reco_cut, "Ch 1, 3, 7 removed");
  legend1->Draw();

  c->cd(2);
  hy_reco->SetXTitle("y");
  hy_reco->SetYTitle("A.e");
  hy_reco->Draw("HIST E");
  hy_reco_cut->Draw("SAME HIST E");
  hy_reco_cut->SetLineColor(2);
  hy_reco->GetXaxis()->SetLabelSize(0.04);
  hy_reco->GetXaxis()->SetTitleSize(0.04);
  hy_reco->GetYaxis()->SetLabelSize(0.04);
  hy_reco->GetYaxis()->SetTitleSize(0.04);
  auto legend2 = new TLegend(0.2, 0.9, 0.7, 0.8); // position x1, position y1, position x2, position y2 ?
  legend2->SetHeader(Form("A.e difference : %.3f (percent)", (acceff_y - acceff_y_cut) * 100));
  legend2->AddEntry(hy_reco, "Perfect detector");
  legend2->AddEntry(hy_reco_cut, "Ch 1, 3, 7 removed");
  legend2->Draw();*/

  double Ae_ratio_y = (acceff_y_cut / acceff_y) * 100;
  double Ae_ratio_pt = (acceff_pt_cut / acceff_pt) * 100;

  if (Ae_ratio_pt > 75. && Ae_ratio_y > 75.){
    hpt_reco_cut->SetLineColor(4);
    hy_reco_cut->SetLineColor(4);
  } else {
    hpt_reco_cut->SetLineColor(2);
    hy_reco_cut->SetLineColor(2);
  }

  c->cd(1);
  hpt_reco_cut->SetXTitle("p_{T} (GeV/c^{2})");
  hpt_reco_cut->SetYTitle("A.#epsilon ratio");
  hpt_reco_cut->SetAxisRange(0., 2., "Y");
  hpt_reco_cut->GetXaxis()->SetLabelSize(0.04);
  hpt_reco_cut->GetXaxis()->SetTitleSize(0.04);
  hpt_reco_cut->GetYaxis()->SetLabelSize(0.04);
  hpt_reco_cut->GetYaxis()->SetTitleSize(0.04);
  hpt_reco_cut->Divide(hpt_reco);
  hpt_reco_cut->Draw("HIST E");
  auto legend3 = new TLegend(0.1, 0.9, 0.6, 0.8);
  //hpt_reco_cut->SetStats(0); \
  legend3->SetHeader(Form("A.e ratio : %.3f (percent)", Ae_ratio_pt));
  legend3->AddEntry(hpt_reco_cut, "Ratio Ch 5 removed");
  legend3->Draw();

  c->cd(2);
  hy_reco_cut->SetXTitle("y");
  hy_reco_cut->SetYTitle("A.#epsilon ratio");
  hpt_reco_cut->SetAxisRange(0., 2., "Y");
  hy_reco_cut->GetXaxis()->SetLabelSize(0.04);
  hy_reco_cut->GetXaxis()->SetTitleSize(0.04);
  hy_reco_cut->GetYaxis()->SetLabelSize(0.04);
  hy_reco_cut->GetYaxis()->SetTitleSize(0.04);
  hy_reco_cut->Divide(hy_reco);
  hy_reco_cut->Draw("HIST E");
  auto legend4 = new TLegend(0.1, 0.9, 0.6, 0.8);
  legend4->SetHeader(Form("A.e ratio : %.3f (percent)", Ae_ratio_y));
  legend4->AddEntry(hy_reco_cut, "Ratio Ch 5 removed");
  legend4->Draw();

  /*c->cd(3);
  htestpt_reco->SetXTitle("p_{T} (GeV/c^{2})");
  htestpt_reco->SetYTitle("A.e");
  htestpt_reco->Draw("HIST E");*/

  /*c->cd(4);
  htesty_reco->SetXTitle("y");
  htesty_reco->SetYTitle("A.e");
  htesty_reco->Draw("HIST E");*/

  /*c->cd(6);
  hmpt_reco->SetXTitle("p_{T} (GeV/c^{2})");
  hmpt_reco->SetYTitle("invariant mass (GeV/c)");
  hmpt_reco->Draw("COLZ");*/

  /*c->cd(7);
  hmy_reco->SetXTitle("y");
  hmy_reco->SetYTitle("invariant mass (GeV/c)");
  hmy_reco->Draw("COLZ");*/
}