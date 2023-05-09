#include "KineReader.h"
#include "Histogrammer.h"
#include "Steer/MCKinematicsReader.h"
#include "MinvTask.h"
#include "SimulationDataFormat/MCTrack.h"
#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <vector>
#include <Math/Vector4D.h>
#include <Math/GenVector/VectorUtil.h>

using namespace std;
using namespace o2::steer;
using namespace o2::mch::eval;

// ===== CONSTRUCTOR =====
KineReader::KineReader() {}

ROOT::Math::PxPyPzMVector getLorentzVector(const o2::MCTrack& t)
{
  return ROOT::Math::PxPyPzMVector(t.Px(), t.Py(), t.Pz(), t.GetMass());
}

// ===== MAIN =====
int main()
{
  //===== INITIALIZATION HISTOS =====
  // unique_ptr<TFile> myFile(TFile::Open("Test.root", "RECREATE"));

  /*TH1F* pt = new TH1F("pt", "Pt", 200, 0., 16.);
  TH1F* eta = new TH1F("eta", "Eta", 100, -5., -2.);
  TH1F* minv = new TH1F("minv", "Minv", 101, 0., 5.);
  TH1F* rap = new TH1F("rap", "Rapidity", 80, -5., -2.);

  TH2F* minv_rap = new TH2F("mrap", "M_{inv} depending of y", 200, -5., -1., 500, 0., 5.);
  TH2F* minv_pt = new TH2F("mpt", "M_{inv} depending of p_{T}", 500, 0., 16., 500, 0., 5.);*/

  // ===== CALLING MCKINEMATICSREADER OBJECT =====
  MCKinematicsReader r("sgn", MCKinematicsReader::Mode::kMCKine); // r is the object

  // ===== TEST =====
  ofstream myfile, myfile2;
  myfile.open("Primary.txt");
  myfile2.open("Second.txt");

  int compteur = 0;
  int compteur2 = 0;

  Histogrammer histogrammer;

  // ===== LOOP FOR MUONS AND J/PSI =====
  for (int evt = 0; evt < r.getNEvents(0); evt++) {
    compteur = 0;
    for (int trk = 0; trk < r.getTracks(evt).size(); trk++) {
      auto t = r.getTrack(evt, trk);
      auto lv = getLorentzVector(*t);
      if (t->GetPdgCode() == 13 || t->GetPdgCode() == -13) {
        histogrammer.fillSingleParticleHistos(lv);
        // pt->Fill(r.getTrack(evt, trk)->GetPt());
        // eta->Fill(r.getTrack(evt, trk)->GetEta());
        compteur += 1;
        compteur2 += 1;
      } else if (t->GetPdgCode() == 443) {
        // rap->Fill(r.getTrack(evt, trk)->GetRapidity());
      } else {
        continue;
      }
      if (t->isPrimary() == 1) {
        myfile << r.getTrack(evt, trk)->GetPdgCode() << "\n";
      } else if (r.getTrack(evt, trk)->isSecondary() == 1) {
        myfile2 << r.getTrack(evt, trk)->GetPdgCode() << "\n";
      }
    }
  }
  myfile.close();
  myfile2.close();

  // ===== LOOP FOR DIMUONS =====
  for (int Evt = 0; Evt < r.getNEvents(0); Evt++) {
    for (int mu1 = 0; mu1 < r.getTracks(Evt).size(); mu1++) {
      auto t1 = r.getTrack(Evt, mu1);
      for (int mu2 = mu1 + 1; mu2 < r.getTracks(Evt).size(); mu2++) {
        auto t2 = r.getTrack(Evt, mu2);
        if (t1->GetPdgCode() == 13 && t2->GetPdgCode() == -13 || t1->GetPdgCode() == -13 && t2->GetPdgCode() == 13) {
          if (t1->isPrimary() == 1 && t2->isPrimary() == 1) {
            auto lv1 = getLorentzVector(*t1);
            auto lv2 = getLorentzVector(*t2);

            /*auto Masses = r.getTrack(Evt, g)->GetMass() * r.getTrack(Evt, g)->GetMass() + r.getTrack(Evt, f)->GetMass() * r.getTrack(Evt, f)->GetMass();
            auto Energies = r.getTrack(Evt, g)->GetEnergy() * r.getTrack(Evt, f)->GetEnergy();
            auto Impulsion = r.getTrack(Evt, g)->Px() * r.getTrack(Evt, f)->Px() + r.getTrack(Evt, g)->Py() * r.getTrack(Evt, f)->Py() + r.getTrack(Evt, g)->Pz() * r.getTrack(Evt, f)->Pz();
            auto Minv = sqrt(Masses + 2 * (Energies - Impulsion));
            auto RAP = -sqrt(r.getTrack(Evt, g)->GetRapidity() * r.getTrack(Evt, g)->GetRapidity() + r.getTrack(Evt, f)->GetRapidity() * r.getTrack(Evt, f)->GetRapidity());
            auto PT = sqrt(r.getTrack(Evt, g)->GetPt() * r.getTrack(Evt, g)->GetPt() + r.getTrack(Evt, f)->GetPt() * r.getTrack(Evt, f)->GetPt());
            ROOT::Math::PxPyPzMVector Lorentz_Gen1(t1->Px(), t1->Py(), t1->Pz(), t1->GetMass());
            ROOT::Math::PxPyPzMVector Lorentz_Gen2(t2->Px(), t2->Py(), t2->Pz(), t2->GetMass());
            auto minv_bis = ROOT::Math::VectorUtil::InvariantMass(Lorentz_Gen1, Lorentz_Gen2);
            auto Lorentz_Gen_Fusion = Lorentz_Gen1 + Lorentz_Gen2;
            auto Lorentz_rap = Lorentz_Gen_Fusion.Rapidity();
            auto Lorentz_pt = Lorentz_Gen_Fusion.Pt();*/
            // minv->Fill(minv_bis);
            // minv_rap->Fill(Lorentz_rap, minv_bis);
            // minv_pt->Fill(Lorentz_pt, minv_bis);

            histogrammer.fillDoubleParticleHistos(lv1, lv2);

          } else {
            continue;
          }
        } else {
          continue;
        }
      }
    }
  }

  std::cout << "Nb muons total : " << compteur2 << std::endl;

  // ===== WRITING HISTOS =====
  /*pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  rap->GetXaxis()->SetTitle("y");
  minv->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");

  minv_rap->GetYaxis()->SetTitle("invariant mass (GeV/c^{2})");
  minv_rap->GetXaxis()->SetTitle("y");
  minv_pt->GetYaxis()->SetTitle("invariant mass (GeV/c^{2})");
  minv_pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  myFile->WriteObject(pt, "pt");
  myFile->WriteObject(eta, "eta");
  myFile->WriteObject(minv, "minv");
  myFile->WriteObject(rap, "rap");
  myFile->WriteObject(minv_rap, "minv_rap");
  myFile->WriteObject(minv_pt, "minv_pt");*/

  histogrammer.save("Histos.root");

  return 0;
}
