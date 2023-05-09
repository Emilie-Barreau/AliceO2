#include "KineReader.h"
#include "Histogrammer.h"
#include "Steer/MCKinematicsReader.h"
#include "MinvTask.h"
#include "SimulationDataFormat/MCTrack.h"
#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <vector>
#include <Math/Vector4D.h>
#include <RooCrystalBall.h>
#include <Math/GenVector/VectorUtil.h>

using namespace std;
using namespace o2::steer;
using namespace o2::mch::eval;

class TH1;
class TH2;
class TF1;
class KineReader;

// ATTENTION ! getTracks et getTrack pas pareil !!!
// Par event (5000), il y a plusieurs traces d'ou le getTrackS
// Avec getTrackS, preciser l'event et on a des infos sur l'ensemble des traces
// Avec getTrack, 1 trace à chaque fois donc preciser l'event et l'id de la trace
// Avec getNEvents, preciser la source (ici par defaut y'en a qu'une)

// ===== CONSTRUCTOR =====
o2::mch::eval::KineReader::KineReader() {}

// ===== TEST CLASSES =====
void createHistos(){
  cout << "truc" << endl;
}

ROOT::Math::PxPyPzMVector getLorentzVector(const o2::MCTrack& t) {
  return ROOT::Math::PxPyPzMVector(t.Px(), t.Py(), t.Pz(), t.GetMass());
}


// ===== MAIN =====
int main()
{
  //===== INITIALIZATION HISTOS =====
  createHistos();
  unique_ptr<TFile> myFile(TFile::Open("Test.root", "RECREATE"));

// histogrammer.save("Test.root");

  TH1F* phi = new TH1F("phi", "Phi", 100, 0., 10.);
  TH1F* pt = new TH1F("pt", "Pt", 200, 0., 16.);
  TH1F* eta = new TH1F("eta", "Eta", 100, -5., -2.);
  TH1F* mass = new TH1F("mass", "Mass", 40, 0., 1.);
  TH1F* minv = new TH1F("minv", "Minv", 101, 0., 5.);
  TH1F* rap = new TH1F("rap", "Rapidity", 80, -5., -2.);
  TH1F* pt_bis = new TH1F("pt", "Pt", 200, 0., 20.);

  TH2F* minv_rap = new TH2F("mrap", "M_{inv} depending of y", 200, -5., -1., 500, 0., 5.);
  TH2F* minv_pt = new TH2F("mpt", "M_{inv} depending of p_{T}", 500, 0., 16., 500, 0., 5.);

  // ===== CALLING MCKINEMATICREADER OBJECT =====
  MCKinematicsReader r("sgn", MCKinematicsReader::Mode::kMCKine); // r is the object
  std::cout << "bonjour le monde " << r.isInitialized() << " " << r.getTracks(0).size() << " " << r.getTrack(0, 1)->GetPdgCode() << std::endl;
  std::cout << "Nb d'events : " << r.getNEvents(0) << std::endl;

  // ===== TEST =====
  ofstream myfile, myfile2;
  myfile.open("Primary.txt");
  myfile2.open("Second.txt");

  ofstream file_Nevent, file_getTracks;
  file_Nevent.open("Events.txt");
  file_Nevent << r.getNEvents(0) << endl;
  file_Nevent.close();

  int Event = r.getNEvents(0);
  o2::mch::eval::KineReader kr;
  kr.Nevents = Event;

  file_getTracks.open("GetTracks.txt");

  int compteur = 0;
  int compteur2 = 0;

  Histogrammer histogrammer;

  // ===== LOOP FOR MUONS AND J/PSI =====
  for (int evt = 0; evt < r.getNEvents(0); evt++) {        
    compteur = 0;
    for (int trk = 0; trk < r.getTracks(evt).size(); trk++) { 
      // cout << "j i : " << j << " " << i << endl;
      auto t = r.getTrack(evt,trk);
      auto lv = getLorentzVector(*t);

      if (t->GetPdgCode() == 13 || t->GetPdgCode() == -13) { 
        histogrammer.fillSingleParticleHisto(lv);
        phi->Fill(r.getTrack(evt, trk)->GetPhi());
        pt->Fill(r.getTrack(evt, trk)->GetPt());
        // rap->Fill(r.getTrack(j, i)->GetRapidity());
        eta->Fill(r.getTrack(evt, trk)->GetEta());
        mass->Fill(r.getTrack(evt, trk)->GetMass());
        compteur += 1;
        compteur2 += 1;
        // cout << "Nb of tracks muons : " << r.getTracks(j).size() << endl;
      } else if (r.getTrack(evt, trk)->GetPdgCode() == 443) {                                  
        if (r.getTrack(evt, trk)->GetRapidity() > -4 && r.getTrack(evt, trk)->GetRapidity() < -2, 5) { 
          rap->Fill(r.getTrack(evt, trk)->GetRapidity());
          pt_bis->Fill(r.getTrack(evt, trk)->GetPt());
        } else {
          continue;
        }
      } else {
        continue;
      }
      if (r.getTrack(evt, trk)->isPrimary() == 1) { // counting of primary and secondary muons
        myfile << r.getTrack(evt, trk)->GetPdgCode() << "\n";
      } else if (r.getTrack(evt, trk)->isSecondary() == 1) {
        myfile2 << r.getTrack(evt, trk)->GetPdgCode() << "\n";
      }
    }
    // cout << "Nb of muons per track : " << compteur << endl;
    // cout << "Nb of tracks : " << r.getTracks(j).size() << endl;
  }
  myfile.close();
  myfile2.close();

  // ===== LOOP FOR DIMUONS =====
  for (int Evt = 0; Evt < r.getNEvents(0); Evt++) {                                                                                                                              
    file_getTracks << r.getTracks(Evt).size() << "\n";
    for (int mu1 = 0; mu1 < r.getTracks(Evt).size(); mu1++) {                                                                                                                       
      for (int mu2 = mu1 + 1; mu2 < r.getTracks(Evt).size(); mu2++) {                                                                                                                   
        if (r.getTrack(Evt, mu1)->GetPdgCode() == 13 && r.getTrack(Evt, mu2)->GetPdgCode() == -13 || r.getTrack(Evt, mu1)->GetPdgCode() == -13 && r.getTrack(Evt, mu2)->GetPdgCode() == 13) { 
          if (r.getTrack(Evt, mu1)->isPrimary() == 1 && r.getTrack(Evt, mu2)->isPrimary() == 1) {                                                

            /*auto Masses = r.getTrack(Evt, g)->GetMass() * r.getTrack(Evt, g)->GetMass() + r.getTrack(Evt, f)->GetMass() * r.getTrack(Evt, f)->GetMass();
            auto Energies = r.getTrack(Evt, g)->GetEnergy() * r.getTrack(Evt, f)->GetEnergy();
            auto Impulsion = r.getTrack(Evt, g)->Px() * r.getTrack(Evt, f)->Px() + r.getTrack(Evt, g)->Py() * r.getTrack(Evt, f)->Py() + r.getTrack(Evt, g)->Pz() * r.getTrack(Evt, f)->Pz();
            auto Minv = sqrt(Masses + 2 * (Energies - Impulsion));
            auto RAP = -sqrt(r.getTrack(Evt, g)->GetRapidity() * r.getTrack(Evt, g)->GetRapidity() + r.getTrack(Evt, f)->GetRapidity() * r.getTrack(Evt, f)->GetRapidity());
            auto PT = sqrt(r.getTrack(Evt, g)->GetPt() * r.getTrack(Evt, g)->GetPt() + r.getTrack(Evt, f)->GetPt() * r.getTrack(Evt, f)->GetPt());*/

            /*TLorentzVector P;
            auto track = *(r.getTrack(h, g));
            track.Get4Momentum(P);*/
            // cout << "Mother ID : " << r.getTrack(h, g)->getMotherTrackId() << endl;

            // test avec Lorentz vectors
            //KineReader lor1, lor2;
            //auto Vec1 = lor1.Lorentz_Gen1(r.getTrack(h, g)->Px(), r.getTrack(h, g)->Py(), r.getTrack(h, g)->Pz(), r.getTrack(h, g)->GetMass());
            //auto Vec2 = lor2.Lorentz_Gen2(r.getTrack(h, f)->Px(), r.getTrack(h, f)->Py(), r.getTrack(h, f)->Pz(), r.getTrack(h, f)->GetMass());

            ROOT::Math::PxPyPzMVector Lorentz_Gen1(r.getTrack(Evt, mu1)->Px(), r.getTrack(Evt, mu1)->Py(), r.getTrack(Evt, mu1)->Pz(), r.getTrack(Evt, mu1)->GetMass());
            ROOT::Math::PxPyPzMVector Lorentz_Gen2(r.getTrack(Evt, mu2)->Px(), r.getTrack(Evt, mu2)->Py(), r.getTrack(Evt, mu2)->Pz(), r.getTrack(Evt, mu2)->GetMass());
            auto minv_bis = ROOT::Math::VectorUtil::InvariantMass(Lorentz_Gen1, Lorentz_Gen2);
            auto Lorentz_Gen_Fusion = Lorentz_Gen1 + Lorentz_Gen2;
            auto Lorentz_rap = Lorentz_Gen_Fusion.Rapidity();
            auto Lorentz_pt = Lorentz_Gen_Fusion.Pt();
            minv->Fill(minv_bis);
            // auto rap_bis = -sqrt(Lorentz_Gen1.Rapidity() * Lorentz_Gen1.Rapidity() + Lorentz_Gen2.Rapidity() * Lorentz_Gen2.Rapidity());
            minv_rap->Fill(Lorentz_rap, minv_bis);
            // auto pt_bis = sqrt(Lorentz_Gen1.Pt() * Lorentz_Gen1.Pt() + Lorentz_Gen2.Pt() * Lorentz_Gen2.Pt());
            minv_pt->Fill(Lorentz_pt, minv_bis);

          } else {
            continue;
          }
        } else {
          continue;
        }
      }
    }
  }

  file_getTracks.close();

  std::cout << "Nb muons total : " << compteur2 << std::endl;

  // ===== WRITING HISTOS =====
  pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  rap->GetXaxis()->SetTitle("y");
  mass->GetXaxis()->SetTitle("mass (GeV/c^{2})");
  minv->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");

  minv_rap->GetYaxis()->SetTitle("invariant mass (GeV/c^{2})");
  minv_rap->GetXaxis()->SetTitle("y");
  minv_pt->GetYaxis()->SetTitle("invariant mass (GeV/c^{2})");
  minv_pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  myFile->WriteObject(phi, "phi");
  myFile->WriteObject(pt, "pt");
  myFile->WriteObject(eta, "eta");
  myFile->WriteObject(mass, "mass");
  myFile->WriteObject(minv, "minv");
  myFile->WriteObject(rap, "rap");
  myFile->WriteObject(pt_bis, "pt_bis");
  myFile->WriteObject(minv_rap, "minv_rap");
  myFile->WriteObject(minv_pt, "minv_pt");

  return 0;
}

